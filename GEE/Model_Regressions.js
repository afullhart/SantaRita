var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid');

var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area(1).divide(10000)); 
                     });

// Safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

var f_date = '2019-05-26', l_date = '2019-05-31';

var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)                     // FILTER LOCATION FIRST
  .filterDate(f_date, l_date);                // FILTER DATE SECOND

// Extract projection from the first image in the collection
var projSent2 = sent2_ic.first().select('B2').projection();

// Mosaic, clip, select bands, and APPLY SCALE FACTOR (0.0001)
var sent2_im = sent2_ic
  .mosaic()
  .clip(v_extent)
  .select(['B2', 'B3', 'B4', 'B5', 'B8'])
  .multiply(0.0001) // Ensures prediction image matches your scaled training data
  .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

// Calculate NDVI using scaled true reflectance
var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');

// Calculate MCARI using scaled true reflectance
var mcari = sent2_im.expression(
    '((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
      'B3': sent2_im.select('B3'), // Green
      'B4': sent2_im.select('B4'), // Red
      'B5': sent2_im.select('B5')  // Red Edge 1
}).rename('MCARI');

// Add both NDVI and MCARI bands to the image
sent2_im = sent2_im.addBands([ndvi, mcari]);

var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI'];


// =========================================================================
// K-FOLDS CROSS VALIDATION (K=5)
// =========================================================================

var k_folds = 5;
var v_list_seeds = ee.List([123, 456, 789, 111, 333]);
var fold_list = ee.List.sequence(0, 4);

// Add a random column and create fold IDs (0, 1, 2, 3, 4)
var fc_folds = fc.randomColumn('random', 123).map(function(ft) {
  return ft.set('fold', ee.Number(ft.get('random')).multiply(k_folds).floor());
});

// Reusable function to run CV for any given target property
function runCV(targetProp) {
  var fold_rmse = fold_list.map(function(fold) {
    var i = ee.Number(fold);
    var current_seed = v_list_seeds.get(i);
    
    // Split into 80% train (out-of-fold) and 20% test (in-fold)
    var train_fc = fc_folds.filter(ee.Filter.neq('fold', i));
    var test_fc = fc_folds.filter(ee.Filter.eq('fold', i));
    
    // Train Model (bagFraction included here)
    var rf = ee.Classifier.smileRandomForest({
      numberOfTrees: 500, 
      minLeafPopulation: 5, 
      bagFraction: 0.8, 
      seed: current_seed
    })
    .setOutputMode('REGRESSION')
    .train({
      features: train_fc, 
      classProperty: targetProp, 
      inputProperties: inputProps
    });
    
    // Classify test set
    var tested = test_fc.classify({
      classifier: rf, 
      outputName: 'predicted'
    });
    
    // Calculate RMSE for this fold
    var mse = tested.map(function(ft) {
      var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
      return ft.set('sq_diff', diff.multiply(diff));
    }).reduceColumns({
      reducer: ee.Reducer.mean(),
      selectors: ['sq_diff']
    }).get('mean');
    
    return ee.Number(mse).sqrt();
  });
  
  // Return the median RMSE across all 5 folds
  return fold_rmse.reduce(ee.Reducer.median());
}

print('Median RMSE (BGR):', runCV('BGR'));
print('Median RMSE (LPI):', runCV('LPI'));
print('Median RMSE (MFT):', runCV('MFT'));


// =========================================================================
// FINAL MODEL TRAINING & PREDICTION
// =========================================================================

// Train Model 1: BGR (bagFraction removed)
var model_bgr = ee.Classifier.smileRandomForest({numberOfTrees: 500, minLeafPopulation: 5, seed: 123})
  .setOutputMode('REGRESSION')
  .train({features: fc, classProperty: 'BGR', inputProperties: inputProps});

// Train Model 2: LPI (bagFraction removed)
var model_lpi = ee.Classifier.smileRandomForest({numberOfTrees: 500, minLeafPopulation: 5, seed: 123})
  .setOutputMode('REGRESSION')
  .train({features: fc, classProperty: 'LPI', inputProperties: inputProps});

// Train Model 3: MFT (bagFraction removed)
var model_mft = ee.Classifier.smileRandomForest({numberOfTrees: 500, minLeafPopulation: 5, seed: 123})
  .setOutputMode('REGRESSION')
  .train({features: fc, classProperty: 'MFT', inputProperties: inputProps});

// Classify the image for each variable and explicitly rename the output bands
var pred_bgr = sent2_im.classify(model_bgr).rename('Pred_BGR');
var pred_lpi = sent2_im.classify(model_lpi).rename('Pred_LPI');
var pred_mft = sent2_im.classify(model_mft).rename('Pred_MFT');

// Combine the three prediction bands into a single image
var combined_preds = ee.Image([pred_bgr, pred_lpi, pred_mft]);



// =========================================================================
// SAMPLING & CSV EXPORT
// =========================================================================

// Convert the grid features to point geometries located at their centers
var fc_centers = fc.map(function(ft) {
  return ft.centroid(1); // 1m max error for centroid calculation
});

// Sample the combined prediction image at these exact center points.
// We also retain the true 'BGR', 'LPI', and 'MFT' values from the training points.
var sampled_data = combined_preds.sampleRegions({
  collection: fc_centers,
  properties: ['BGR', 'LPI', 'MFT'], // Retain the True values
  scale: projSent2.nominalScale(), // Matches the 10m Sentinel-2 scale
  tileScale: 4
});

// Format the collection for a clean CSV export
var export_csv = sampled_data.map(function(ft) {
  // Return a Feature with null geometry for a clean tabular CSV
  return ee.Feature(null, { 
    'True_BGR': ft.get('BGR'),
    'Predicted_BGR': ft.get('Pred_BGR'),
    'True_LPI': ft.get('LPI'),
    'Predicted_LPI': ft.get('Pred_LPI'),
    'True_MFT': ft.get('MFT'),
    'Predicted_MFT': ft.get('Pred_MFT')
  });
});

// Export the table to Google Drive as a CSV
Export.table.toDrive({
  collection: export_csv,
  description: 'SRER_Metrics_True_vs_Predicted',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});



// =========================================================================
// VISUALIZATION
// =========================================================================

var pred_im = pred_bgr;

// Generate hillshade
var dem = ee.Image('USGS/3DEP/10m').clip(v_extent);
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
Map.addLayer(hillshade, {min: 0, max: 255}, 'Hillshade');


//HSV Hillshade method
// Turn prediction into a baked RGB image (outputs values 0-255)
var pred_visualized = pred_im.visualize({
  min: 0, 
  max: 80, 
  palette: ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']
});

// Normalize the hillshade to a 0.0 - 1.0 scale
var hillshade_norm = hillshade.divide(255.0);

// Normalize the RGB image to 0.0 - 1.0, THEN convert to HSV
var hsv_image = pred_visualized.divide(255.0).rgbToHsv();

// Recombine: Boost saturation and blend the hillshade with the original value
var draped_hsv = ee.Image.cat([
  hsv_image.select('hue'),
  
  // Boost the color saturation by 20% and clamp it so it doesn't exceed 1.0
  hsv_image.select('saturation').multiply(1.2).clamp(0, 1), 
  
  // Blend the original color's brightness (40%) with the hillshade (60%)
  hsv_image.select('value').multiply(0.4).add(hillshade_norm.multiply(0.6)) 
]);

// Convert back to RGB (outputs 0.0 - 1.0 range, which Map.addLayer handles perfectly)
var final_draped_rgb = draped_hsv.hsvToRgb();

Map.addLayer(final_draped_rgb, {}, 'Draped BGR (HSV Method)');

// Visualizing just the BGR prediction for reference
//Map.addLayer(pred_im, {min:5, max:80, palette:['#487d4a', '#3EB489', '#B8EF80', '#FAC05B', '#964B00']}, 'Predicted');
//Map.addLayer(pred_im, {min:5, max:80, palette:['#00204d', '#31446b', '#666970', '#958f78', '#cbb67a', '#ffea46']}, 'Predicted');
//Map.addLayer(pred_im, {min:5, max:80, palette:['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']}, 'Predicted');
Map.addLayer(pred_im, {min:15, max:75, palette:['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']}, 'Predicted');

Map.addLayer(bounds_fc);
