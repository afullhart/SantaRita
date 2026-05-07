// =========================================================================
// USER INPUTS & DATES
// =========================================================================

// Dry Season (Pre-Monsoon) Dates
var may_start = '2019-05-26';
var may_end   = '2019-05-31';

// Post-Monsoon Dates
var sep_start = '2019-09-24';
var sep_end   = '2019-09-30';

// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

// Safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

// =========================================================================
// SENTINEL-2 PREDICTOR EXTRACTION (REUSABLE FUNCTION)
// =========================================================================

// Extract projection from a reference image using the May dates defined above
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .filterDate(may_start, may_end)
  .first().select('B2').projection();

// Reusable function to build the predictor stack for any date range
function buildS2Composite(startDate, endDate) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(v_extent)                     
    .filterDate(startDate, endDate);                

  var sent2_im = sent2_ic
    .mosaic()
    .clip(v_extent)
    .select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12']) 
    .multiply(0.0001) 
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');

  var mcari = sent2_im.expression(
      '((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
        'B3': sent2_im.select('B3'), 
        'B4': sent2_im.select('B4'), 
        'B5': sent2_im.select('B5')  
  }).rename('MCARI');

  var bsi = sent2_im.expression(
      '((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
        'B2':  sent2_im.select('B2'),
        'B4':  sent2_im.select('B4'),
        'B8':  sent2_im.select('B8'),
        'B11': sent2_im.select('B11')
  }).rename('BSI');

  var nbr2 = sent2_im.normalizedDifference(['B11', 'B12']).rename('NBR2');

  return sent2_im.addBands([ndvi, mcari, bsi, nbr2]);
}

// Build the two distinct seasonal predictor maps
var sent2_may = buildS2Composite(may_start, may_end);
var sent2_sep = buildS2Composite(sep_start, sep_end);

var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2'];


// =========================================================================
// K-FOLDS CROSS VALIDATION (K=5) - UNIFIED MODEL, STRATIFIED ERROR
// =========================================================================
var k_folds = 5;
var v_list_seeds = ee.List([123, 456, 789, 111, 333]);
var fold_list = ee.List.sequence(0, 4);

var fc_folds = fc.randomColumn('random', 123).map(function(ft) {
  return ft.set('fold', ee.Number(ft.get('random')).multiply(k_folds).floor());
});

// Reusable function to run CV and return stratified metrics
function runCV(targetProp) {
  var fold_metrics = fold_list.map(function(fold) {
    var i = ee.Number(fold);
    var current_seed = v_list_seeds.get(i);
    
    // Train on the COMBINED dataset
    var train_fc = fc_folds.filter(ee.Filter.neq('fold', i));
    var test_fc = fc_folds.filter(ee.Filter.eq('fold', i));
    
    var gtb = ee.Classifier.smileGradientTreeBoost({
      numberOfTrees: 300, 
      shrinkage: 0.01,
      samplingRate: 0.7, 
      maxNodes: 12,       
      seed: current_seed
    })
    .setOutputMode('REGRESSION')
    .train({
      features: train_fc, 
      classProperty: targetProp, 
      inputProperties: inputProps
    });
    
    // Predict the entire COMBINED test set
    var tested = test_fc.classify({
      classifier: gtb, 
      outputName: 'predicted'
    });
    
    // Calculate the squared differences for all points
    var with_sq_err = tested.map(function(ft) {
      var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
      return ft.set('sq_diff', diff.multiply(diff));
    });

    // Splitting the errors by month to return two separate metrics per fold
    var may_mse = with_sq_err.filter(ee.Filter.eq('Month', 'May'))
      .reduceColumns({reducer: ee.Reducer.mean(), selectors: ['sq_diff']}).get('mean');
      
    var sep_mse = with_sq_err.filter(ee.Filter.eq('Month', 'Sept'))
      .reduceColumns({reducer: ee.Reducer.mean(), selectors: ['sq_diff']}).get('mean');
    
    return ee.Dictionary({
      'May_RMSE': ee.Number(may_mse).sqrt(),
      'Sept_RMSE': ee.Number(sep_mse).sqrt()
    });
  });
  
  // Extract lists of the 5 RMSEs for each month
  var may_rmse_list = fold_metrics.map(function(d) { return ee.Dictionary(d).get('May_RMSE'); });
  var sep_rmse_list = fold_metrics.map(function(d) { return ee.Dictionary(d).get('Sept_RMSE'); });
  
  // Return the median across the 5 folds
  return ee.Dictionary({
    'May_Median': ee.List(may_rmse_list).reduce(ee.Reducer.median()),
    'Sept_Median': ee.List(sep_rmse_list).reduce(ee.Reducer.median())
  });
}

print('--- UNIFIED K-FOLDS RMSE (STRATIFIED TESTING ERROR) ---');
var cv_bgr = runCV('BGR');
print('BGR RMSE -> May:', cv_bgr.get('May_Median'), '| Sept:', cv_bgr.get('Sept_Median'));

var cv_lpi = runCV('LPI');
print('LPI RMSE -> May:', cv_lpi.get('May_Median'), '| Sept:', cv_lpi.get('Sept_Median'));

var cv_mft = runCV('MFT');
print('MFT RMSE -> May:', cv_mft.get('May_Median'), '| Sept:', cv_mft.get('Sept_Median'));
print('-------------------------------------------------------');


// =========================================================================
// FINAL MODEL TRAINING (UNIFIED)
// =========================================================================

var regularized_params = {
  numberOfTrees: 300, 
  shrinkage: 0.01,
  samplingRate: 0.7, 
  maxNodes: 12, 
  seed: 123
};

// Train ONE set of models using the combined dataset (fc)
var model_bgr = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'BGR', inputProperties: inputProps});

var model_lpi = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'LPI', inputProperties: inputProps});

var model_mft = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'MFT', inputProperties: inputProps});


// -------------------------------------------------------------------------
// CALCULATE & PRINT FINAL MODEL (STRATIFIED TRAINING ERROR)
// -------------------------------------------------------------------------

// Helper function to calculate the fitted RMSE, filtered by a specific month
function getStratifiedTrainingRMSE(trainedModel, dataset, targetProp, monthFilter) {
  
  // Filter the dataset down to the month of interest BEFORE calculating error
  var filtered_dataset = dataset.filter(ee.Filter.eq('Month', monthFilter));
  
  var classified = filtered_dataset.classify({
    classifier: trainedModel, 
    outputName: 'predicted'
  });
  
  var mse = classified.map(function(ft) {
    var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
    return ft.set('sq_diff', diff.multiply(diff));
  }).reduceColumns({
    reducer: ee.Reducer.mean(),
    selectors: ['sq_diff']
  }).get('mean');
  
  return ee.Number(mse).sqrt();
}

print('--- UNIFIED FITTED MODEL RMSE (STRATIFIED TRAINING ERROR) ---');
print('Fitted BGR RMSE -> May:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'Sept'));
print('Fitted LPI RMSE -> May:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'Sept'));
print('Fitted MFT RMSE -> May:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'Sept'));
print('------------------------------------------------');


// -------------------------------------------------------------------------
// CLASSIFY SEASONAL IMAGES
// -------------------------------------------------------------------------

// Classify May and Sept composite images using the SAME unified models
var p_bgr_may = sent2_may.classify(model_bgr).rename('Pred_BGR');
var p_lpi_may = sent2_may.classify(model_lpi).rename('Pred_LPI');
var p_mft_may = sent2_may.classify(model_mft).rename('Pred_MFT');
var combined_preds_may = ee.Image([p_bgr_may, p_lpi_may, p_mft_may]);

var p_bgr_sep = sent2_sep.classify(model_bgr).rename('Pred_BGR');
var p_lpi_sep = sent2_sep.classify(model_lpi).rename('Pred_LPI');
var p_mft_sep = sent2_sep.classify(model_mft).rename('Pred_MFT');
var combined_preds_sep = ee.Image([p_bgr_sep, p_lpi_sep, p_mft_sep]);


// =========================================================================
// SEASON-AWARE SAMPLING & CSV EXPORT
// =========================================================================

// Convert the grid features to point geometries located at their centers
var fc_centers = fc.map(function(ft) {
  return ft.centroid(1); 
});

// 1. Re-split the geometric centers by month
var centers_may = fc_centers.filter(ee.Filter.eq('Month', 'May'));
var centers_sep = fc_centers.filter(ee.Filter.eq('Month', 'Sept'));

// 2. Sample the respective prediction maps using the accurately dated points
var sampled_may = combined_preds_may.sampleRegions({
  collection: centers_may,
  properties: ['BGR', 'LPI', 'MFT', 'Month'],
  scale: projSent2.nominalScale(),
  tileScale: 4
});

var sampled_sep = combined_preds_sep.sampleRegions({
  collection: centers_sep,
  properties: ['BGR', 'LPI', 'MFT', 'Month'], 
  scale: projSent2.nominalScale(), 
  tileScale: 4
});

// 3. Merge them back together for a clean, unified export
var sampled_data_merged = sampled_may.merge(sampled_sep);

// 4. Format for CSV
var export_csv = sampled_data_merged.map(function(ft) {
  return ee.Feature(null, { 
    'Month': ft.get('Month'),
    'True_BGR': ft.get('BGR'),
    'Predicted_BGR': ft.get('Pred_BGR'),
    'True_LPI': ft.get('LPI'),
    'Predicted_LPI': ft.get('Pred_LPI'),
    'True_MFT': ft.get('MFT'),
    'Predicted_MFT': ft.get('Pred_MFT')
  });
});

Export.table.toDrive({
  collection: export_csv,
  description: 'SRER_Metrics_True_vs_Predicted_Unified',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});


// =========================================================================
// VISUALIZATION
// =========================================================================

// Generate hillshade (This is the only layer set to show by default)
var dem = ee.Image('USGS/3DEP/10m').clip(v_extent);
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
Map.addLayer(hillshade, {min: 0, max: 255}, 'Hillshade');

var hillshade_norm = hillshade.divide(255.0);

// Visualize May Predictions
var pred_im_may = combined_preds_may.select('Pred_BGR');
var hsv_image_may = pred_im_may.visualize({
  min: 0, max: 80, palette: ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']
}).divide(255.0).rgbToHsv();

var draped_hsv_may = ee.Image.cat([
  hsv_image_may.select('hue'),
  hsv_image_may.select('saturation').multiply(1.2).clamp(0, 1), 
  hsv_image_may.select('value').multiply(0.4).add(hillshade_norm.multiply(0.6)) 
]);

Map.addLayer(draped_hsv_may.hsvToRgb(), {}, 'Draped BGR (May)', false);
Map.addLayer(pred_im_may, {min:15, max:75, palette:['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']}, 'Predicted BGR (May)', false);

// Visualize September Predictions
var pred_im_sep = combined_preds_sep.select('Pred_BGR');
var hsv_image_sep = pred_im_sep.visualize({
  min: 0, max: 80, palette: ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']
}).divide(255.0).rgbToHsv();

var draped_hsv_sep = ee.Image.cat([
  hsv_image_sep.select('hue'),
  hsv_image_sep.select('saturation').multiply(1.2).clamp(0, 1), 
  hsv_image_sep.select('value').multiply(0.4).add(hillshade_norm.multiply(0.6)) 
]);

Map.addLayer(draped_hsv_sep.hsvToRgb(), {}, 'Draped BGR (Sept)', false);
Map.addLayer(pred_im_sep, {min:15, max:75, palette:['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027']}, 'Predicted BGR (Sept)', false);

// Boundary Layer
Map.addLayer(bounds_fc, {}, 'SR_bounds', false);
