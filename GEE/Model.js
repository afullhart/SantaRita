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

// Mosaic, clip, select bands (added B5), and APPLY SCALE FACTOR (0.0001)
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


var v_model = ee.Classifier.smileRandomForest({
  numberOfTrees: 500,
  minLeafPopulation: 5,
  seed: 123,
  }).setOutputMode('REGRESSION')
  .train({
    features: fc,
      classProperty: 'BGR',
      inputProperties: ['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI'] 
  });

var pred_im = sent2_im.classify(v_model);

// Generate hillshade
var dem = ee.Image('USGS/3DEP/10m').clip(v_extent);
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
Map.addLayer(hillshade, {min: 0, max: 255}, 'Hillshade');


Map.addLayer(pred_im, {min:0, max:80, palette:['#487d4a', '#3EB489', '#B8EF80', '#FAC05B', '#964B00']});
//Map.addLayer(pred_im, {min:0, max:1, palette:['#487d4a', '#3EB489', '#FAC05B', '#964B00']});





