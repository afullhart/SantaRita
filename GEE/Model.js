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

// Mosaic, clip, and select the bands we need right away (including B8 for NDVI)
var sent2_im = sent2_ic
  .mosaic()
  .clip(v_extent)
  .select(['B2', 'B3', 'B4', 'B8'])
  .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

// Calculate NDVI using B8 (NIR) and B4 (Red) and add it as a new band
var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');
sent2_im = sent2_im.addBands(ndvi);


var v_model = ee.Classifier.smileRandomForest({
  numberOfTrees: 500,
  minLeafPopulation: 5,
  seed: 123,
  }).setOutputMode('REGRESSION')
  .train({
    features: fc,
      classProperty: 'BGR', 
      inputProperties: ['B2', 'B3', 'B4', 'B8', 'NDVI']
  });

print(v_model);

var pred_im = sent2_im.classify(v_model);

Map.addLayer(pred_im, {min:0, max:80, palette:['#487d4a', '#3EB489', '#FAC05B', '#964B00']});
//Map.addLayer(pred_im, {min:0, max:1, palette:['#487d4a', '#3EB489', '#FAC05B', '#964B00']});



