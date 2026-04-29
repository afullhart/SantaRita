var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area().divide(10000)); 
});

// Because bounds_geom is a Geometry, .bounds() safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)                     // FILTER LOCATION FIRST
  .filterDate('2019-05-01', '2019-05-31');    // FILTER DATE SECOND

var v_sent2_joined_grids = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_grid_joined');

print(v_sent2_joined_grids.size());

//function calculateLPI(classified_img, grid_collection, period_name) {
var classified_img = v_classified_may;
var grid_collection = v_sent2_joined_grids;
var period_name = ee.String('May');

var binary = classified_img.eq(3).selfMask();

//function gridcell_iter_fn(){
var native_proj = v_classified_may.projection();

function feat_iteration_fn(ft){

  // Multiply the binary mask with the area of each pixel in square meters
  var v_area_image = binary.multiply(ee.Image.pixelArea());
  var v_area_ft = ft.area();
  // Reduce the image to calculate the total area covered by the specific class
  var v_area = v_area_image.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: ft.geometry(),
    scale: 0.05,
    maxPixels: 1e13
   }).get('classification');
  var v_pct_area = ee.Number(v_area).divide(v_area_ft).multiply(100);

  // 1. Convert the binary 5cm raster mask directly into vector polygons
  var patch_vectors = binary.reduceToVectors({
    reducer: ee.Reducer.countEvery(),
    geometry: ft.geometry(),        // The area to process
    scale: 0.05,               // Force native 5cm scale
    crs: native_proj.crs(),
    geometryType: 'polygon',
    eightConnected: true,
    labelProperty: 'class_val',
    maxPixels: 1e13,           // Allow processing over a massive area
    tileScale: 16              // Prevent memory crashes during tile extraction
  });

  // 2. Calculate the true physical area of every native-resolution polygon
  var patches_with_area = patch_vectors.map(function(ft){
    return ft.set('area_sqm', ft.area(0.05)); // Calculate area with 5cm error margin
  });

  // Extract all 'area_sqm' values directly into an ee.List of numbers
  var areas = patches_with_area.aggregate_array('area_sqm');
  
  // Find the max, defaulting to 0 if no patches intersect
  var max_area = ee.Algorithms.If(
    areas.length().gt(0),
    areas.reduce(ee.Reducer.max()),
    0
  );

  function Get_Mean_Nearest_Bground_Pixel(v_image, v_points) {
    var v_distance = v_image.fastDistanceTransform().sqrt().multiply(ee.Image.pixelArea().sqrt()).rename("distance");
    v_points = v_distance.reduceRegions(v_points, ee.Reducer.first().setOutputs(["distance"]));
    return ee.Number(v_points.reduceColumns(ee.Reducer.mean(),['distance']).get('mean'));
  }

  var N_PTS = 1000;
  var v_rnd = ee.FeatureCollection.randomPoints(ft.geometry(), N_PTS, 1234, 0.05);
  var v_nearestMeanValues = Get_Mean_Nearest_Bground_Pixel(binary, v_rnd);

  return ft.set('LPI', max_area, 'BGR', v_pct_area, 'MFT', v_nearestMeanValues);
}

var final_grids = v_sent2_joined_grids.map(feat_iteration_fn);

print(final_grids.first());

// 5. Export the results
Export.table.toDrive({
  collection: final_grids,
  description: 'Native_5cm_LPI_Vectors',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});


Export.table.toAsset({
  collection: final_grids,
  assetId: 'projects/ee-andrewfullhart/assets/SR_s2_model_grid',
  description: 'ftv_sentinel2_grid_srer_slud'
});

