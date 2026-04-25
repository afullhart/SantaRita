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
// Create binary mask globally (Assuming 3 is bare ground)
var binary = classified_img.eq(3).selfMask();

Map.addLayer(binary, null, 'binary');

var connected = binary.connectedComponents({
  connectedness:ee.Kernel.square(1),
  maxSize:1024
});

Map.addLayer(connected, null, 'connected');

//Calculate the sizes of those patches globally (in pixels)
var patch_sizes = connected.reduceConnectedComponents({
  reducer: ee.Reducer.count(),
  maxSize:1024
});

//Native image
//var native_proj = v_classified_may.projection();
//var native_patch_sizes_reproj = patch_sizes.reproject(native_proj);
//var patch_sizes = native_patch_sizes_reproj;

var visParamsMagma = {
  min: 1,
  max: 30000,
  palette: ['440154', '414487', '2A788E', '22A884', '7AD151', 'FDE725']};

Map.addLayer(patch_sizes, visParamsMagma, 'patch_sizes');

// Convert global pixel count to physical area
var patch_area = patch_sizes.multiply(ee.Image.pixelArea()).rename('max_patch_area');

var visParamsMagma = {
  min: 1,
  max: 10000,
  palette: ['440154', '414487', '2A788E', '22A884', '7AD151', 'FDE725']};

Map.addLayer(patch_area, visParamsMagma, 'patch_area');

// --- 2. The Crucial Step: reduceRegions ---
// This calculates the max patch area for every grid cell efficiently in parallel
var lpi_reduced = patch_area.reduceRegions({
  collection:grid_collection,
  reducer:ee.Reducer.max(), 
  scale:0.05,               
  tileScale:16  // This prevents memory errors during reduction
});

print(lpi_reduced.first());


