// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area(1).divide(10000)); 
                     });

// =========================================================================
// PART 1: STATIC GRID GENERATION
// =========================================================================

// Extract projection from Sentinel-2 to build the 10m grid
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .first()
  .select('B2')
  .projection();

var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// Removed the v_srer_polys filter here so it covers all drone prints
var focused_grid = sent2_grid.filterBounds(v_foot_prints);

// Mask only the drone footprints
var valid_area_mask = ee.Image.constant(0).paint(v_foot_prints, 1);

var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1, 
  crs: projSent2.crs(),
  tileScale: 4
});

var final_grid = grid_overlap.filter(ee.Filter.gte('mean', 0.99));

var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});

// Use saveFirst with outer: true to keep grid cells that don't match a polygon (Left Outer Join)
var v_saveFirstJoin = ee.Join.saveFirst({
  matchKey: 'poly_match',
  outer: true
});

// The clean base grid conditionally setting polygon attributes
var base_grid = v_saveFirstJoin.apply(final_grid, v_srer_polys, v_spatial_filter)
  .map(function (ft){
    var match = ft.get('poly_match');
    
    // If there is a match, get the property. If not, return a blank string.
    var plantComm = ee.Algorithms.If(match, ee.Feature(match).get('Plant_Comm'), '');
    var pasture = ee.Algorithms.If(match, ee.Feature(match).get('Pasture'), '');
    var sDesc = ee.Algorithms.If(match, ee.Feature(match).get('S_Desc'), '');

    return ee.Feature(ft.geometry()).set({
      'Plant_Comm': plantComm,
      'Pasture': pasture,
      'S_Desc': sDesc
    });
  });

// =========================================================================
// PART 2: SEASONAL DUPLICATION & SHAPEFILE EXPORT
// =========================================================================

// Create the May cells
var may_grid = base_grid.map(function(f) {
  return f.set('month', 'May');
});

// Create the September cells
var sep_grid = base_grid.map(function(f) {
  return f.set('month', 'Sept');
});

// Merge them into one master dataset
var final_shapefile_grid = may_grid.merge(sep_grid);

// Transform geometries to UTM Zone 12N
var final_utm_grid = final_shapefile_grid.map(function(f) { 
  return f.transform('EPSG:32612', 0.05); 
});

// Export DIRECTLY to Google Drive as a .shp
Export.table.toDrive({
  collection: final_utm_grid,
  description: 'SRER_Grid',
  folder: 'GEE_Downloads',
  fileFormat: 'SHP'
});
