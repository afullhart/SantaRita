var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

// Classified drone images for month of May and Sept
var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area(1).divide(10000)); 
                     });

// Safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

// =========================================================================
// STATIC GRID GENERATION (Run once for both months)
// =========================================================================

// Extract projection from a single image in the collection to build the grid
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .first()
  .select('B2')
  .projection();

// Generate the base Sentinel-2 pixel grid
var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// Pre-filter the grid so we don't calculate overlap on empty space
var focused_grid = sent2_grid.filterBounds(v_srer_polys).filterBounds(v_foot_prints);

// Create a high-resolution binary mask (1) where the footprints and polygons overlap
var footprint_mask = ee.Image.constant(0).paint(v_foot_prints, 1);
var poly_mask = ee.Image.constant(0).paint(v_srer_polys, 1);
var valid_area_mask = footprint_mask.and(poly_mask);

// Calculate the exact percentage of overlap for every grid cell.
var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1, 
  crs: projSent2.crs(),
  tileScale: 4
});

// Filter for strict 100% overlap
var final_grid = grid_overlap.filter(ee.Filter.gte('mean', 0.99));

// Transfer the polygon attributes to your perfectly overlapping grid cells
var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});
var v_saveAllJoin = ee.Join.saveAll({matchesKey:'polys'});

var v_sent2_joined_grids = v_saveAllJoin.apply(final_grid, v_srer_polys, v_spatial_filter)
  .map(function (ft){
    var ft1 = ee.Feature(ee.List(ft.get('polys')).get(0));
    return ft.set({
      'Plant_Comm': ft1.get('Plant_Comm'),
      'Pasture': ft1.get('Pasture'),
      'Transect': ft1.get('Transect'),
      'Utility': ft1.get('Utility'),
      'S_Desc': ft1.get('S_Desc'),
      'Exclosure': ft1.get('Exclosure'),
      'area_ha': ft1.get('area_ha'),
      'polys': null // clean up the temporary join list
    });
  });

// =========================================================================
// FUNCTION TO PROCESS S2 DATA BY DATE RANGE
// =========================================================================

function extractS2Data(startDate, endDate, monthLabel) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(v_extent)                     
    .filterDate(startDate, endDate);                

  // Mosaic, clip, select bands, and APPLY SCALE FACTOR (0.0001)
  var sent2_im = sent2_ic
    .mosaic()
    .clip(v_extent)
    .select(['B2', 'B3', 'B4', 'B5', 'B8'])
    .multiply(0.0001) 
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  // Calculate NDVI 
  var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');

  // Calculate MCARI 
  var mcari = sent2_im.expression(
      '((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
        'B3': sent2_im.select('B3'), 
        'B4': sent2_im.select('B4'), 
        'B5': sent2_im.select('B5')  
  }).rename('MCARI');

  // Append both computed indices to the image
  sent2_im = sent2_im.addBands([ndvi, mcari]);
  
  var bandsToExtract = sent2_im.select(['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI']);

  // Add the "Month" property to the grid features BEFORE reducing
  var gridWithMonth = v_sent2_joined_grids.map(function(feat) {
    return feat.set('Month', monthLabel);
  });

  // Extract the data using the grid that now contains the month identifier
  var extracted = bandsToExtract.reduceRegions({
    collection: gridWithMonth,
    reducer: ee.Reducer.mean(),
    scale: projSent2.nominalScale(),
    crs: projSent2.crs(),
    tileScale: 4
  });

  // Clean up any potential grids that might have fallen on masked image pixels (null values)
  return extracted.filter(ee.Filter.notNull(['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI']));
}

// =========================================================================
// EXECUTE AND MERGE
// =========================================================================

// Run the function for May
var may_data = extractS2Data('2019-05-26', '2019-05-31', 'May');

// Run the function for September (Adjust these dates to match your exact window)
var sep_data = extractS2Data('2019-09-01', '2019-09-30', 'Sept');

// Merge the two collections 
var combined_data = may_data.merge(sep_data);

// =========================================================================
// EXPORT
// =========================================================================

Export.table.toAsset({
  collection: combined_data,
  assetId: 'projects/ee-andrewfullhart/assets/SR_s2_grid',
  description: 'ftv_sentinel2_grid_srer_slud_may_sep'
});
