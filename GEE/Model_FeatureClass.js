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
  .multiply(0.0001) // Converts integer DNs to true surface reflectance (0.0 - 1.0)
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

// Append both computed indices to the image
sent2_im = sent2_im.addBands([ndvi, mcari]);

// Generate the base Sentinel-2 pixel grid
var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// =========================================================================
// FAST FRACTIONAL OVERLAP WORKFLOW
// =========================================================================

// Pre-filter the grid so we don't calculate overlap on empty space
var focused_grid = sent2_grid.filterBounds(v_srer_polys).filterBounds(v_foot_prints);

// Create a high-resolution binary mask (1) where the footprints and polygons overlap
var footprint_mask = ee.Image.constant(0).paint(v_foot_prints, 1);
var poly_mask = ee.Image.constant(0).paint(v_srer_polys, 1);
var valid_area_mask = footprint_mask.and(poly_mask);

// Calculate the exact percentage of overlap for every grid cell.
// By reducing a 1m resolution mask over the grid, the 'mean' equals the area fraction.
var grid_overlap = valid_area_mask.reduceRegions({
  collection: focused_grid,
  reducer: ee.Reducer.mean(),
  scale: 1, // 1m scale ensures highly accurate sub-pixel area math
  crs: projSent2.crs(),
  tileScale: 4
});

// Filter for strict 100% overlap (using 0.99 to account for tiny floating-point rounding)
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

// Extract the Sentinel-2 bands, NDVI, and MCARI for these grid polygons
var bandsToExtract = sent2_im.select(['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI']);

var v_sent2_with_bands = bandsToExtract.reduceRegions({
  collection: v_sent2_joined_grids,
  reducer: ee.Reducer.mean(),
  scale: projSent2.nominalScale(),
  crs: projSent2.crs(),
  tileScale: 4
});

// Clean up any potential grids that might have fallen on masked image pixels (null values)
v_sent2_with_bands = v_sent2_with_bands.filter(ee.Filter.notNull(['B2', 'B3', 'B4', 'B5', 'B8', 'NDVI', 'MCARI']));

// =========================================================================
// EXPORT
// =========================================================================

Export.table.toAsset({
  collection: v_sent2_with_bands,
  assetId: 'projects/ee-andrewfullhart/assets/SR_s2_grid_joined',
  description: 'ftv_sentinel2_grid_srer_slud'
});

