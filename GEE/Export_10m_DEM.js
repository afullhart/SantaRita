// =========================================================================
// SETUP BOUNDS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// PURE DEM RESAMPLING (NO SMOOTHING)
// =========================================================================
// 1. Grab the target standard 10m projection system to ensure grid alignment
var proj10m = ee.Image('USGS/3DEP/10m').projection();

// 2. Load the high-resolution 1m DEM ImageCollection and filter to the boundary
var dem_1m_col = ee.ImageCollection('USGS/3DEP/1m').filterBounds(bounds_geom);

// 3. Extract the native 1m projection system
var proj1m = dem_1m_col.first().projection();

// 4. Mosaic the collection and explicitly apply the native 1m projection
var dem_1m = dem_1m_col.mosaic().setDefaultProjection(proj1m);

// 5. Aggregate 1m to 10m using a true spatial mean reduction (No smoothing filters)
var dem_10m_resampled = dem_1m
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 1024
  })
  .reproject({
    crs: proj10m
  })
  .clip(bounds_geom);

// =========================================================================
// ON-THE-FLY HILLSHADE VISUALIZATION
// =========================================================================
// Compute the standard terrain hillshade directly from the raw 10m resampled dataset
var hillshade = ee.Terrain.hillshade(dem_10m_resampled, 270, 45);

// Center the interactive map framework over your study bounds
Map.centerObject(bounds_geom, 12);

// Display the true aggregated elevation data using a multi-color gradient
Map.addLayer(dem_10m_resampled, {
  min: 1000,
  max: 2000,
  palette: ['#1a9850', '#ffffbf', '#d73027']
}, 'Raw Resampled 10m DEM (Elevation)');

// Display the un-smoothed hillshade terrain profile
Map.addLayer(hillshade, {min: 0, max: 255}, 'Raw Resampled 10m DEM (Hillshade)');

// =========================================================================
// EXPORT TO SINGLE-BAND ASSET
// =========================================================================

Export.image.toAsset({
  image: dem_10m_resampled.rename('elevation'),
  description: 'SR_10m_DEM_Resampled',
  assetId: 'projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled',
  region: bounds_geom,
  scale: 10,
  crs: proj10m.crs().getInfo(),
  maxPixels: 1e13
});
