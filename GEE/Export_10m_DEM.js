// =========================================================================
// SETUP BOUNDS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// DEM PURE RESAMPLING (NO SMOOTHING)
// =========================================================================
// Load the high-resolution 1m DEM ImageCollection and filter to the boundary
var dem_1m_col = ee.ImageCollection('USGS/3DEP/1m').filterBounds(bounds_geom);

// Extract the native 1m projection system
var proj1m = dem_1m_col.first().projection();

// Mosaic the collection and explicitly enforce the native 1m projection
var dem_1m = dem_1m_col.mosaic().setDefaultProjection(proj1m);

// Aggregate 1m to 10m using a true spatial mean reduction
// CRITICAL: We skip manual re-projections here to block grid aliasing errors.
var dem_10m_resampled = dem_1m
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 1024
  })
  .clip(bounds_geom);

// =========================================================================
// LIVE MAP VISUALIZATION (QUALITY CHECK)
// =========================================================================
var hillshade_check = ee.Terrain.hillshade(dem_10m_resampled, 270, 45);

Map.centerObject(bounds_geom, 12);
Map.addLayer(hillshade_check, {min: 0, max: 255}, 'Stripe-Free Asset Hillshade Preview');

// =========================================================================
// BULLETPROOF EXPORT WITH EXPLICIT CRS DEFINITION
// =========================================================================
Export.image.toAsset({
  image: dem_10m_resampled.rename('elevation'),
  description: 'SR_10m_DEM_Pristine',
  assetId: 'projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled',
  region: bounds_geom,
  scale: 10,
  // THE CRITICAL FIX: Direct string definitions ensure pristine backend translation
  crs: 'EPSG:4269', 
  pyramidingPolicy: { 'elevation': 'mean' }
});
