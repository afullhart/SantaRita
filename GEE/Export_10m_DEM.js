// =========================================================================
// SETUP BOUNDS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// DEM PRE-PROCESSING & RESAMPLING
// =========================================================================
// 1. Grab the target 10m projection from the original DEM to ensure alignment
var proj10m = ee.Image('USGS/3DEP/10m').projection();

// 2. Load the 1m DEM Collection, filter to the ROI, and mosaic
var dem_1m_col = ee.ImageCollection('USGS/3DEP/1m').filterBounds(bounds_geom);

// 3. Extract the native 1m projection and apply it to the mosaic
var proj1m = dem_1m_col.first().projection();
var dem_1m = dem_1m_col.mosaic().setDefaultProjection(proj1m);

// 4. Aggregate 1m to 10m using spatial averaging and reproject
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
// EXPORT TO ASSET
// =========================================================================
Export.image.toAsset({
  image: dem_10m_resampled,
  description: 'SR_10m_DEM',
  assetId: 'projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled',
  region: bounds_geom,
  scale: 10,
  crs: proj10m.crs().getInfo(),
  maxPixels: 1e13
});
