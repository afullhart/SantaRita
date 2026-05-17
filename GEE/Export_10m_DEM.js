// =========================================================================
// SETUP BOUNDS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// DEM PRE-PROCESSING & STRIPING REMOVAL
// =========================================================================
// 1. Grab the target 10m projection from the original DEM to ensure alignment
var proj10m = ee.Image('USGS/3DEP/10m').projection();

// 2. Load the 1m DEM Collection and filter to the study area
var dem_1m_col = ee.ImageCollection('USGS/3DEP/1m').filterBounds(bounds_geom);

// 3. Extract the native 1m projection
var proj1m = dem_1m_col.first().projection();

// 4. Flatten the swaths using a mean reduction (better than flat mosaic for striping)
var dem_1m_raw = dem_1m_col.mean().setDefaultProjection(proj1m);

// 5. Apply a Gaussian low-pass filter to blur out swath edge calibration artifacts
var smoothKernel = ee.Kernel.gaussian({
  radius: 2,       // 2-meter radius to blend sharp vertical edges
  units: 'meters',
  sigma: 1
});
var dem_1m_clean = dem_1m_raw.convolve(smoothKernel);

// 6. Aggregate the stripe-free 1m DEM to 10m using spatial averaging
var dem_10m_resampled = dem_1m_clean
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
