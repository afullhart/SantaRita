// =========================================================================
// REGIONAL TIME-SERIES: EXACT 10m EXPORT TASK
// =========================================================================

var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_S2');
var bounds_geom = bounds_fc.first().geometry().bounds();

var dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled'); 
var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

var core_training_fc = fc.filter(ee.Filter.notNull(inputProps));
var hwr_training_fc = core_training_fc.filter(ee.Filter.notNull(['Herb_Woody_Ratio'])).map(function(ft) {
  var log_hwr = ee.Number(ft.get('Herb_Woody_Ratio')).add(1).log(); 
  return ft.set('Log_HWR', log_hwr);
});

var hyperpars = { numberOfTrees: 400, shrinkage: 0.05, samplingRate: 0.7, maxNodes: 32, loss: 'Huber', seed: 123 };

var model_bgr = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'MFT', inputProperties: inputProps});
var model_hwr = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: hwr_training_fc, classProperty: 'Log_HWR', inputProperties: inputProps});

var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(bounds_geom).first().select('B2').projection();

function buildS2Composite(startDate, endDate) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(bounds_geom)                    
    .filterDate(startDate, endDate)
    .map(function(img) { return img.set('has_aux_bands', img.bandNames().contains('MSK_CLDPRB')); })
    .filter(ee.Filter.eq('has_aux_bands', true))
    .map(function(img) {
      var probMask = img.select('MSK_CLDPRB').lt(20);
      var scl = img.select('SCL');
      var sclMask = scl.neq(8).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(3));
      var blueMask = img.select('B2').lt(2500);
      var masterMask = probMask.and(sclMask).and(blueMask);
      var maskedImg = img.updateMask(masterMask);

      var sz_num = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180);
      var sa_num = ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180);
      var cosZ = ee.Image.constant(sz_num.cos());
      var sinZ = ee.Image.constant(sz_num.sin());
      var sa_img = ee.Image.constant(sa_num);
      var cosS = terrain_slope.cos();
      var sinS = terrain_slope.sin();
      var cosAzAsp = sa_img.subtract(terrain_aspect).cos();
      var illumination = cosZ.multiply(cosS).add(sinZ.multiply(sinS).multiply(cosAzAsp)).rename('illumination');
            
      return maskedImg.addBands(illumination);
    });                

  var sent2_median = sent2_ic.median().clip(bounds_geom);
  var optical_bands = sent2_median.select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12'])
    .multiply(0.0001).setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = optical_bands.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mcari = optical_bands.expression('((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {'B3': optical_bands.select('B3'), 'B4': optical_bands.select('B4'), 'B5': optical_bands.select('B5')}).rename('MCARI');
  var bsi = optical_bands.expression('((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {'B2': optical_bands.select('B2'), 'B4': optical_bands.select('B4'), 'B8': optical_bands.select('B8'), 'B11': optical_bands.select('B11')}).rename('BSI');
  var nbr2 = optical_bands.normalizedDifference(['B11', 'B12']).rename('NBR2');

  return optical_bands.addBands([ndvi, mcari, bsi, nbr2, sent2_median.select('illumination'), terrain_slope.rename('slope'), terrain_aspect.rename('aspect')]);
}

print('Building the 10m export task... Check your Tasks tab.');

// -------------------------------------------------------------------------
// EXTRACT 10M TIME-SERIES FOR EXPORT
// -------------------------------------------------------------------------
var regionalTimeSeriesData = cloud_windows.map(function(window) {
  var sDate = ee.String(window.get('Start_Date'));
  var eDate = ee.Date(sDate).advance(7, 'day');
  var s2_img = buildS2Composite(sDate, eDate);

  var p_bgr = s2_img.classify(model_bgr).rename('BGR_pct');
  var p_lpi = s2_img.classify(model_lpi).rename('LPI_pct');
  var p_mft = s2_img.classify(model_mft).rename('Fetch_m');
  var p_hwr = s2_img.classify(model_hwr).rename('Log_HWR');

  var combined_preds = ee.Image.cat([p_bgr, p_lpi, p_mft, p_hwr]);

  // EXACT 10m REDUCTION
  var regional_mean = combined_preds.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: bounds_geom,
    scale: 10,       // Full native resolution
    tileScale: 16,   // Max memory chunking
    maxPixels: 1e13
  });

  return ee.Feature(null, {
    'Date': window.get('Start_Date'),
    'System_Time': window.get('system:time_start'),
    'Mean_BGR_pct': regional_mean.get('BGR_pct'),
    'Mean_LPI_pct': regional_mean.get('LPI_pct'),
    'Mean_Fetch_m': regional_mean.get('Fetch_m'),
    'Mean_Log_HWR': regional_mean.get('Log_HWR') 
  });
});

var validRegionalData = regionalTimeSeriesData.filter(ee.Filter.notNull(['Mean_BGR_pct']));

// EXPORT DIRECTLY TO DRIVE
Export.table.toDrive({
  collection: validRegionalData,
  description: 'SRER_Regional_Average_TimeSeries',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV',
  selectors: ['Date', 'System_Time', 'Mean_BGR_pct', 'Mean_LPI_pct', 'Mean_Fetch_m', 'Mean_Log_HWR']
});
