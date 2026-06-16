// =========================================================================
// REGIONAL TIME-SERIES: EXACT 30m EXPORT TASK (LANDSAT)
// =========================================================================

var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_landsat_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat');
var bounds_geom = bounds_fc.first().geometry().bounds();

var dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled'); 
var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// Landsat-specific predictor bands
var inputProps = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

// Filter for core valid data, explicitly ensuring Herb and Woody percent values exist
var core_training_fc = fc.filter(ee.Filter.notNull(inputProps.concat(['Herb_pct', 'Woody_pct'])));
var hwr_training_fc = core_training_fc.filter(ee.Filter.notNull(['Herb_Woody_Ratio'])).map(function(ft) {
  var log_hwr = ee.Number(ft.get('Herb_Woody_Ratio')).add(1).log(); 
  return ft.set('Log_HWR', log_hwr);
});

var hyperpars = { numberOfTrees: 400, shrinkage: 0.05, samplingRate: 0.7, maxNodes: 32, loss: 'Huber', seed: 123 };

var model_bgr = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'MFT', inputProperties: inputProps});
var model_herb = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'Herb_pct', inputProperties: inputProps});
var model_woody = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'Woody_pct', inputProperties: inputProps});
var model_hwr = ee.Classifier.smileGradientTreeBoost(hyperpars).setOutputMode('REGRESSION').train({features: hwr_training_fc, classProperty: 'Log_HWR', inputProperties: inputProps});

var projLandsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).first().select('SR_B2').projection();

// =========================================================================
// LANDSAT EXTRACTION PIPELINE & COMPOSITE GENERATOR
// =========================================================================
function maskClouds(qa) {
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  return dilatedCloud.and(cirrus).and(cloud).and(shadow);
}

function prepOLI(img) {
  var scaled = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2); 
  var qaMask = maskClouds(img.select('QA_PIXEL'));
  var finalMask = qaMask.and(scaled.select('SR_B5').gt(0.15)).focal_min({radius: 30, kernelType: 'square', units: 'meters'});

  var sz_num = ee.Number(img.get('SUN_ELEVATION')).subtract(90).multiply(-1).multiply(Math.PI / 180);
  var sa_num = ee.Number(img.get('SUN_AZIMUTH')).multiply(Math.PI / 180);
  var illumination = ee.Image.constant(sz_num.cos()).multiply(terrain_slope.cos()).add(ee.Image.constant(sz_num.sin()).multiply(terrain_slope.sin()).multiply(ee.Image.constant(sa_num).subtract(terrain_aspect).cos())).rename('illumination');

  var optical = scaled.rename(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']).addBands(illumination);
  return optical.updateMask(finalMask);
}

function prepTM(img) {
  var scaled = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2); 
  var qaMask = maskClouds(img.select('QA_PIXEL'));
  var finalMask = qaMask.and(scaled.select('SR_B4').gt(0.15)).focal_min({radius: 30, kernelType: 'square', units: 'meters'});

  var sz_num = ee.Number(img.get('SUN_ELEVATION')).subtract(90).multiply(-1).multiply(Math.PI / 180);
  var sa_num = ee.Number(img.get('SUN_AZIMUTH')).multiply(Math.PI / 180);
  var illumination = ee.Image.constant(sz_num.cos()).multiply(terrain_slope.cos()).add(ee.Image.constant(sz_num.sin()).multiply(terrain_slope.sin()).multiply(ee.Image.constant(sa_num).subtract(terrain_aspect).cos())).rename('illumination');

  var optical = scaled.rename(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']).addBands(illumination);
  return optical.updateMask(finalMask);
}

function buildLandsatComposite(valid_ids_str, year, month) {
  var fullIdList = ee.String(valid_ids_str).split(',');
  var indexList = fullIdList.map(function(id) {
    var parts = ee.String(id).split('/');
    return parts.get(parts.length().subtract(1));
  });

  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI);
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI);
  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM);
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM);

  var combined = l9.merge(l8).merge(l7).merge(l5);
  var ls_median = combined.median().clip(bounds_geom);

  // SERVER-SAFE SYNTHETIC MEDIAN REPLACEMENT FOR SEPTEMBER 2019
  var isTargetCond = ee.Number(year).eq(2019).and(ee.Number(month).eq(9));
  
  var aug31_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(bounds_geom)
    .filterDate('2019-08-31', '2019-09-01')
    .map(prepOLI); 
    
  var oct02_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(bounds_geom)
    .filterDate('2019-10-02', '2019-10-03')
    .map(prepOLI);

  var synthetic_sept = aug31_col.merge(oct02_col).median().clip(bounds_geom);

  // If the target condition is met, completely substitute the synthetic median
  ls_median = ee.Image(ee.Algorithms.If(isTargetCond, synthetic_sept, ls_median));

  var optical_bands = ls_median.setDefaultProjection({crs: projLandsat.crs(), scale: projLandsat.nominalScale()});

  var ndvi = optical_bands.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var bsi = optical_bands.expression(
    '((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))', {
      'Blue':  optical_bands.select('Blue'),
      'Red':   optical_bands.select('Red'),
      'NIR':   optical_bands.select('NIR'),
      'SWIR1': optical_bands.select('SWIR1')
    }).rename('BSI');
  var nbr2 = optical_bands.normalizedDifference(['SWIR1', 'SWIR2']).rename('NBR2');

  return optical_bands.addBands([
    ndvi, bsi, nbr2,
    terrain_slope.rename('slope'),
    terrain_aspect.rename('aspect')
  ]);
}

print('Building the 30m Landsat export task... Check your Tasks tab.');

// -------------------------------------------------------------------------
// EXTRACT 30M TIME-SERIES FOR EXPORT
// -------------------------------------------------------------------------
var regionalTimeSeriesData = cloud_windows.map(function(window) {
  var valid_ids_str = ee.String(window.get('Valid_IDs'));
  var year = ee.Number(window.get('Year'));
  var month = ee.Number(window.get('Month'));
  
  var ls_img = buildLandsatComposite(valid_ids_str, year, month);

  var p_bgr = ls_img.classify(model_bgr).rename('BGR_pct');
  var p_lpi = ls_img.classify(model_lpi).rename('LPI_pct');
  var p_mft = ls_img.classify(model_mft).rename('Fetch_m');
  var p_herb = ls_img.classify(model_herb).rename('Herb_pct');
  var p_woody = ls_img.classify(model_woody).rename('Woody_pct');
  var p_hwr = ls_img.classify(model_hwr).rename('Log_HWR');

  var combined_preds = ee.Image.cat([p_bgr, p_lpi, p_mft, p_herb, p_woody, p_hwr]);

  // EXACT 30m REDUCTION (Landsat Native Resolution)
  var regional_mean = combined_preds.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: bounds_geom,
    scale: 30,       // Full native resolution for Landsat
    tileScale: 16,   // Max memory chunking
    maxPixels: 1e13
  });

  return ee.Feature(null, {
    'Date': window.get('Start_Date'), // Extracted from the Cloud FeatureClass
    'System_Time': window.get('system:time_start'),
    'Mean_BGR_pct': regional_mean.get('BGR_pct'),
    'Mean_LPI_pct': regional_mean.get('LPI_pct'),
    'Mean_Fetch_m': regional_mean.get('Fetch_m'),
    'Mean_Herb_pct': regional_mean.get('Herb_pct'),
    'Mean_Woody_pct': regional_mean.get('Woody_pct'),
    'Mean_Log_HWR': regional_mean.get('Log_HWR') 
  });
});

var validRegionalData = regionalTimeSeriesData.filter(ee.Filter.notNull(['Mean_BGR_pct']));

// EXPORT DIRECTLY TO DRIVE
Export.table.toDrive({
  collection: validRegionalData,
  description: 'SRER_Regional_Average_TimeSeries_Landsat', // Renamed to prevent overwriting Sentinel-2 output
  folder: 'GEE_Downloads',
  fileFormat: 'CSV',
  selectors: ['Date', 'System_Time', 'Mean_BGR_pct', 'Mean_LPI_pct', 'Mean_Fetch_m', 'Mean_Herb_pct', 'Mean_Woody_pct', 'Mean_Log_HWR']
});
