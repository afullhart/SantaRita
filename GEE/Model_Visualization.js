// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_landsat_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// TOPOGRAPHIC PRE-PROCESSING & HILLSHADE VISUALS
// =========================================================================
var dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled'); 

var hillshade = ee.Terrain.hillshade(dem, 270, 45);
var hillshade_norm = hillshade.divide(255.0);

var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// =========================================================================
// BACKGROUND MODEL TRAINING (LANDSAT PREDICTORS)
// =========================================================================
var inputProps = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

// 1. CORE DATASET: Keep all valid points for BGR, LPI, MFT
var core_training_fc = fc.filter(ee.Filter.notNull(inputProps));

// 2. HWR DATASET: Drop nulls and apply Log Transform specifically for HWR
var hwr_training_fc = core_training_fc.filter(ee.Filter.notNull(['Herb_Woody_Ratio'])).map(function(ft) {
  var log_hwr = ee.Number(ft.get('Herb_Woody_Ratio')).add(1).log(); 
  return ft.set('Log_HWR', log_hwr);
});

var hyperpars = {
  numberOfTrees: 400,
  shrinkage: 0.05,
  samplingRate: 0.7,
  maxNodes: 32,
  loss: 'Huber',
  seed: 123
};

var model_bgr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'MFT', inputProperties: inputProps});
var model_hwr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: hwr_training_fc, classProperty: 'Log_HWR', inputProperties: inputProps});

// =========================================================================
// LANDSAT EXTRACTION PIPELINE & COMPOSITE GENERATOR
// =========================================================================
var projLandsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(bounds_geom).first().select('SR_B2').projection();

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

  // =====================================================================
  // SERVER-SAFE SPECIAL CASE PATCH FOR SEPTEMBER 2019
  // =====================================================================
  // Generates a boolean image (1 if target date, 0 otherwise)
  var isTarget = ee.Image.constant(ee.Number(year).eq(2019).and(ee.Number(month).eq(9)));
  
  var aug31_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(bounds_geom)
    .filterDate('2019-08-31', '2019-09-01')
    .map(prepOLI); 
    
  // If not Sept 2019, the mask zeros out the patch and nothing is filled
  var aug31_patch = aug31_col.median().clip(bounds_geom).updateMask(isTarget);
  ls_median = ls_median.unmask(aug31_patch);
  // =====================================================================

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

function drapeHillshade(image, minVal, maxVal, customPalette) {
  var palette = customPalette || ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027'];
  var hsv = image.visualize({min: minVal, max: maxVal, palette: palette}).divide(255.0).rgbToHsv();
  
  var draped = ee.Image.cat([
    hsv.select('hue'),
    hsv.select('saturation').multiply(1.2).clamp(0, 1), 
    hsv.select('value').multiply(0.4).add(hillshade_norm.multiply(0.6)) 
  ]);
  return draped.hsvToRgb();
}

// =========================================================================
// UI WIDGETS & INTERACTIVITY
// =========================================================================
var mainPanel = ui.Panel({style: {width: '400px', padding: '15px', backgroundColor: '#f8f9fa'}});

var title = ui.Label('SRER Predictive Model (Landsat 1984-2025)', {fontWeight: 'bold', fontSize: '20px', margin: '0 0 15px 0', backgroundColor: '#f8f9fa'});
var desc = ui.Label('Select a date to render predictive maps. Click anywhere on the map to extract local time-series data. A time series chart will appear at the bottom of this panel upon clicking.', {fontSize: '12px', color: '#555', margin: '0 0 20px 0', backgroundColor: '#f8f9fa'});

var yearLabel = ui.Label('Select Year (1984 - 2025):', {fontWeight: 'bold', backgroundColor: '#f8f9fa'});
var yearSlider = ui.Slider({min: 1984, max: 2025, value: 2019, step: 1, style: {width: '90%'}});
var monthLabel = ui.Label('Select Month (1 - 12):', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var monthSlider = ui.Slider({min: 1, max: 12, value: 9, step: 1, style: {width: '90%'}});

// Memory Management Chart Selector
var chartSelectLabel = ui.Label('Select Chart to Render (Saves Memory):', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var chartSelect = ui.Select({
  items: ['Core Metrics (BGR, LPI, MFT)', 'Herb-to-Woody Ratio (Log HWR)'],
  value: 'Core Metrics (BGR, LPI, MFT)',
  style: {width: '90%'}
});

var statusLabel = ui.Label('Ready to render.', {color: 'blue', fontWeight: 'bold', margin: '20px 0', backgroundColor: '#f8f9fa'});
var renderBtn = ui.Button({label: 'Generate Predictive Maps', onClick: updateMap, style: {stretch: 'horizontal', margin: '20px 0 0 0'}});
var chartPanel = ui.Panel({style: {margin: '20px 0 0 0', backgroundColor: '#f8f9fa'}});

mainPanel.add(title).add(desc).add(yearLabel).add(yearSlider).add(monthLabel).add(monthSlider);
mainPanel.add(chartSelectLabel).add(chartSelect); 
mainPanel.add(renderBtn).add(statusLabel).add(chartPanel);

ui.root.insert(0, mainPanel);

Map.centerObject(bounds_geom, 12); 
Map.setOptions('ROADMAP'); 
Map.style().set('cursor', 'crosshair');

// =========================================================================
// RENDER MAP LAYERS
// =========================================================================
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  statusLabel.setValue('Searching for optimal imagery window...');

  var window_query = cloud_windows.filter(ee.Filter.eq('Year', y)).filter(ee.Filter.eq('Month', m));
  
  window_query.size().evaluate(function(count) {
    if (count === 0) {
      statusLabel.setValue('No valid data found for ' + y + '-' + m);
      return;
    }
    
    var optimal_window = ee.Feature(window_query.first());
    
    optimal_window.evaluate(function(feat) {
      var valid_ids_str = feat.properties.Valid_IDs;
      var label = feat.properties.Window_Label;
      
      // Pass the year and month so the September 2019 patch can trigger if needed
      var ls_img = buildLandsatComposite(valid_ids_str, y, m);

      var p_bgr = ls_img.classify(model_bgr).rename('Pred_BGR');
      var p_lpi = ls_img.classify(model_lpi).rename('Pred_LPI');
      var p_mft = ls_img.classify(model_mft).rename('Pred_MFT');
      var p_hwr = ls_img.classify(model_hwr).rename('Pred_Log_HWR');

      var standardPalette = ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027'];
      var reversedPalette = ['#d73027', '#fc8d59', '#fee08b', '#ffffbf', '#d9ef8b', '#91cf60', '#1a9850'];

      var bgr_draped = drapeHillshade(p_bgr, 0, 75, standardPalette);
      var lpi_draped = drapeHillshade(p_lpi, 0, 80, standardPalette);
      var mft_draped = drapeHillshade(p_mft, 0, 0.5, standardPalette);
      var hwr_draped = drapeHillshade(p_hwr, 0, 1.5, reversedPalette);

      var markerLayer = null;
      Map.layers().forEach(function(layer) {
        if (layer.getName() === 'Selected Point') { markerLayer = layer; }
      });
      Map.layers().reset();

      Map.addLayer(hillshade, {min: 0, max: 255}, 'Terrain Hillshade', false);
      Map.addLayer(ls_img, {bands: ['Red', 'Green', 'Blue'], min: 0.0, max: 0.3, gamma: 1.4}, 'Landsat RGB (' + y + '-' + m + ')', false);
      Map.addLayer(hwr_draped, {}, 'Predicted Log HWR (' + y + '-' + m + ')', false);
      Map.addLayer(mft_draped, {}, 'Predicted MFT (' + y + '-' + m + ')', false);
      Map.addLayer(lpi_draped, {}, 'Predicted LPI (' + y + '-' + m + ')', false);
      Map.addLayer(bgr_draped, {}, 'Predicted BGR (' + y + '-' + m + ')');
      
      var empty = ee.Image().byte();
      var outline = empty.paint({featureCollection: bounds_fc, color: 1, width: 2});
      Map.addLayer(outline, {palette: '#000000'}, 'SRER Boundary');

      if (markerLayer) { Map.layers().add(markerLayer); }
      
      statusLabel.setValue('Displaying Window: \n' + label);
    });
  });
}

// =========================================================================
// RENDER DYNAMIC TIME SERIES CHARTS
// =========================================================================
Map.onClick(function(coords) {
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var dot = ui.Map.Layer(point, {color: 'FF0000'}, 'Selected Point');
  
  var layers = Map.layers();
  for (var i = 0; i < layers.length(); i++) {
    if (layers.get(i).getName() === 'Selected Point') {
      layers.remove(layers.get(i));
      break;
    }
  }
  Map.layers().add(dot);

  statusLabel.setValue('Extracting 40-year time-series... please wait.');
  var selectedChart = chartSelect.getValue(); 

  var timeSeriesRaw = cloud_windows.map(function(window) {
    var valid_ids_str = ee.String(window.get('Valid_IDs'));
    var year = ee.Number(window.get('Year'));
    var month = ee.Number(window.get('Month'));
    
    // Pass the server-side year and month to ensure the Sept 2019 data is patched for the chart
    var ls_img = buildLandsatComposite(valid_ids_str, year, month);
    
    // Extacts an exact 60m regional average around the clicked point
    var sampled = ls_img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.buffer(30), 
      scale: 30
    });

    return ee.Feature(null, sampled).set('system:time_start', window.get('system:time_start'));
  });

  var validTimeSeries = timeSeriesRaw.filter(ee.Filter.notNull(inputProps));

  if (selectedChart === 'Core Metrics (BGR, LPI, MFT)') {
    var classifiedTimeSeries = validTimeSeries
      .classify({classifier: model_bgr, outputName: 'BGR_pct'})
      .classify({classifier: model_lpi, outputName: 'LPI_pct'})
      .classify({classifier: model_mft, outputName: 'Fetch_m'})
      .sort('system:time_start');
      
    var chart = ui.Chart.feature.byFeature({
      features: classifiedTimeSeries,
      xProperty: 'system:time_start',
      yProperties: ['BGR_pct', 'LPI_pct', 'Fetch_m'] 
    })
    .setChartType('LineChart')
    .setOptions({
      title: 'Predicted Core Metrics Over Time',
      hAxis: {title: 'Date', format: 'MMM yyyy'},
      series: {
        0: {targetAxisIndex: 0, color: '#d73027', lineWidth: 2, pointSize: 3, labelInLegend: 'BGR (%)'},
        1: {targetAxisIndex: 0, color: '#fc8d59', lineWidth: 2, pointSize: 3, labelInLegend: 'LPI (%)'},
        2: {targetAxisIndex: 1, color: '#1a9850', lineWidth: 2, pointSize: 3, labelInLegend: 'Mean Fetch (m)'}
      },
      vAxes: {
        0: {title: 'Percentage (%)'},
        1: {title: 'Mean Fetch Distance (m)'} 
      },
      interpolateNulls: true,
      backgroundColor: '#f8f9fa'
    });

    chartPanel.clear();
    chartPanel.add(chart);
    
  } else {
    var hwrTimeSeries = validTimeSeries
      .classify({classifier: model_hwr, outputName: 'Log_HWR'})
      .sort('system:time_start');
      
    var hwrChart = ui.Chart.feature.byFeature({
      features: hwrTimeSeries,
      xProperty: 'system:time_start',
      yProperties: ['Log_HWR']
    })
    .setChartType('LineChart')
    .setOptions({
      title: 'Predicted Herb-to-Woody Ratio',
      hAxis: {title: 'Date', format: 'MMM yyyy'},
      series: {
        0: {color: '#4575b4', lineWidth: 2, pointSize: 3, labelInLegend: 'Log HWR'} 
      },
      vAxis: {title: 'Log Ratio'},
      interpolateNulls: true,
      backgroundColor: '#f8f9fa'
    });

    chartPanel.clear();
    chartPanel.add(hwrChart);
  }

  statusLabel.setValue('Time-series chart loaded successfully.');
});

updateMap();
