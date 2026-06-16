// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_S2');
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
// BACKGROUND MODEL TRAINING
// =========================================================================
var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

// 1. CORE DATASET: Keep all valid points for BGR, LPI, MFT, Herb_pct, and Woody_pct
var core_training_fc = fc.filter(ee.Filter.notNull(inputProps.concat(['Herb_pct', 'Woody_pct'])));

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
var model_herb = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'Herb_pct', inputProperties: inputProps});
var model_woody = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: core_training_fc, classProperty: 'Woody_pct', inputProperties: inputProps});
var model_hwr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: hwr_training_fc, classProperty: 'Log_HWR', inputProperties: inputProps});

// =========================================================================
// SENTINEL-2 EXTRACTION
// =========================================================================
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(bounds_geom).first().select('B2').projection();

function buildS2Composite(startDate, endDate) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(bounds_geom)                    
    .filterDate(startDate, endDate)
    .map(function(img) {
      return img.set('has_aux_bands', img.bandNames().contains('MSK_CLDPRB'));
    })
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
    .multiply(0.0001)
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = optical_bands.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mcari = optical_bands.expression('((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
    'B3': optical_bands.select('B3'), 'B4': optical_bands.select('B4'), 'B5': optical_bands.select('B5')  
  }).rename('MCARI');
  var bsi = optical_bands.expression('((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
    'B2': optical_bands.select('B2'), 'B4': optical_bands.select('B4'), 'B8': optical_bands.select('B8'), 'B11': optical_bands.select('B11')
  }).rename('BSI');
  var nbr2 = optical_bands.normalizedDifference(['B11', 'B12']).rename('NBR2');

  return optical_bands.addBands([
    ndvi, mcari, bsi, nbr2, 
    sent2_median.select('illumination'),
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

var title = ui.Label('SRER Predictive Model', {fontWeight: 'bold', fontSize: '20px', margin: '0 0 15px 0', backgroundColor: '#f8f9fa'});
var desc = ui.Label('Select a date to render predictive maps. Click anywhere on the map to extract local time-series data. A time series chart will appear at the bottom of this panel upon clicking.', {fontSize: '12px', color: '#555', margin: '0 0 20px 0', backgroundColor: '#f8f9fa'});

var yearLabel = ui.Label('Select Year (2018 - 2025):', {fontWeight: 'bold', backgroundColor: '#f8f9fa'});
var yearSlider = ui.Slider({min: 2018, max: 2025, value: 2019, step: 1, style: {width: '90%'}});
var monthLabel = ui.Label('Select Month (1 - 12):', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var monthSlider = ui.Slider({min: 1, max: 12, value: 5, step: 1, style: {width: '90%'}});

// Added new Absolute Cover Chart Option to Memory Management
var chartSelectLabel = ui.Label('Select Chart to Render:', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var chartSelect = ui.Select({
  items: ['Core Metrics (BGR, LPI, MFT)', 'Absolute Cover (Herb & Woody)', 'Herb-to-Woody Ratio (Log HWR)'],
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
// MAP RENDERING
// =========================================================================
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  statusLabel.setValue('Searching for optimal imagery window...');

  var window_query = cloud_windows.filter(ee.Filter.eq('Year', y)).filter(ee.Filter.eq('Month', m));
  var optimal_window = ee.Feature(window_query.first());
  var sDate = ee.String(optimal_window.get('Start_Date'));
  var eDate = ee.Date(sDate).advance(7, 'day'); 
  var s2_img = buildS2Composite(sDate, eDate);

  var p_bgr = s2_img.classify(model_bgr).rename('Pred_BGR');
  var p_lpi = s2_img.classify(model_lpi).rename('Pred_LPI');
  var p_mft = s2_img.classify(model_mft).rename('Pred_MFT');
  var p_herb = s2_img.classify(model_herb).rename('Pred_Herb');
  var p_woody = s2_img.classify(model_woody).rename('Pred_Woody');
  var p_hwr = s2_img.classify(model_hwr).rename('Pred_Log_HWR');

  // Custom Palettes for Visual Differentiation
  var standardPalette = ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027'];
  var reversedPalette = ['#d73027', '#fc8d59', '#fee08b', '#ffffbf', '#d9ef8b', '#91cf60', '#1a9850'];
  var herbPalette = ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#006837', '#004529'];
  var woodyPalette = ['#f5f5f5', '#dfc27d', '#bf812d', '#8c510a', '#543005'];

  var bgr_draped = drapeHillshade(p_bgr, 0, 75, standardPalette);
  var lpi_draped = drapeHillshade(p_lpi, 0, 80, standardPalette);
  var mft_draped = drapeHillshade(p_mft, 0, 0.5, standardPalette);
  var herb_draped = drapeHillshade(p_herb, 0, 80, herbPalette);
  var woody_draped = drapeHillshade(p_woody, 0, 60, woodyPalette);
  var hwr_draped = drapeHillshade(p_hwr, 0, 2, reversedPalette);

  var markerLayer = null;
  Map.layers().forEach(function(layer) {
    if (layer.getName() === 'Selected Point') { markerLayer = layer; }
  });
  Map.layers().reset();

  Map.addLayer(hillshade, {min: 0, max: 255}, 'Terrain Hillshade', false);
  Map.addLayer(s2_img, {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3, gamma: 1.4}, 'Sentinel-2 RGB (' + y + '-' + m + ')', false);
  Map.addLayer(hwr_draped, {}, 'Predicted Log HWR (' + y + '-' + m + ')', false);
  Map.addLayer(woody_draped, {}, 'Predicted Woody Cover (' + y + '-' + m + ')', false);
  Map.addLayer(herb_draped, {}, 'Predicted Herb Cover (' + y + '-' + m + ')', false);
  Map.addLayer(mft_draped, {}, 'Predicted MFT (' + y + '-' + m + ')', false);
  Map.addLayer(lpi_draped, {}, 'Predicted LPI (' + y + '-' + m + ')', false);
  Map.addLayer(bgr_draped, {}, 'Predicted BGR (' + y + '-' + m + ')');
  
  var empty = ee.Image().byte();
  var outline = empty.paint({featureCollection: bounds_fc, color: 1, width: 2});
  Map.addLayer(outline, {palette: '#000000'}, 'SRER Boundary');

  if (markerLayer) { Map.layers().add(markerLayer); }
  
  optimal_window.get('Window_Label').evaluate(function(label) {
    if (label) {
      statusLabel.setValue('Displaying Window: \n' + label);
    } else {
      statusLabel.setValue('No valid data found for this date.');
    }
  });
}

// =========================================================================
// INTERACTIVE TIME-SERIES
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

  statusLabel.setValue('Extracting time-series... please wait.');
  var selectedChart = chartSelect.getValue(); 

  var timeSeriesRaw = cloud_windows.map(function(window) {
    var sDate = ee.String(window.get('Start_Date'));
    var eDate = ee.Date(sDate).advance(7, 'day');
    var s2_img = buildS2Composite(sDate, eDate);
    
    var sampled = s2_img.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: point,
      scale: 10
    });

    return ee.Feature(null, sampled).set('system:time_start', window.get('system:time_start'));
  });

  var validTimeSeries = timeSeriesRaw.filter(ee.Filter.notNull(inputProps));

  var classifiedTimeSeries = validTimeSeries
    .classify({classifier: model_bgr, outputName: 'BGR_pct'})
    .classify({classifier: model_lpi, outputName: 'LPI_pct'})
    .classify({classifier: model_mft, outputName: 'Fetch_m'})
    .classify({classifier: model_herb, outputName: 'Herb_pct'})
    .classify({classifier: model_woody, outputName: 'Woody_pct'})
    .classify({classifier: model_hwr, outputName: 'Log_HWR'})
    .sort('system:time_start');
  
  if (selectedChart === 'Core Metrics (BGR, LPI, MFT)') {
    var coreChart = ui.Chart.feature.byFeature({
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
    chartPanel.add(coreChart);
    
  } else if (selectedChart === 'Absolute Cover (Herb & Woody)') {
    var hwCoverChart = ui.Chart.feature.byFeature({
      features: classifiedTimeSeries,
      xProperty: 'system:time_start',
      yProperties: ['Herb_pct', 'Woody_pct'] 
    })
    .setChartType('LineChart')
    .setOptions({
      title: 'Predicted Absolute Cover Over Time',
      hAxis: {title: 'Date', format: 'MMM yyyy'},
      series: {
        0: {color: '#91cf60', lineWidth: 2, pointSize: 3, labelInLegend: 'Herbaceous (%)'},
        1: {color: '#8c510a', lineWidth: 2, pointSize: 3, labelInLegend: 'Woody (%)'}
      },
      vAxis: {title: 'Cover Percentage (%)', viewWindow: {min: 0, max: 100}},
      interpolateNulls: true,
      backgroundColor: '#f8f9fa'
    });
    chartPanel.clear();
    chartPanel.add(hwCoverChart);
    
  } else {
    var hwrChart = ui.Chart.feature.byFeature({
      features: classifiedTimeSeries,
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
