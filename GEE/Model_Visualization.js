// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// TOPOGRAPHIC PRE-PROCESSING & HILLSHADE VISUALS
// =========================================================================
// Load your newly exported, pre-calculated 10m resampled DEM asset instead of the catalog
var dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled'); 

// Compute hillshade and terrain parameters instantly without on-the-fly aggregation math
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
var hillshade_norm = hillshade.divide(255.0);

var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// =========================================================================
// BACKGROUND MODEL TRAINING
// =========================================================================
var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

// --- Hyperparameters ---
var goldilocks_params = {
  numberOfTrees: 300,
  shrinkage: 0.05,
  samplingRate: 0.5,
  maxNodes: 24,
  seed: 123
};

var model_bgr = ee.Classifier.smileGradientTreeBoost(goldilocks_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(goldilocks_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(goldilocks_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'MFT', inputProperties: inputProps});

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
      // --- Defense 1, 2, 3 ---
      var probMask = img.select('MSK_CLDPRB').lt(20);
      var scl = img.select('SCL');
      var sclMask = scl.neq(8).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(3));
      var blueMask = img.select('B2').lt(2500);
      var masterMask = probMask.and(sclMask).and(blueMask);
      var maskedImg = img.updateMask(masterMask);

      // --- Calculate Illumination Condition (NO DIVISION) ---
      var sz_num = ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180);
      var sa_num = ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180);
            
      var cosZ = ee.Image.constant(sz_num.cos());
      var sinZ = ee.Image.constant(sz_num.sin());
      var sa_img = ee.Image.constant(sa_num);
            
      var cosS = terrain_slope.cos();
      var sinS = terrain_slope.sin();
      var cosAzAsp = sa_img.subtract(terrain_aspect).cos();
            
      // Calculate illumination and append it directly as a band
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

  // --- Stack optical bands, indices, and all 3 topographic predictors ---
  return optical_bands.addBands([
    ndvi, mcari, bsi, nbr2, 
    sent2_median.select('illumination'),
    terrain_slope.rename('slope'),
    terrain_aspect.rename('aspect')
  ]);
}

function drapeHillshade(image, minVal, maxVal) {
  var palette = ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027'];
  var hsv = image.visualize({min: minVal, max: maxVal, palette: palette}).divide(255.0).rgbToHsv();
  
  var draped = ee.Image.cat([
    hsv.select('hue'),
    hsv.select('saturation').multiply(1.2).clamp(0, 1), 
    hsv.select('value').multiply(0.4).add(hillshade_norm.multiply(0.6)) 
  ]);
  return draped.hsvToRgb();
}

// =========================================================================
// REGIONAL TIME-SERIES CHART (PRINT TO CONSOLE)
// =========================================================================

// Map over all optimal cloud windows to generate the regional predictions
var regionalTimeSeriesData = cloud_windows.map(function(window) {
  var sDate = ee.String(window.get('Start_Date'));
  var eDate = ee.Date(sDate).advance(7, 'day');
  var s2_img = buildS2Composite(sDate, eDate);
  var p_bgr = s2_img.classify(model_bgr).rename('BGR');
  var p_lpi = s2_img.classify(model_lpi).rename('LPI');
  var p_mft = s2_img.classify(model_mft).rename('MFT');
  var combined_preds = ee.Image.cat([p_bgr, p_lpi, p_mft]);

  var regional_mean = combined_preds.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: bounds_geom,
    scale: 60, 
    maxPixels: 1e13 
  });

  return ee.Feature(null, {
    'system:time_start': window.get('system:time_start'),
    'BGR_pct': regional_mean.get('BGR'),
    'LPI_pct': regional_mean.get('LPI'),
    'Fetch_m': regional_mean.get('MFT')
  });
});

regionalTimeSeriesData = regionalTimeSeriesData
  .filter(ee.Filter.notNull(['BGR_pct', 'LPI_pct', 'Fetch_m']))
  .sort('system:time_start');

var regionalChart = ui.Chart.feature.byFeature({
  features: regionalTimeSeriesData,
  xProperty: 'system:time_start',
  yProperties: ['BGR_pct', 'LPI_pct', 'Fetch_m']
})
.setChartType('LineChart')
.setOptions({
  title: 'SRER Regional Average: Predicted Metrics Over Time',
  hAxis: {title: 'Date', format: 'MMM yyyy'},
  series: {
    0: {targetAxisIndex: 0, color: '#d73027', lineWidth: 2, pointSize: 3, labelInLegend: 'Mean BGR (%)'},
    1: {targetAxisIndex: 0, color: '#fc8d59', lineWidth: 2, pointSize: 3, labelInLegend: 'Mean LPI (%)'},
    2: {targetAxisIndex: 1, color: '#1a9850', lineWidth: 2, pointSize: 3, labelInLegend: 'Mean Fetch (m)'}
  },
  vAxes: {
    0: {title: 'Percentage (%)'},
    1: {title: 'Mean Fetch Distance (m)'}
  },
  interpolateNulls: true
});

print(regionalChart);

// =========================================================================
// UI WIDGETS & INTERACTIVITY
// =========================================================================
var mainPanel = ui.Panel({style: {width: '400px', padding: '15px', backgroundColor: '#f8f9fa'}});

var title = ui.Label('SRER Predictive Model', {fontWeight: 'bold', fontSize: '20px', margin: '0 0 15px 0', backgroundColor: '#f8f9fa'});
var desc = ui.Label('Select a date to render predictive maps. Click anywhere on the map to generate a time-series chart for that location.', {fontSize: '12px', color: '#555', margin: '0 0 20px 0', backgroundColor: '#f8f9fa'});

var yearLabel = ui.Label('Select Year (2018 - 2025):', {fontWeight: 'bold', backgroundColor: '#f8f9fa'});
var yearSlider = ui.Slider({min: 2018, max: 2025, value: 2019, step: 1, style: {width: '90%'}});
var monthLabel = ui.Label('Select Month (1 - 12):', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var monthSlider = ui.Slider({min: 1, max: 12, value: 5, step: 1, style: {width: '90%'}});
var statusLabel = ui.Label('Ready to render.', {color: 'blue', fontWeight: 'bold', margin: '20px 0', backgroundColor: '#f8f9fa'});
var renderBtn = ui.Button({label: 'Generate Predictive Maps', onClick: updateMap, style: {stretch: 'horizontal', margin: '20px 0 0 0'}});
var chartPanel = ui.Panel({style: {margin: '20px 0 0 0', backgroundColor: '#f8f9fa'}});

mainPanel.add(title); 
mainPanel.add(desc); 
mainPanel.add(yearLabel);
mainPanel.add(yearSlider); 
mainPanel.add(monthLabel);
mainPanel.add(monthSlider); 
mainPanel.add(renderBtn);
mainPanel.add(statusLabel); 
mainPanel.add(chartPanel); 

ui.root.insert(0, mainPanel);

Map.centerObject(bounds_geom, 12); 
Map.setOptions('ROADMAP'); 
Map.style().set('cursor', 'crosshair');

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

  var bgr_draped = drapeHillshade(p_bgr, 0, 75);
  var lpi_draped = drapeHillshade(p_lpi, 0, 80);
  var mft_draped = drapeHillshade(p_mft, 0, 0.5);

  var markerLayer = null;
  Map.layers().forEach(function(layer) {
    if (layer.getName() === 'Selected Point') { markerLayer = layer; }
  });
  Map.layers().reset();

  Map.addLayer(hillshade, {min: 0, max: 255}, 'Terrain Hillshade', false);
  Map.addLayer(s2_img, {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3, gamma: 1.4}, 'Sentinel-2 RGB (' + y + '-' + m + ')', false);
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

  var timeSeriesData = cloud_windows.map(function(window) {
    var sDate = ee.String(window.get('Start_Date'));
    var eDate = ee.Date(sDate).advance(7, 'day');
    var s2_img = buildS2Composite(sDate, eDate);
    var p_bgr = s2_img.classify(model_bgr).rename('BGR');
    var p_lpi = s2_img.classify(model_lpi).rename('LPI');
    var p_mft = s2_img.classify(model_mft).rename('MFT');
    var combined_preds = ee.Image.cat([p_bgr, p_lpi, p_mft]);
    
    var sampled = combined_preds.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: point,
      scale: 10
    });

    return ee.Feature(null, {
      'system:time_start': window.get('system:time_start'),
      'BGR_pct': sampled.get('BGR'),
      'LPI_pct': sampled.get('LPI'),
      'Fetch_m': sampled.get('MFT')
    });
  });

  timeSeriesData = timeSeriesData
    .filter(ee.Filter.notNull(['BGR_pct', 'LPI_pct', 'Fetch_m']))
    .sort('system:time_start');
  
  var chart = ui.Chart.feature.byFeature({
    features: timeSeriesData,
    xProperty: 'system:time_start',
    yProperties: ['BGR_pct', 'LPI_pct', 'Fetch_m']
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'Predicted Metrics Over Time',
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
});

updateMap();
