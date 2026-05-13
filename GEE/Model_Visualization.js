// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass');
var bounds_geom = bounds_fc.first().geometry().bounds();

var dem = ee.Image('USGS/3DEP/10m').clip(bounds_geom);
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
var hillshade_norm = hillshade.divide(255.0);

// =========================================================================
// BACKGROUND MODEL TRAINING
// =========================================================================
var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2'];
var regularized_params = {numberOfTrees: 300, shrinkage: 0.01, samplingRate: 0.7, maxNodes: 12, seed: 123};

var model_bgr = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'MFT', inputProperties: inputProps});

// =========================================================================
// SENTINEL-2 EXTRACTION (OPTION A: MEDIAN COMPOSITE)
// =========================================================================
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(bounds_geom).first().select('B2').projection();

function buildS2Composite(startDate, endDate) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(bounds_geom)                     
    .filterDate(startDate, endDate)
    // OPTION A: Mask clouds > 20% probability
    .map(function(img) {
      var cloudMask = img.select('MSK_CLDPRB').lt(20);
      return img.updateMask(cloudMask);
    });                

  var sent2_im = sent2_ic
    // OPTION A: Median reducer to blend clear pixels
    .median() 
    .clip(bounds_geom)
    .select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12']) 
    .multiply(0.0001) 
    .setDefaultProjection({crs: projSent2.crs(), scale: projSent2.nominalScale()});

  var ndvi = sent2_im.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mcari = sent2_im.expression('((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / B4)', {
    'B3': sent2_im.select('B3'), 'B4': sent2_im.select('B4'), 'B5': sent2_im.select('B5')  
  }).rename('MCARI');
  var bsi = sent2_im.expression('((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
    'B2': sent2_im.select('B2'), 'B4': sent2_im.select('B4'), 'B8': sent2_im.select('B8'), 'B11': sent2_im.select('B11')
  }).rename('BSI');
  var nbr2 = sent2_im.normalizedDifference(['B11', 'B12']).rename('NBR2');

  return sent2_im.addBands([ndvi, mcari, bsi, nbr2]);
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
// UI WIDGETS & INTERACTIVITY
// =========================================================================
var mainPanel = ui.Panel({style: {width: '350px', padding: '15px', backgroundColor: '#f8f9fa'}});

var title = ui.Label('SRER Predictive Model', {fontWeight: 'bold', fontSize: '20px', margin: '0 0 15px 0', backgroundColor: '#f8f9fa'});
var desc = ui.Label('Select a year and month. The model will extract the optimal 7-day cloud-free window for that month and generate predictive maps.', {fontSize: '12px', color: '#555', margin: '0 0 20px 0', backgroundColor: '#f8f9fa'});

var yearLabel = ui.Label('Select Year (2018 - 2025):', {fontWeight: 'bold', backgroundColor: '#f8f9fa'});
var yearSlider = ui.Slider({min: 2018, max: 2025, value: 2019, step: 1, style: {width: '90%'}});
var monthLabel = ui.Label('Select Month (1 - 12):', {fontWeight: 'bold', backgroundColor: '#f8f9fa', margin: '15px 0 0 0'});
var monthSlider = ui.Slider({min: 1, max: 12, value: 5, step: 1, style: {width: '90%'}});
var statusLabel = ui.Label('Ready to render.', {color: 'blue', fontWeight: 'bold', margin: '20px 0', backgroundColor: '#f8f9fa'});

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

  var bgr_draped = drapeHillshade(p_bgr, 0, 80);
  var lpi_draped = drapeHillshade(p_lpi, 0, 60);
  var mft_draped = drapeHillshade(p_mft, 0, 0.5);

  Map.layers().reset();
  Map.addLayer(hillshade, {min: 0, max: 255}, 'Terrain Hillshade', false);
  Map.addLayer(s2_img, {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3, gamma: 1.4}, 'Sentinel-2 RGB (' + y + '-' + m + ')', false);
  Map.addLayer(mft_draped, {}, 'Predicted MFT (' + y + '-' + m + ')', false);
  Map.addLayer(lpi_draped, {}, 'Predicted LPI (' + y + '-' + m + ')', false);
  Map.addLayer(bgr_draped, {}, 'Predicted BGR (' + y + '-' + m + ')');
  
  var empty = ee.Image().byte();
  var outline = empty.paint({featureCollection: bounds_fc, color: 1, width: 2});
  Map.addLayer(outline, {palette: '#000000'}, 'SRER Boundary');

  optimal_window.get('Window_Label').evaluate(function(label) {
    if (label) {
      statusLabel.setValue('Displaying Window: \n' + label);
    } else {
      statusLabel.setValue('No valid data found for this date.');
    }
  });
}

var renderBtn = ui.Button({label: 'Generate Predictive Maps', onClick: updateMap, style: {stretch: 'horizontal', margin: '20px 0 0 0'}});

mainPanel.add(title);
mainPanel.add(desc);
mainPanel.add(yearLabel);
mainPanel.add(yearSlider);
mainPanel.add(monthLabel);
mainPanel.add(monthSlider);
mainPanel.add(renderBtn);
mainPanel.add(statusLabel);

ui.root.insert(0, mainPanel);
Map.centerObject(bounds_geom, 12);
Map.setOptions('ROADMAP'); 
updateMap();
