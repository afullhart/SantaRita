// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass');

// Clean extraction of geometry bounds to prevent compiler token errors
var bounds_geom = bounds_fc.first().geometry().bounds();

print(fc.first());

// Global memory variables for the Pixel Inspector
var current_rap_img = null;
var current_srer_img = null;

// =========================================================================
// TOPOGRAPHIC PRE-PROCESSING
// =========================================================================
// Load your pre-calculated 10m resampled DEM asset instead of the catalog
var dem = ee.Image('projects/ee-andrewfullhart/assets/SR_10m_DEM_Resampled');

var terrain_slope = ee.Terrain.slope(dem).multiply(Math.PI / 180);
var terrain_aspect = ee.Terrain.aspect(dem).multiply(Math.PI / 180);

// =========================================================================
// BACKGROUND MODEL TRAINING (ONLY BGR REQUIRED)
// =========================================================================
var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'];

// We kept null features in the asset to align with the shapefile geometry,
// but the classifier requires valid numbers. We drop them here right before training.
var training_fc = fc.filter(ee.Filter.notNull(inputProps));

// --- Hyperparameters ---
var hyperpars = {
  numberOfTrees: 400,
  shrinkage: 0.05,
  samplingRate: 0.7,
  maxNodes: 32,
  loss: 'Huber',
  seed: 123
};

var model_bgr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION')
  .train({features: training_fc, classProperty: 'BGR', inputProperties: inputProps});

// =========================================================================
// SENTINEL-2 EXTRACTION (WITH TOPOGRAPHIC PREDICTORS)
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

// --- Dynamic HSV Draping Function ---
var hillshade = ee.Terrain.hillshade(dem, 270, 45);
var hillshade_norm = hillshade.divide(255.0);

function drapeHillshade(image, minVal, maxVal, customPalette) {
  var hsv = image.visualize({min: minVal, max: maxVal, palette: customPalette}).divide(255.0).rgbToHsv();
  
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
print('Generating regional averages for the console... (This may take ~30 seconds)');

var years = ee.List.sequence(2018, 2025);

var regionalTimeSeriesData = ee.FeatureCollection(years.map(function(y) {
  var year = ee.Number(y);
  
  // 1. RAP Annual Mean for the region (Selecting BGR and LTR)
  var rap_col = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-10m')
    .filterDate(year.format('%d').cat('-01-01'), year.format('%d').cat('-12-31'))
    .filterBounds(bounds_geom);
    
  var rap_annual_yr = rap_col.mosaic().select(['BGR', 'LTR']);
  
  // Create sum band
  var rap_sum = rap_annual_yr.select('BGR').add(rap_annual_yr.select('LTR')).rename('RAP_Sum');
  var rap_combined = ee.Image.cat([rap_annual_yr.select('BGR').rename('RAP_BGR'), rap_sum]);
  
  var rap_mean = rap_combined.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: bounds_geom,
    scale: 60, 
    maxPixels: 1e9
  });

  // 2. SRER Model Annual Mean for the region
  var year_wins = cloud_windows.filter(ee.Filter.eq('Year', year)).toList(12);
  var model_col = ee.ImageCollection(year_wins.map(function(win) {
    var feature = ee.Feature(win);
    var sDate = ee.String(feature.get('Start_Date'));
    var eDate = ee.Date(sDate).advance(7, 'day');
    return buildS2Composite(sDate, eDate).classify(model_bgr).rename('Model_BGR');
  }));
  var srer_annual_yr = model_col.mean();
  
  var srer_mean = srer_annual_yr.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: bounds_geom,
    scale: 60, 
    maxPixels: 1e9
  }).get('Model_BGR');

  return ee.Feature(null, {
    'Year': year.format('%d'),
    'RAP_BGR': rap_mean.get('RAP_BGR'),
    'RAP_Sum': rap_mean.get('RAP_Sum'),
    'Model_BGR': srer_mean
  });
})).filter(ee.Filter.notNull(['RAP_BGR', 'Model_BGR', 'RAP_Sum']));

var regionalChart = ui.Chart.feature.byFeature({
  features: regionalTimeSeriesData,
  xProperty: 'Year',
  yProperties: ['RAP_BGR', 'Model_BGR', 'RAP_Sum']
})
.setChartType('LineChart')
.setOptions({
  title: 'SRER Regional Average: Bare Ground % (2018-2025)',
  hAxis: {title: 'Year'},
  vAxis: {title: 'Mean Bare Ground (%)'},
  series: {
    0: {color: '#1a9850', lineWidth: 2, pointSize: 4, labelInLegend: 'RAP Bare Ground (10m)'}, 
    1: {color: '#d73027', lineWidth: 2, pointSize: 4, labelInLegend: 'SRER Model'},
    2: {color: '#984ea3', lineWidth: 2, pointSize: 4, labelInLegend: 'RAP BGR + LTR Sum', lineDashStyle: [4, 4]} 
  },
  interpolateNulls: true
});

print(regionalChart);

// =========================================================================
// UI WIDGETS & INTERACTIVITY
// =========================================================================
var mainPanel = ui.Panel({style: {width: '400px', padding: '15px', backgroundColor: '#f8f9fa'}});

var title = ui.Label('SRER Model vs. RAP (10m)', {fontWeight: 'bold', fontSize: '20px', margin: '0 0 10px 0', backgroundColor: '#f8f9fa'});
var desc = ui.Label('Aggregates 12 monthly predictions from the SRER model and compares the annual mean against the 10m Rangeland Analysis Platform (RAP). Click anywhere on the map to inspect raw pixel values and generate a time-series chart.', {fontSize: '12px', color: '#555', margin: '0 0 20px 0', backgroundColor: '#f8f9fa'});

var yearLabel = ui.Label('Select Year:', {fontWeight: 'bold', backgroundColor: '#f8f9fa'});
var yearSlider = ui.Slider({min: 2018, max: 2025, value: 2019, step: 1, style: {width: '90%'}}); 
var statusLabel = ui.Label('Ready to process.', {color: 'blue', fontWeight: 'bold', margin: '20px 0', backgroundColor: '#f8f9fa'});
var renderBtn = ui.Button({label: 'Compare RAP vs. Model', onClick: updateMap, style: {stretch: 'horizontal', margin: '10px 0'}});

// Pixel Inspector Panel
var inspectorPanel = ui.Panel({style: {margin: '20px 0 0 0', padding: '10px', border: '1px solid #ccc', backgroundColor: '#ffffff'}});
var inspectorTitle = ui.Label('Pixel Inspector', {fontWeight: 'bold', margin: '0 0 5px 0', backgroundColor: '#ffffff'});
var inspectorCoords = ui.Label('Click on the map to inspect...', {fontSize: '12px', color: '#777', margin: '0 0 10px 0', backgroundColor: '#ffffff'});
var rapLabel = ui.Label('RAP (10m): --', {backgroundColor: '#ffffff'});
var srerLabel = ui.Label('SRER Model: --', {backgroundColor: '#ffffff'});
var diffLabel = ui.Label('Difference: --', {backgroundColor: '#ffffff'});

inspectorPanel.add(inspectorTitle).add(inspectorCoords).add(rapLabel).add(srerLabel).add(diffLabel);

// Chart Panel
var chartPanel = ui.Panel({style: {margin: '15px 0 0 0', backgroundColor: '#f8f9fa'}});

// Legend Definitions
var bgrPalette = ['#1a9850', '#91cf60', '#d9ef8b', '#ffffbf', '#fee08b', '#fc8d59', '#d73027'];
var diffPalette = ['#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020']; 

mainPanel.add(title).add(desc).add(yearLabel).add(yearSlider).add(renderBtn).add(statusLabel).add(inspectorPanel).add(chartPanel);
ui.root.insert(0, mainPanel);

Map.centerObject(bounds_geom, 12);
Map.setOptions('HYBRID'); 
Map.style().set('cursor', 'crosshair'); 

// =========================================================================
// RENDER FUNCTION
// =========================================================================
function updateMap() {
  var y = yearSlider.getValue();
  statusLabel.setValue('Processing 12 monthly composites...\n(This will take ~10-15 seconds)');
  statusLabel.style().set('color', 'orange');

  var rap_collection = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-10m')
    .filterDate(y + '-01-01', y + '-12-31')
    .filterBounds(bounds_geom);
  
  var rap_annual = rap_collection.mosaic().select('BGR').clip(bounds_geom);

  var year_windows_list = cloud_windows.filter(ee.Filter.eq('Year', y)).toList(12);
  
  var monthly_predictions = ee.ImageCollection(year_windows_list.map(function(win) {
    var feature = ee.Feature(win);
    var sDate = ee.String(feature.get('Start_Date'));
    var eDate = ee.Date(sDate).advance(7, 'day');
    
    var s2_img = buildS2Composite(sDate, eDate);
    return s2_img.classify(model_bgr).rename('Model_BGR');
  }));
  
  var srer_annual = monthly_predictions.mean().clip(bounds_geom);
  var difference = srer_annual.subtract(rap_annual).rename('Diff_BGR');

  // --- Save raw images to global memory for the Pixel Inspector ---
  current_rap_img = rap_annual;
  current_srer_img = srer_annual;

  // Apply HSV Terrain Draping with 0-80 scaling
  var rap_draped = drapeHillshade(rap_annual, 0, 80, bgrPalette);
  var srer_draped = drapeHillshade(srer_annual, 0, 80, bgrPalette);
  var diff_draped = drapeHillshade(difference, -20, 20, diffPalette);

  // Preserve the clicked point marker if it exists
  var markerLayer = null;
  Map.layers().forEach(function(layer) {
    if (layer.getName() === 'Selected Point') { markerLayer = layer; }
  });

  Map.layers().reset();
  
  var empty = ee.Image().byte();
  var outline = empty.paint({featureCollection: bounds_fc, color: 1, width: 2});
  Map.addLayer(outline, {palette: '#000000'}, 'SRER Boundary');
  
  Map.addLayer(hillshade, {min: 0, max: 255}, 'Terrain Hillshade', false);
  Map.addLayer(rap_draped, {}, 'RAP 10m Bare Ground (' + y + ')', false);
  Map.addLayer(srer_draped, {}, 'SRER Model Bare Ground (' + y + ')', false);
  Map.addLayer(diff_draped, {}, 'Difference (Model minus RAP)');

  if (markerLayer) { Map.layers().add(markerLayer); }

  statusLabel.setValue('SUCCESS!\nDisplaying Difference Map for ' + y);
  statusLabel.style().set('color', 'green');
  
  // Reset inspector labels for the new year
  inspectorCoords.setValue('Click on the map to inspect...');
  rapLabel.setValue('RAP (10m): --');
  srerLabel.setValue('SRER Model: --');
  diffLabel.setValue('Difference: --');
}

// =========================================================================
// PIXEL INSPECTOR & CHART GENERATOR (ON CLICK)
// =========================================================================
Map.onClick(function(coords) {
  // Prevent inspection if layers haven't been rendered yet
  if (!current_rap_img || !current_srer_img) {
    inspectorCoords.setValue('Please generate the maps first.');
    return;
  }

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

  // Update UI to show loading state
  inspectorCoords.setValue('Loading point data... (' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4) + ')');
  rapLabel.setValue('RAP (10m): --');
  srerLabel.setValue('SRER Model: --');
  diffLabel.setValue('Difference: --');
  
  chartPanel.clear();
  chartPanel.add(ui.Label('Generating 2018-2025 time-series chart...\n(This may take a few seconds)', {color: '#777', fontSize: '12px', whiteSpace: 'pre-wrap', backgroundColor: '#f8f9fa'}));

  // 1. EVALUATE SINGLE YEAR FOR THE INSPECTOR PANEL
  var combined = ee.Image.cat([
    current_rap_img.rename('RAP'), 
    current_srer_img.rename('SRER')
  ]);
  
  var sampled = combined.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  });

  sampled.evaluate(function(result) {
    if (result && result.RAP !== null && result.SRER !== null) {
      inspectorCoords.setValue('Lon: ' + coords.lon.toFixed(4) + ', Lat: ' + coords.lat.toFixed(4));
      rapLabel.setValue('RAP (10m): ' + result.RAP.toFixed(2) + '%');
      srerLabel.setValue('SRER Model: ' + result.SRER.toFixed(2) + '%');
      
      var diff = result.SRER - result.RAP;
      var diffText = diff > 0 ? '+' + diff.toFixed(2) + '%' : diff.toFixed(2) + '%';
      diffLabel.setValue('Difference: ' + diffText);
    } else {
      inspectorCoords.setValue('No data at this location.');
      rapLabel.setValue('RAP (10m): N/A');
      srerLabel.setValue('SRER Model: N/A');
      diffLabel.setValue('Difference: N/A');
    }
  });

  // 2. GENERATE TIME-SERIES CHART ACROSS ALL YEARS
  var years = ee.List.sequence(2018, 2025);
  
  var timeSeriesData = ee.FeatureCollection(years.map(function(y) {
    var year = ee.Number(y);
    
    // RAP Annual for this year (Selecting BGR and LTR)
    var rap_col = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-10m')
      .filterDate(year.format('%d').cat('-01-01'), year.format('%d').cat('-12-31'))
      .filterBounds(point);
      
    var rap_annual_yr = rap_col.mosaic().select(['BGR', 'LTR']);
    var rap_sum = rap_annual_yr.select('BGR').add(rap_annual_yr.select('LTR')).rename('RAP_Sum');
    var rap_combined_yr = ee.Image.cat([rap_annual_yr.select('BGR').rename('RAP_BGR'), rap_sum]);

    // SRER Model Annual for this year
    var year_wins = cloud_windows.filter(ee.Filter.eq('Year', year)).toList(12);
    var model_col = ee.ImageCollection(year_wins.map(function(win) {
      var feature = ee.Feature(win);
      var sDate = ee.String(feature.get('Start_Date'));
      var eDate = ee.Date(sDate).advance(7, 'day');
      return buildS2Composite(sDate, eDate).classify(model_bgr).rename('Model_BGR');
    }));
    var srer_annual_yr = model_col.mean();

    // Sample the point
    var combined_yr = ee.Image.cat([rap_combined_yr, srer_annual_yr.rename('Model_BGR')]);
    var sampled_yr = combined_yr.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: point,
      scale: 10
    });

    return ee.Feature(null, {
      'Year': year.format('%d'),
      'RAP_BGR': sampled_yr.get('RAP_BGR'),
      'RAP_Sum': sampled_yr.get('RAP_Sum'),
      'Model_BGR': sampled_yr.get('Model_BGR')
    });
  })).filter(ee.Filter.notNull(['RAP_BGR', 'Model_BGR', 'RAP_Sum']));

  // Render the Chart
  var chart = ui.Chart.feature.byFeature({
    features: timeSeriesData,
    xProperty: 'Year',
    yProperties: ['RAP_BGR', 'Model_BGR', 'RAP_Sum']
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'Point Trend: Bare Ground % (2018-2025)',
    hAxis: {title: 'Year'},
    vAxis: {title: 'Bare Ground (%)'},
    series: {
      0: {color: '#1a9850', lineWidth: 2, pointSize: 4, labelInLegend: 'RAP Bare Ground (10m)'}, 
      1: {color: '#d73027', lineWidth: 2, pointSize: 4, labelInLegend: 'SRER Model'},
      2: {color: '#984ea3', lineWidth: 2, pointSize: 4, labelInLegend: 'RAP BGR + LTR Sum', lineDashStyle: [4, 4]}
    },
    interpolateNulls: true,
    backgroundColor: '#ffffff'
  });

  chartPanel.clear();
  chartPanel.add(chart);
});

// Trigger initial load
updateMap();
