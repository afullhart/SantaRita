// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

// LOAD THE PRE-CALCULATED S2 CLOUD ASSET
var cloud_windows = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_S2');

// The exact visualization parameters you requested
var rgbVis = {
  bands: ['B4', 'B3', 'B2'],
  min: 0.0,
  max: 0.3, 
  gamma: 1.4
};

// =========================================================================
// CONSOLE CHART: TIME-SERIES OF ALL OPTIMAL WINDOWS (INSTANT RENDER)
// =========================================================================
// Because we are reading an asset, this chart now renders instantly
var sorted_windows = cloud_windows.sort('system:time_start');

var cloudChart = ui.Chart.feature.byFeature({
  features: sorted_windows,
  xProperty: 'system:time_start',
  yProperties: ['Mean_Cloud_Prob']
})
.setChartType('LineChart')
.setOptions({
  title: 'Optimal Window Cloud Probability (2018 - 2025)',
  hAxis: {title: 'Date', format: 'MMM yyyy'},
  vAxis: {title: 'Mean Cloud Probability (%)', viewWindow: {min: 0, max: 100}},
  colors: ['#1f77b4'],
  lineWidth: 2,
  pointSize: 3
});

// Print directly to the console
print('Time-Series: Optimal Window Cloud Probabilities', cloudChart);


// =========================================================================
// UI APP & INTERACTIVITY
// =========================================================================
// Create a main panel to hold the UI
var panel = ui.Panel({
  style: {width: '350px', padding: '15px'}
});
ui.root.insert(0, panel);

// UI Elements
var title = ui.Label('Optimal S2 Window Explorer', {fontWeight: 'bold', fontSize: '20px'});
var desc = ui.Label('Select a Year and Month. The script will fetch the pre-calculated optimal 7-day window from the asset and display it.');

var yearLabel = ui.Label('Select Year:', {fontWeight: 'bold'});
var yearSlider = ui.Slider({
  min: 2018, max: 2025, value: 2019, step: 1,
  style: {stretch: 'horizontal'}
});

var monthLabel = ui.Label('Select Month:', {fontWeight: 'bold'});
var monthSlider = ui.Slider({
  min: 1, max: 12, value: 5, step: 1,
  style: {stretch: 'horizontal'}
});

var statusBox = ui.Label({
  value: 'Ready. Move sliders to load imagery...',
  style: {color: 'blue', margin: '20px 0', whiteSpace: 'pre-wrap'}
});

// Add elements to the panel
panel.add(title);
panel.add(desc);
panel.add(yearLabel).add(yearSlider);
panel.add(monthLabel).add(monthSlider);
panel.add(statusBox);

// Center map on the study area
Map.centerObject(bounds_geom, 13);
Map.setOptions('ROADMAP'); 

// =========================================================================
// RENDER FUNCTION (ASSET LOOKUP)
// =========================================================================
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  
  statusBox.setValue('Querying asset index for optimal window...');
  statusBox.style().set('color', 'orange');
  
  // 1. Filter the asset directly to the chosen year and month
  var window_query = cloud_windows.filter(ee.Filter.eq('Year', y)).filter(ee.Filter.eq('Month', m));
  
  // 2. Check if a valid record exists
  window_query.size().evaluate(function(count) {
    if (count === 0) {
      statusBox.setValue('No record found in the asset for ' + y + '-' + m + '.');
      statusBox.style().set('color', 'red');
      Map.layers().reset(); 
      Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
      return;
    }
    
    // 3. Extract the single best window feature
    var bestWindowFeature = ee.Feature(window_query.first());
    
    bestWindowFeature.evaluate(function(feature) {
      var props = feature.properties;
      
      if (props.Image_Count === 0) {
        statusBox.setValue('No imagery found for ' + y + '-' + m + '.\n(Try a different date)');
        statusBox.style().set('color', 'red');
        Map.layers().reset(); 
        Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
        return;
      }
      
      // Update the UI Text using data straight from the asset columns
      var statusText = 'SUCCESS!' + 
                       '\nOptimal Dates: ' + props.Window_Label + 
                       '\nCloud Probability: ' + props.Mean_Cloud_Prob.toFixed(2) + '%' +
                       '\nImages in Window: ' + props.Image_Count;
      statusBox.setValue(statusText);
      statusBox.style().set('color', 'green');
      
      // 4. Reconstruct the image collection using the specific dates found in the asset
      var startDate = ee.Date(props.Start_Date);
      var endDateFilter = startDate.advance(7, 'day');
      
      var s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(bounds_geom)
        .filterDate(startDate, endDateFilter);
        
      // Triple Defense Masking
      var s2_mosaic = s2_collection
        .map(function(img) {
          // Defense 1: Cloud Probability < 20%
          var probMask = img.select('MSK_CLDPRB').lt(20);
          
          // Defense 2: SCL Explicit Rejection
          var scl = img.select('SCL');
          var sclMask = scl.neq(8).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)).and(scl.neq(3));
          
          // Defense 3: Hard Physical Brightness Limit
          var blueMask = img.select('B2').lt(2500);
          
          // Combine all three
          var masterMask = probMask.and(sclMask).and(blueMask);
          return img.updateMask(masterMask);
        })
        .median()
        .clip(v_extent)
        .multiply(0.0001); 
        
      // Update the Map Layers
      Map.layers().reset(); 
      
      Map.layers().set(0, ui.Map.Layer(s2_mosaic, rgbVis, 'Optimal S2: ' + props.Window_Label));
      Map.layers().set(1, ui.Map.Layer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds'));
    });
  });
}

// Attach the update function to the sliders
yearSlider.onChange(updateMap);
monthSlider.onChange(updateMap);

// Run it once on script start to load the default
updateMap();
