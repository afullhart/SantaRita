// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// DATA PREPARATION
// =========================================================================

// 1. Separate the drone data by season
var fc_may = fc.filter(ee.Filter.eq('Month', 'May'));
var fc_sept = fc.filter(ee.Filter.eq('Month', 'Sept'));

// 2. Extract RAP 2019
var rap_col = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-10m')
  .filterDate('2019-01-01', '2019-12-31')
  .filterBounds(bounds_geom);

var rap_img = rap_col.mosaic().select('BGR');

// 3. THE FIX: Safely sample RAP and force the output name
var sampled_rap = rap_img.reduceRegions({
  collection: fc,
  // Force the output column name so it doesn't overwrite your existing 'mean' column
  reducer: ee.Reducer.mean().setOutputs(['RAP_BGR']),
  scale: 10,
  tileScale: 4
}).filter(ee.Filter.notNull(['RAP_BGR']));

// =========================================================================
// CUMULATIVE DISTRIBUTION (CDF) GENERATOR WITH FAILSAFES
// =========================================================================

function getCDF(collection, property) {
  var histDict = collection.reduceColumns({
    reducer: ee.Reducer.fixedHistogram(0, 101, 101),
    selectors: [property]
  });

  // THE FIX: Provide a safe mathematical fallback if the collection is perfectly empty
  var safeHist = ee.Algorithms.If(
    ee.Algorithms.IsEqual(histDict.get('histogram'), null),
    ee.Array(ee.List.repeat([0, 0], 101)), // Dummy array to prevent crashes
    ee.Array(histDict.get('histogram'))
  );
  
  var histArray = ee.Array(safeHist);
  
  // Extract just the counts (column index 1)
  var counts = histArray.slice(1, 1, 2).project([0]);
  
  // Calculate total features
  var total = ee.Number(counts.reduce(ee.Reducer.sum(), [0]).get([0]));
  
  // Protect against division by zero
  var safeTotal = ee.Algorithms.If(total.eq(0), 1, total);
  
  // Accumulate (cumulative sum) and divide by total for probabilities
  var cdf = counts.accum(0).divide(safeTotal);
  
  return cdf.toList();
}

// Calculate the CDFs
var may_cdf = getCDF(fc_may, 'BGR');
var sept_cdf = getCDF(fc_sept, 'BGR');
var rap_cdf = getCDF(sampled_rap, 'RAP_BGR');

// =========================================================================
// CHART FORMATTING & EXPORT
// =========================================================================

var xValues = ee.List.sequence(0, 100);

var cdfFC = ee.FeatureCollection(xValues.map(function(x) {
  var i = ee.Number(x);
  return ee.Feature(null, {
    'Threshold': i,
    'May 2019 (Drone)': may_cdf.get(i),
    'Sept 2019 (Drone)': sept_cdf.get(i),
    'RAP 2019 (10m)': rap_cdf.get(i)
  });
}));

var cdfChart = ui.Chart.feature.byFeature({
  features: cdfFC,
  xProperty: 'Threshold',
  yProperties: ['May 2019 (Drone)', 'Sept 2019 (Drone)', 'RAP 2019 (10m)']
})
.setChartType('LineChart')
.setOptions({
  title: 'Cumulative Probability Distribution: Bare Ground %',
  hAxis: {title: 'Bare Ground Percentage (%)', viewWindow: {min: 0, max: 100}},
  vAxis: {title: 'Cumulative Probability', viewWindow: {min: 0, max: 1}},
  lineWidth: 3,
  series: {
    0: {color: '#1f78b4'}, // Blue
    1: {color: '#ff7f00'}, // Orange
    2: {color: '#33a02c'}  // Green
  }
});

print(cdfChart);
