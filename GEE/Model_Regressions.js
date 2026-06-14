// =========================================================================
// SETUP & ASSETS
// =========================================================================
// Point directly to the newly generated Landsat Training FeatureClass asset
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_landsat_model_grid_utm');

// =========================================================================
// PREDICTOR VARIABLES (LANDSAT ALIGNED)
// =========================================================================
// Updated to match the exact bands and indices extracted in the Landsat script.
// (MCARI is omitted because Landsat lacks the necessary Red Edge band).
var inputProps = [
  'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 
  'NDVI', 'BSI', 'NBR2', 'slope', 'illumination', 'aspect'
];

// =========================================================================
// DROP NULLS & LOG TRANSFORM SKEWED RATIO DATA
// =========================================================================
// We kept null features in the asset to align with the shapefile geometry,
// but the classifier requires valid numbers. We drop them here right before training.
var validPropsCheck = inputProps.concat(['Herb_Woody_Ratio']);
fc = fc.filter(ee.Filter.notNull(validPropsCheck));

// --- Log Transform HWR Values to avoid exploding values ---
fc = fc.map(function(ft) {
  var raw_hwr = ee.Number(ft.get('Herb_Woody_Ratio'));
  
  // Natural Log (ln) transformation. We add 1 to handle cases where raw_hwr is 0.
  var log_hwr = raw_hwr.add(1).log(); 
  
  return ft.set('Log_HWR', log_hwr);
});

// =========================================================================
// K-FOLDS CROSS VALIDATION (K=5) - UNIFIED MODEL
// =========================================================================
var k_folds = 5;
var v_list_seeds = ee.List([123, 456, 789, 111, 333]);
var fold_list = ee.List.sequence(0, 4);

var hyperpars = {
  numberOfTrees: 400,
  shrinkage: 0.05,  
  samplingRate: 0.7,
  maxNodes: 32, 
  loss: 'Huber',
  seed: null
};

var fc_folds = fc.randomColumn('random', 123).map(function(ft) {
  return ft.set('fold', ee.Number(ft.get('random')).multiply(k_folds).floor());
});

function runCV(targetProp) {
  var fold_metrics = fold_list.map(function(fold) {
    var i = ee.Number(fold);
    var current_seed = v_list_seeds.get(i);
    
    var train_fc = fc_folds.filter(ee.Filter.neq('fold', i));
    var test_fc = fc_folds.filter(ee.Filter.eq('fold', i));
    
    var gtb = ee.Classifier.smileGradientTreeBoost({
      numberOfTrees: hyperpars.numberOfTrees, 
      shrinkage: hyperpars.shrinkage, 
      samplingRate: hyperpars.samplingRate, 
      maxNodes: hyperpars.maxNodes, 
      loss: hyperpars.loss, 
      seed: current_seed
    }).setOutputMode('REGRESSION').train({features: train_fc, classProperty: targetProp, inputProperties: inputProps});
    
    var tested = test_fc.classify({classifier: gtb, outputName: 'predicted'});
    
    var with_sq_err = tested.map(function(ft) {
      var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
      return ft.set('sq_diff', diff.multiply(diff));
    });

    var may_mse = with_sq_err.filter(ee.Filter.eq('Month', 'May'))
      .randomColumn('limit_sort').sort('limit_sort').limit(1000) 
      .reduceColumns({reducer: ee.Reducer.mean(), selectors: ['sq_diff']}).get('mean');
      
    var sep_mse = with_sq_err.filter(ee.Filter.eq('Month', 'Sept'))
      .randomColumn('limit_sort').sort('limit_sort').limit(1000) 
      .reduceColumns({reducer: ee.Reducer.mean(), selectors: ['sq_diff']}).get('mean');
    
    return ee.Dictionary({
      'May_RMSE': ee.Number(may_mse).sqrt(),
      'Sept_RMSE': ee.Number(sep_mse).sqrt()
    });
  });
  
  var may_rmse_list = fold_metrics.map(function(d) { return ee.Dictionary(d).get('May_RMSE'); });
  var sep_rmse_list = fold_metrics.map(function(d) { return ee.Dictionary(d).get('Sept_RMSE'); });
  
  return ee.Dictionary({
    'May_Median': ee.List(may_rmse_list).reduce(ee.Reducer.median()),
    'Sept_Median': ee.List(sep_rmse_list).reduce(ee.Reducer.median())
  });
}

print('--- K-FOLDS RMSE (STRATIFIED TESTING ERROR) ---');
var cv_bgr = runCV('BGR');
print('BGR RMSE -> May:', cv_bgr.get('May_Median'), '| Sept:', cv_bgr.get('Sept_Median'));
var cv_lpi = runCV('LPI');
print('LPI RMSE -> May:', cv_lpi.get('May_Median'), '| Sept:', cv_lpi.get('Sept_Median'));
var cv_mft = runCV('MFT');
print('MFT RMSE -> May:', cv_mft.get('May_Median'), '| Sept:', cv_mft.get('Sept_Median'));
var cv_hwr = runCV('Log_HWR');
print('Log_HWR RMSE -> May:', cv_hwr.get('May_Median'), '| Sept:', cv_hwr.get('Sept_Median'));
print('-------------------------------------------------------');

// =========================================================================
// FINAL MODEL TRAINING
// =========================================================================

hyperpars.seed = 123; 

var model_bgr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'MFT', inputProperties: inputProps});
var model_hwr = ee.Classifier.smileGradientTreeBoost(hyperpars)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'Log_HWR', inputProperties: inputProps});

// -------------------------------------------------------------------------
// CALCULATE & PRINT FINAL MODEL (STRATIFIED TRAINING ERROR)
// -------------------------------------------------------------------------
function getStratifiedTrainingRMSE(trainedModel, dataset, targetProp, monthFilter) {
  var filtered_dataset = dataset.filter(ee.Filter.eq('Month', monthFilter));
  var classified = filtered_dataset.classify({
    classifier: trainedModel, 
    outputName: 'predicted'
  });
  var mse = classified.map(function(ft) {
    var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
    return ft.set('sq_diff', diff.multiply(diff));
  })
  .randomColumn('limit_sort').sort('limit_sort').limit(1000) 
  .reduceColumns({
    reducer: ee.Reducer.mean(),
    selectors: ['sq_diff']
  }).get('mean');
  
  return ee.Number(mse).sqrt();
}

print('--- FITTED MODEL RMSE (STRATIFIED TRAINING ERROR) ---');
print('Fitted BGR RMSE -> May:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'Sept'));
print('Fitted LPI RMSE -> May:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'Sept'));
print('Fitted MFT RMSE -> May:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'Sept'));
print('Fitted Log_HWR RMSE -> May:', getStratifiedTrainingRMSE(model_hwr, fc, 'Log_HWR', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_hwr, fc, 'Log_HWR', 'Sept'));
print('-----------------------------------------------------');

var fc_may = fc.filter(ee.Filter.eq('Month', 'May'));
var fc_sep = fc.filter(ee.Filter.eq('Month', 'Sept'));

function predictFeatures(dataset) {
  var p_bgr = dataset.classify({classifier: model_bgr, outputName: 'Pred_BGR'});
  var p_lpi = p_bgr.classify({classifier: model_lpi, outputName: 'Pred_LPI'});
  var p_mft = p_lpi.classify({classifier: model_mft, outputName: 'Pred_MFT'});
  return p_mft.classify({classifier: model_hwr, outputName: 'Pred_Log_HWR'});
}

var final_predictions = predictFeatures(fc_may).merge(predictFeatures(fc_sep));

var export_csv = final_predictions.map(function(ft) {
  return ee.Feature(null, { 
    'Month': ft.get('Month'),
    'True_BGR': ft.get('BGR'), 'Predicted_BGR': ft.get('Pred_BGR'),
    'True_LPI': ft.get('LPI'), 'Predicted_LPI': ft.get('Pred_LPI'),
    'True_MFT': ft.get('MFT'), 'Predicted_MFT': ft.get('Pred_MFT'),
    'True_Log_HWR': ft.get('Log_HWR'), 'Predicted_Log_HWR': ft.get('Pred_Log_HWR')
  });
});

// =========================================================================
// SCATTER PLOTS & R-SQUARED VISUALIZATION
// =========================================================================
print('Generating Predicted vs. True scatter plots (n=4999 subset)...');

function makeScatterChart(dataset, trueProp, predProp, title, colorCode) {
  var sampledDataset = dataset.randomColumn('chart_sort').sort('chart_sort').limit(4999);

  var chart = ui.Chart.feature.byFeature({
    features: sampledDataset,
    xProperty: trueProp,
    yProperties: [predProp]
  })
  .setChartType('ScatterChart')
  .setOptions({
    title: title + ': Predicted vs. True (Sampled)',
    hAxis: {title: 'True ' + title},
    vAxis: {title: 'Predicted ' + title},
    pointSize: 2,
    colors: [colorCode],
    trendlines: { 
      0: { 
        type: 'linear', 
        color: 'black', 
        lineWidth: 2, 
        opacity: 0.8, 
        showR2: true, 
        visibleInLegend: true 
      } 
    }
  });
  print(chart);
}

makeScatterChart(export_csv, 'True_BGR', 'Predicted_BGR', 'Bare Ground (BGR %)', '#d73027'); 
makeScatterChart(export_csv, 'True_LPI', 'Predicted_LPI', 'Large Patch Index (LPI %)', '#fc8d59'); 
makeScatterChart(export_csv, 'True_MFT', 'Predicted_MFT', 'Mean Fetch (MFT m)', '#1a9850'); 
makeScatterChart(export_csv, 'True_Log_HWR', 'Predicted_Log_HWR', 'Log Herb-to-Woody Ratio', '#4575b4'); 

// =========================================================================
// EXPORT TO DRIVE
// =========================================================================
Export.table.toDrive({
  collection: export_csv,
  description: 'SRER_Metrics_True_vs_Predicted_Landsat',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});
