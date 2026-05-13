// =========================================================================
// USER INPUTS & DATES
// =========================================================================

// Dry Season (Pre-Monsoon) Dates
var may_start = '2019-05-26';
var may_end   = '2019-05-31';

// Post-Monsoon Dates
var sep_start = '2019-09-10';
var sep_end   = '2019-09-20';

// =========================================================================
// SETUP & ASSETS
// =========================================================================
var fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

// =========================================================================
// SENTINEL-2 PREDICTOR EXTRACTION (OPTION A: MEDIAN COMPOSITE)
// =========================================================================
var projSent2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)
  .filterDate(may_start, may_end)
  .first().select('B2').projection();

function buildS2Composite(startDate, endDate) {
  var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(v_extent)                     
    .filterDate(startDate, endDate)
    // OPTION A: Mask clouds > 20% probability
    .map(function(img) {
      var cloudMask = img.select('MSK_CLDPRB').lt(20);
      return img.updateMask(cloudMask);
    });                

  var sent2_im = sent2_ic
    // OPTION A: Median reducer to blend clear pixels
    .median() 
    .clip(v_extent)
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

var sent2_may = buildS2Composite(may_start, may_end);
var sent2_sep = buildS2Composite(sep_start, sep_end);

var inputProps = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12', 'NDVI', 'MCARI', 'BSI', 'NBR2'];

// =========================================================================
// K-FOLDS CROSS VALIDATION (K=5) - UNIFIED MODEL
// =========================================================================
var k_folds = 5;
var v_list_seeds = ee.List([123, 456, 789, 111, 333]);
var fold_list = ee.List.sequence(0, 4);

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
      numberOfTrees: 300, shrinkage: 0.01, samplingRate: 0.7, maxNodes: 12, seed: current_seed
    }).setOutputMode('REGRESSION').train({features: train_fc, classProperty: targetProp, inputProperties: inputProps});
    
    var tested = test_fc.classify({classifier: gtb, outputName: 'predicted'});
    
    var with_sq_err = tested.map(function(ft) {
      var diff = ee.Number(ft.get('predicted')).subtract(ee.Number(ft.get(targetProp)));
      return ft.set('sq_diff', diff.multiply(diff));
    });

    var may_mse = with_sq_err.filter(ee.Filter.eq('Month', 'May'))
      .reduceColumns({reducer: ee.Reducer.mean(), selectors: ['sq_diff']}).get('mean');
      
    var sep_mse = with_sq_err.filter(ee.Filter.eq('Month', 'Sept'))
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

print('--- UNIFIED K-FOLDS RMSE (STRATIFIED TESTING ERROR) ---');
var cv_bgr = runCV('BGR');
print('BGR RMSE -> May:', cv_bgr.get('May_Median'), '| Sept:', cv_bgr.get('Sept_Median'));
var cv_lpi = runCV('LPI');
print('LPI RMSE -> May:', cv_lpi.get('May_Median'), '| Sept:', cv_lpi.get('Sept_Median'));
var cv_mft = runCV('MFT');
print('MFT RMSE -> May:', cv_mft.get('May_Median'), '| Sept:', cv_mft.get('Sept_Median'));
print('-------------------------------------------------------');

// =========================================================================
// FINAL MODEL TRAINING
// =========================================================================
var regularized_params = {numberOfTrees: 300, shrinkage: 0.01, samplingRate: 0.7, maxNodes: 12, seed: 123};

var model_bgr = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'BGR', inputProperties: inputProps});
var model_lpi = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'LPI', inputProperties: inputProps});
var model_mft = ee.Classifier.smileGradientTreeBoost(regularized_params)
  .setOutputMode('REGRESSION').train({features: fc, classProperty: 'MFT', inputProperties: inputProps});


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
  }).reduceColumns({
    reducer: ee.Reducer.mean(),
    selectors: ['sq_diff']
  }).get('mean');
  
  return ee.Number(mse).sqrt();
}

print('--- UNIFIED FITTED MODEL RMSE (STRATIFIED TRAINING ERROR) ---');
print('Fitted BGR RMSE -> May:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_bgr, fc, 'BGR', 'Sept'));
print('Fitted LPI RMSE -> May:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_lpi, fc, 'LPI', 'Sept'));
print('Fitted MFT RMSE -> May:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'May'), '| Sept:', getStratifiedTrainingRMSE(model_mft, fc, 'MFT', 'Sept'));
print('-------------------------------------------------------');


// =========================================================================
// SEASON-AWARE EXPORT (MATCHES CONSOLE EXACTLY)
// =========================================================================
var fc_may = fc.filter(ee.Filter.eq('Month', 'May'));
var fc_sep = fc.filter(ee.Filter.eq('Month', 'Sept'));

function predictFeatures(dataset) {
  var p_bgr = dataset.classify({classifier: model_bgr, outputName: 'Pred_BGR'});
  var p_lpi = p_bgr.classify({classifier: model_lpi, outputName: 'Pred_LPI'});
  return p_lpi.classify({classifier: model_mft, outputName: 'Pred_MFT'});
}

var final_predictions = predictFeatures(fc_may).merge(predictFeatures(fc_sep));

var export_csv = final_predictions.map(function(ft) {
  return ee.Feature(null, { 
    'Month': ft.get('Month'),
    'True_BGR': ft.get('BGR'), 'Predicted_BGR': ft.get('Pred_BGR'),
    'True_LPI': ft.get('LPI'), 'Predicted_LPI': ft.get('Pred_LPI'),
    'True_MFT': ft.get('MFT'), 'Predicted_MFT': ft.get('Pred_MFT')
  });
});

Export.table.toDrive({
  collection: export_csv,
  description: 'SRER_Metrics_True_vs_Predicted_Unified_Fixed',
  folder: 'GEE_Downloads',
  fileFormat: 'CSV'
});
