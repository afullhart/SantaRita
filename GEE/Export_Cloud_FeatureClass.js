// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry().bounds();

// =========================================================================
// CLOUD MASKING LOGIC & HARMONIZATION
// =========================================================================
function maskClouds(qa) {
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  return dilatedCloud.and(cirrus).and(cloud).and(shadow);
}

// Harmonize the NIR band across sensors
function addNirL89(img) {
  return img.addBands(img.select('SR_B5').rename('NIR'));
}

function addNirL57(img) {
  return img.addBands(img.select('SR_B4').rename('NIR'));
}

// =========================================================================
// MONTHLY LANDSAT INVENTORY (1984 - 2025)
// =========================================================================
print('Scanning and filtering 40 years of Landsat imagery (Pixel-Level Median Strategy)...');

var start_year = 1984; 
var end_year = 2025;   
var years = ee.List.sequence(start_year, end_year);
var months = ee.List.sequence(1, 12);

// Map over every Year
var allMonthlyDataList = years.map(function(y) {
  y = ee.Number(y);
  
  // Map over every Month
  var monthlyData = months.map(function(m) {
    m = ee.Number(m);
    
    var startDate = ee.Date.fromYMD(y, m, 1);
    var endDate = startDate.advance(1, 'month'); 
    
    var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate).map(addNirL89);
    var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate).map(addNirL89);
    var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate).map(addNirL57);
    var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(bounds_geom).filterDate(startDate, endDate).map(addNirL57);
      
    var combined_month_col = l9.merge(l8).merge(l7).merge(l5);
    
    // --- CALCULATE CLEAR PIXEL FRACTION ---
    var scored_col = combined_month_col.map(function(img) {
      var qa = img.select('QA_PIXEL');
      var clear_mask = maskClouds(qa);
      
      var nir_scaled = img.select('NIR').multiply(0.0000275).add(-0.2);
      var dark_mask = nir_scaled.gt(0.15); 
      var native_mask = img.select('SR_B2').mask(); 
      
      // Combine the masks into a base validity mask
      var base_mask = clear_mask.and(dark_mask).and(native_mask);
      
      // --- THE FIX: Apply a 1-pixel spatial buffer to expand the masked areas ---
      var master_mask = base_mask.focal_min({radius: 1, kernelType: 'square', units: 'pixels'}).unmask(0).rename('Quality');
      
      var fraction = master_mask.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: bounds_geom,
        scale: 30,
        maxPixels: 1e9
      }).get('Quality');
      
      fraction = ee.Algorithms.If(ee.Algorithms.IsEqual(fraction, null), 0, fraction);
      return img.set('Local_Clear_Fraction', fraction);
    });
    
    var good_col = scored_col.filter(ee.Filter.gte('Local_Clear_Fraction', 0.20));
    var imgCount = good_col.size();
    
    var valid_ids = good_col.aggregate_array('system:id').join(',');
    
    var windowLabel = ee.Algorithms.If(
      imgCount.eq(0),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (No Data)'),
      ee.String(y.format('%d')).cat('-').cat(m.format('%02d')).cat(' (Median of ').cat(imgCount.format('%d')).cat(' imgs)')
    );

    return ee.Feature(bounds_geom, {
      'Year': y,
      'Month': m,
      'Start_Date': startDate.format('YYYY-MM-dd'),
      'Window_Label': windowLabel,
      'Image_Count': imgCount,
      'Valid_IDs': valid_ids, 
      'system:time_start': startDate.millis() 
    });
  });
  
  return monthlyData;
});

var monthlyIndexFc = ee.FeatureCollection(allMonthlyDataList.flatten())
  .filter(ee.Filter.neq('Image_Count', 0)) 
  .sort('system:time_start'); 

Export.table.toAsset({
  collection: monthlyIndexFc,
  description: 'Export_Landsat_Monthly_Index_Asset',
  assetId: 'projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat'
});
