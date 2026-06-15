// =========================================================================
// SETUP & ASSETS
// =========================================================================
var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_geom = bounds_fc.first().geometry();
var v_extent = bounds_geom.bounds();

var monthly_index = ee.FeatureCollection('projects/ee-andrewfullhart/assets/Cloud_FeatureClass_Landsat');

var rgbVis = {
  bands: ['red', 'green', 'blue'],
  min: 0.0,
  max: 0.3, 
  gamma: 1.2
};

// =========================================================================
// LANDSAT HARMONIZATION & CLOUD MASKING FUNCTIONS
// =========================================================================
function maskClouds(qa) {
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var shadow = qa.bitwiseAnd(1 << 4).eq(0);
  return dilatedCloud.and(cirrus).and(cloud).and(shadow);
}

function prepOLI(img) {
  var scaled = img.select(['SR_B5', 'SR_B4', 'SR_B3', 'SR_B2'])
                  .multiply(0.0000275).add(-0.2); 
                  
  var qa = img.select('QA_PIXEL');
  var qaMask = maskClouds(qa);
  
  var darkMask = scaled.select('SR_B5').gt(0.15); 
  var baseMask = qaMask.and(darkMask);
  
  var finalMask = baseMask.focal_min({radius: 30, kernelType: 'square', units: 'meters'});

  var optical = scaled.select(['SR_B4', 'SR_B3', 'SR_B2'])
                      .rename(['red', 'green', 'blue']);
                      
  return optical.updateMask(finalMask).copyProperties(img, ['system:time_start']);
}

function prepTM(img) {
  var scaled = img.select(['SR_B4', 'SR_B3', 'SR_B2', 'SR_B1'])
                  .multiply(0.0000275).add(-0.2); 
                  
  var qa = img.select('QA_PIXEL');
  var qaMask = maskClouds(qa);
  
  var darkMask = scaled.select('SR_B4').gt(0.15); 
  var baseMask = qaMask.and(darkMask);

  var finalMask = baseMask.focal_min({radius: 30, kernelType: 'square', units: 'meters'});

  var optical = scaled.select(['SR_B3', 'SR_B2', 'SR_B1'])
                      .rename(['red', 'green', 'blue']);
                      
  return optical.updateMask(finalMask).copyProperties(img, ['system:time_start']);
}

// =========================================================================
// UI APP & INTERACTIVITY
// =========================================================================
var panel = ui.Panel({style: {width: '350px', padding: '15px'}});
ui.root.insert(0, panel);

var title = ui.Label('Landsat Monthly Median Viewer', {fontWeight: 'bold', fontSize: '20px'});
var desc = ui.Label('Generates a pristine, cloud-masked monthly median composite using only high-quality images dynamically filtered from the historical archive.', {fontSize: '12px', color: '#555'});

var yearLabel = ui.Label('Select Year:', {fontWeight: 'bold'});
var yearSlider = ui.Slider({min: 1984, max: 2025, value: 2019, step: 1, style: {stretch: 'horizontal'}});

var monthLabel = ui.Label('Select Month:', {fontWeight: 'bold'});
var monthSlider = ui.Slider({min: 1, max: 12, value: 9, step: 1, style: {stretch: 'horizontal'}});

var statusBox = ui.Label({value: 'Ready. Move sliders to load imagery...', style: {color: 'blue', margin: '20px 0', whiteSpace: 'pre-wrap'}});

panel.add(title).add(desc).add(yearLabel).add(yearSlider).add(monthLabel).add(monthSlider).add(statusBox);

Map.centerObject(bounds_geom, 13);
Map.setOptions('HYBRID'); 

// =========================================================================
// RENDER FUNCTION
// =========================================================================
function updateMap() {
  var y = yearSlider.getValue();
  var m = monthSlider.getValue();
  
  statusBox.setValue('Checking asset index for pristine imagery...');
  statusBox.style().set('color', 'orange');
  
  var index_query = monthly_index
    .filter(ee.Filter.eq('Year', y))
    .filter(ee.Filter.eq('Month', m));
    
  index_query.size().evaluate(function(count) {
    if (count === 0) {
      statusBox.setValue('Asset indicates No Clean Data available for ' + y + '-' + m);
      statusBox.style().set('color', 'red');
      Map.layers().reset(); 
      Map.addLayer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds');
      return;
    }
    
    var metadata = ee.Feature(index_query.first());
    
    metadata.evaluate(function(feature) {
      var props = feature.properties;
      var label = props.Window_Label;
      var valid_ids_str = props.Valid_IDs; 
      
      statusBox.setValue('SUCCESS!\n' + label + '\nFetching exact images...');
      statusBox.style().set('color', 'green');
      
      var fullIdList = ee.String(valid_ids_str).split(',');
      
      var indexList = fullIdList.map(function(id) {
        var parts = ee.String(id).split('/');
        return parts.get(parts.length().subtract(1));
      });
      
      var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI);
      var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepOLI);
      var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM);
      var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filter(ee.Filter.inList('system:index', indexList)).map(prepTM);
        
      var combined_month_col = l9.merge(l8).merge(l7).merge(l5);
      var monthly_median = combined_month_col.median().clip(v_extent);
      
      // =====================================================================
      // COMPLETELY REPLACE THE IMAGE WITH A SYNTHETIC MEDIAN
      // =====================================================================
      if (y === 2019 && m === 9) {
        statusBox.setValue(statusBox.getValue() + '\nReplacing L7 stripes with Aug 31/Oct 2 synthetic median...');
        
        var aug31_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
          .filterBounds(bounds_geom)
          .filterDate('2019-08-31', '2019-09-01')
          .map(prepOLI); 
          
        var oct02_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
          .filterBounds(bounds_geom)
          .filterDate('2019-10-02', '2019-10-03')
          .map(prepOLI);

        monthly_median = aug31_col.merge(oct02_col).median().clip(v_extent);
      }
      // =====================================================================

      Map.layers().reset(); 
      Map.layers().set(0, ui.Map.Layer(monthly_median, rgbVis, 'Median ' + y + '-' + m));
      Map.layers().set(1, ui.Map.Layer(bounds_fc.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'SRER Bounds'));
    });
  });
}

yearSlider.onChange(updateMap);
monthSlider.onChange(updateMap);

updateMap();
