//Can't use projections for shapefile exports
//(No UTM)

var model_grid_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_s2_model_grid_utm');

// Export the FeatureCollection to Google Drive as a Shapefile
Export.table.toDrive({
  collection: model_grid_fc,
  description: 'SRER_Model_Grid_Export',
  folder: 'GEE_Downloads', 
  fileFormat: 'SHP'
});
