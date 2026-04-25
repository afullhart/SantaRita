var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');

// Extract the raw Geometry from the first feature directly
var bounds_geom = bounds_fc.first().geometry();

var v_classified_may = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm');
var v_classified_sep = ee.Image('users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm');

var v_foot_prints = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_drone_footprints');

var v_srer_polys = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_ecological_states')
                     .map(function(ft){
                       return ft.set('area_ha', ft.area(1).divide(10000)); 
});

// Because bounds_geom is a Geometry, .bounds() safely returns an ee.Geometry bounding box
var v_extent = bounds_geom.bounds();

var sent2_ic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(v_extent)                     // FILTER LOCATION FIRST
  .filterDate('2019-05-01', '2019-05-31');    // FILTER DATE SECOND

var f_date = '2019-05-26', l_date = '2019-05-31';

var projSent2 = sent2_ic.filterDate(f_date, l_date)
                          .filterBounds(v_extent).first().select('B2').projection();

var sent2_im = sent2_ic.filterDate(f_date, l_date)
                          .filterBounds(v_extent)
                          .mosaic()
                          .clip(v_extent)
                          .setDefaultProjection({crs:projSent2.crs(), scale:projSent2.nominalScale()});
                           
var sent2_grid = v_extent.coveringGrid(projSent2, projSent2.nominalScale());

// This is one option to perform a spatial join
var v_spatial_filter = ee.Filter.intersects({leftField:'.geo', rightField:'.geo', maxError:1});

// Define a save all join.
var v_saveAllJoin = ee.Join.saveAll({matchesKey:'polys',});

// Perform spatial join of grids and polygons and keep only those grids that overlap at least 50% of footprint area.
var v_sensor_grids = sent2_grid;             // Sentinel2 
var v_desc = 'sent2_grids_joined';

// Apply the join.
var v_intersect1 = v_saveAllJoin.apply(v_sensor_grids, v_foot_prints, v_spatial_filter);

// Spatial join and filtering out grids with less than 100% area of overlap.
var v_fprint_filtered = v_sensor_grids.filterBounds(v_srer_polys);
var v_intersect = v_intersect1.map(function(feature){
  var v_poly = feature.geometry();
  var v_intersection = v_poly.intersection(v_foot_prints.geometry(), ee.ErrorMargin(1))
                            .intersection(v_srer_polys.geometry(), ee.ErrorMargin(1));
  var v_totalArea = v_poly.area(1);
  var v_overlapped = v_intersection.area(1).divide(v_totalArea);
  return feature.set({'overlapped': v_overlapped});
}).map(function (ft) {  
        return ft.set('polys',null);
});    

v_intersect = v_intersect.filter(ee.Filter.gte("overlapped", 1.0));

// Adding a unique sequenced id to featureCollection
var v_indexes = ee.List(v_intersect.aggregate_array('system:index'));
var v_ids = ee.List.sequence(1, v_intersect.size());
var v_idByIndex = ee.Dictionary.fromLists(v_indexes, v_ids);
var v_datasetWithId = v_intersect.map(function(feature){
  return feature.set('id', ee.Number(v_idByIndex.get(feature.get('system:index'))).toInt());
});

var v_sent2_joined_grids = v_datasetWithId;

var v_sent2_joined_grids = v_saveAllJoin.apply(v_sent2_joined_grids, v_srer_polys, v_spatial_filter)
                                        .map(function (ft){  // Set attributes to each grid out of the polygons surveyed
                                              var ft1 = ee.Feature(ee.List(ft.get('polys')).get(0));
                                              return ft.set({'Plant_Comm': ft1.get('Plant_Comm'),
                                                            'Pasture': ft1.get('Pasture'),
                                                            'Transect': ft1.get('Transect'),
                                                            'Utility' : ft1.get('Utility'),
                                                            'S_Desc' : ft1.get('S_Desc'),
                                                            'Exclosure' : ft1.get('Exclosure'),
                                                            'area_ha': ft1.get('area_ha'),
                                                            'polys':null,
                                                            'overlapped':null
                                              }); 
                                        });  



Export.table.toAsset({
  collection:v_sent2_joined_grids,
  assetId:'projects/ee-andrewfullhart/assets/SR_s2_grid_joined',
  description:'ftv_sentinel2_grid_srer_slud'
});


