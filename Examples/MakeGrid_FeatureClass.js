var bounds_fc = ee.FeatureCollection('projects/ee-andrewfullhart/assets/SR_bounds');
var bounds_ft = ee.Feature(bounds_fc.first());

///////////////////////////////////////////////////////
//Create Landsat Grid Outline
///////////////////////////////////

var roi = bounds_ft.geometry();
Map.addLayer(bounds_ft);

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
           .filterBounds(roi)
           .first();

var landsatProj = l8.select(0).projection();

var vectorGrid = roi.coveringGrid(landsatProj, 30);

var gridStyle = {
  color: 'red',
  fillColor: '00000000', 
  width: 1               
};

Map.setCenter(-110.85, 31.80, 14);
var visParams = {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 7000, max: 12000};
Map.addLayer(l8, visParams, 'Landsat Image');
Map.addLayer(vectorGrid.style(gridStyle), {}, 'Vector 30m Grid');


