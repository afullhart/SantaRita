
var roi = ee.Geometry.BBox(-111.00, 31.70, -110.75, 31.95);

// var naip_ic = ee.ImageCollection('USDA/NAIP/DOQQ')
//               .filter(ee.Filter.date('2017-05-01', '2017-05-31'))
//               .filterBounds(roi);

var naip_ic = ee.ImageCollection('USDA/NAIP/DOQQ')
               .filter(ee.Filter.date('2019-06-01', '2019-06-30'))
               .filterBounds(roi);

print(naip_ic.size());

var naip_im = naip_ic.mosaic();

Map.setCenter(-110.90, 31.75, 15);
print(naip_im);

var ndvi = ee.Image(naip_im.select('N').subtract(naip_im.select('R'))).divide(naip_im.select('N').add(naip_im.select('R')));

// Define the visualization parameters
var ndviVis = {
  min: -0.1, 
  max: 0.8,
  palette: [
    '#d73027', // Red (Water, rock, bare soil)
    '#f46d43', // Orange
    '#fdae61', // Light Orange
    '#fee08b', // Yellow (Very sparse vegetation)
    '#d9ef8b', // Light Green
    '#a6d96a', // Green
    '#66bd63', // Darker Green
    '#1a9850'  // Dark Green (Dense, healthy vegetation)
  ]
};

Map.addLayer(ndvi, ndviVis, 'NAIP NDVI');
