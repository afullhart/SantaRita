// 1. Define the input assets
var class_may_2019 = ee.Image("users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_may_2019_5cm");
var class_sep_2019 = ee.Image("users/gponce/usda_ars/assets/images/aes/srer/suas/2019/full_ortho_classified_sep_2019_5cm");

// 2. Define the exact export parameters
var export_crs = 'EPSG:32612'; // WGS 1984 UTM Zone 12N
var export_scale = 0.05;       // 5cm pixel resolution

// 3. Export the May 2019 Image to Google Drive
Export.image.toDrive({
  image: class_may_2019,
  description: 'SRER_Classified_May_2019_UTM12N',
  folder: 'GEE_Downloads', // Change this to your preferred Drive folder
  scale: export_scale,
  crs: export_crs,
  region: class_may_2019.geometry(), // Automatically bounds to the image footprint
  maxPixels: 1e13 // Critical for 5cm data to prevent "too many pixels" errors
});

// 4. Export the September 2019 Image to Google Drive
Export.image.toDrive({
  image: class_sep_2019,
  description: 'SRER_Classified_Sep_2019_UTM12N',
  folder: 'GEE_Downloads',
  scale: export_scale,
  crs: export_crs,
  region: class_sep_2019.geometry(),
  maxPixels: 1e13
});
