import arcpy
import os
import numpy as np
from scipy.ndimage import label, distance_transform_edt
import concurrent.futures
import multiprocessing

# ====================================================================
# USER INPUTS
# ====================================================================
# Using the original raw shapefile to prevent tectonic shifting
input_gee_shapefile = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Model_Grid_Input_Shp\SRER_Grid.shp'

# Two 5cm class mosaics
may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'

output_fc = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_Grid_FeatureClass'

# The exact name of the field in your shapefile containing the month data
month_field_name = 'month' 

# ====================================================================
# WORKER FUNCTION (Runs on multiple cores simultaneously)
# ====================================================================
def process_grid_cell(data_packet):
  """This function runs isolated on a separate CPU core"""
  oid, xmin, ymin, xmax, ymax, tiff_path, c_size = data_packet
  
  # Class values
  target_value = 3  # Bare ground
  herb_value = 1    # Herb
  woody_value = 2   # Woody
  
  try:
    lower_left = arcpy.Point(xmin, ymin)
    ncols = int(round((xmax - xmin) / c_size))
    nrows = int(round((ymax - ymin) / c_size))
    
    main_array = arcpy.RasterToNumPyArray(tiff_path, lower_left, ncols, nrows)
    total_pixels = main_array.size
    
    if total_pixels == 0:
      return oid, None

    is_bare = (main_array == target_value)

    # --- 1. BGR ---
    bgr_percent = (np.sum(is_bare) / total_pixels) * 100

    # --- 2. LPI ---
    # FIX: Force 8-connectivity to match Earth Engine's "eightConnected: true"
    structure = np.ones((3, 3), dtype=int)
    labeled_array, num_features = label(is_bare, structure=structure)
    if num_features > 0:
      counts = np.bincount(labeled_array.ravel())
      max_patch_pixels = np.max(counts[1:]) if len(counts) > 1 else 0
      lpi_percent = (max_patch_pixels / total_pixels) * 100
    else:
      lpi_percent = 0.0

    # --- 3. Mean Fetch (Includes Obstacles as 0) ---
    # distance_transform_edt calculates distance from True (bare) to nearest False (obstacle)
    # False pixels (obstacles) automatically receive a distance of 0.0
    dist_array = distance_transform_edt(is_bare) * c_size
    
    # Calculate mean across ALL pixels in the cell matrix
    mean_fetch_exact = np.mean(dist_array)

    # --- 4. HERB-TO-WOODY RATIO ---
    herb_pixels = np.sum(main_array == herb_value)
    woody_pixels = np.sum(main_array == woody_value)
    herb_woody_ratio = (float(herb_pixels) / float(woody_pixels)) if woody_pixels > 0 else None

    # --- 5. CANOPY GAP FRACTIONS ---
    horizontal_transects = [main_array[i, :] for i in range(nrows)]
    vertical_transects = [main_array[:, j] for j in range(ncols)]
    all_transects = horizontal_transects + vertical_transects
    total_transect_length_m = len(all_transects) * (ncols * c_size)

    all_gap_lengths = []
    for t in all_transects:
      is_gap = (t == target_value)
      padded = np.concatenate(([False], is_gap, [False]))
      diffs = np.diff(padded.astype(int))
      starts = np.where(diffs == 1)[0]
      ends = np.where(diffs == -1)[0]
      all_gap_lengths.extend((ends - starts) * c_size)

    all_gap_lengths = np.array(all_gap_lengths)

    f_0_24   = (np.sum(all_gap_lengths[(all_gap_lengths < 0.25)]) / total_transect_length_m) * 100
    f_25_50  = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.25) & (all_gap_lengths <= 0.50)]) / total_transect_length_m) * 100
    f_51_100 = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.51) & (all_gap_lengths <= 1.00)]) / total_transect_length_m) * 100
    f_101_200= (np.sum(all_gap_lengths[(all_gap_lengths >= 1.01) & (all_gap_lengths <= 2.00)]) / total_transect_length_m) * 100
    f_gt_200 = (np.sum(all_gap_lengths[(all_gap_lengths > 2.00)]) / total_transect_length_m) * 100

    return oid, [bgr_percent, lpi_percent, mean_fetch_exact, herb_woody_ratio, f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200]
    
  except Exception as e:
    return oid, None

# ====================================================================
# MAIN EXECUTION BLOCK
# ====================================================================
if __name__ == '__main__':
  arcpy.ResetEnvironments() 
  arcpy.env.overwriteOutput = True

  print('Projecting and prepping data to WGS 84 UTM Zone 12N...')
  # Force WGS 84 UTM 12N to match GEE math exactly
  spatial_ref_wgs84 = arcpy.SpatialReference(32612) 
  projected_grid = 'in_memory\\Projected_Grid'
  
  arcpy.management.Project(
    in_dataset=input_gee_shapefile, 
    out_dataset=projected_grid, 
    out_coor_system=spatial_ref_wgs84 
  )
  arcpy.management.CopyFeatures(projected_grid, output_fc)

  new_fields = ['BGR_pct', 'LPI_pct', 'Fetch_m', 'Herb_Woody_Ratio', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
  for field in new_fields:
    arcpy.management.AddField(output_fc, field, 'DOUBLE')

  # Assuming both mosaics are 5cm, we just grab the cell size from one of them
  cell_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))

  # --- 1. Gather all tasks and assign correct seasonal TIFF ---
  print('Reading grid boundaries and seasonal routing...')
  tasks = []
  
  # We now read the OID, Geometry, AND the Month field
  with arcpy.da.SearchCursor(output_fc, ['OID@', 'SHAPE@', month_field_name]) as cursor:
    for oid, geom, month_val in cursor:
      
      # Flexible routing
      month_str = str(month_val).strip().lower()
      
      if month_str in ['5', '05', 'may']:
        active_tiff = may_tiff
      elif month_str in ['9', '09', 'sep', 'sept', 'september']:
        active_tiff = sep_tiff
      else:
        print(f"Warning: Unknown month '{month_val}' for OID {oid}. Skipping.")
        continue

      ext = geom.extent
      
      # THE SENTINEL-2 SNAP: Erase the GEE export drift by rounding to the nearest perfect 10m interval
      min_x = round(ext.XMin / 10.0) * 10.0
      min_y = round(ext.YMin / 10.0) * 10.0
      max_x = round(ext.XMax / 10.0) * 10.0
      max_y = round(ext.YMax / 10.0) * 10.0

      # Package the exact right TIFF path for this perfectly snapped cell
      tasks.append((oid, min_x, min_y, max_x, max_y, active_tiff, cell_size))

  total_cells = len(tasks)
  print(f'Found {total_cells} cells to process.')
  
  # --- 2. MULTIPROCESSING ---
  print(f'Spinning up {multiprocessing.cpu_count()} CPU cores. Hold on tight...')
  results_dict = {}
  
  with concurrent.futures.ProcessPoolExecutor() as executor:
    for count, result in enumerate(executor.map(process_grid_cell, tasks), 1):
      oid, metrics = result
      if metrics:
        results_dict[oid] = metrics
      
      if count % 1000 == 0:
        print(f'Processed {count} / {total_cells} cells...')

  # --- 3. Write results back ---
  print('Writing results back to the Feature Class attribute table...')
  update_fields = ['OID@'] + new_fields
  with arcpy.da.UpdateCursor(output_fc, update_fields) as cursor:
    for row in cursor:
      oid = row[0]
      if oid in results_dict:
        # Changed slice from row[1:9] to row[1:10] to account for 9 calculated fields
        row[1:10] = results_dict[oid]
        cursor.updateRow(row)

  print('\nProcessing Complete! Seasonal routing, geodetic snapping, and math fixes were successful.')
