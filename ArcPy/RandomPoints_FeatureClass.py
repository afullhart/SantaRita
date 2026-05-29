import arcpy
import os
import numpy as np
from scipy.ndimage import label, distance_transform_edt
import concurrent.futures
import multiprocessing

# ====================================================================
# USER INPUTS
# ====================================================================
# Original footprints to pull the "month" attribute from
input_footprints = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\srer_drone_veg_monitoring_footprints_2019\srer_drone_veg_monitoring_footprints_2019.shp'

# The generated feature classes from our previous steps
out_gdb = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb'
output_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m')
output_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m')
output_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m')

# TIFFs
may_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif'
sep_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_Sep_2019_UTM12N_Mosaic.tif'

month_field_name = 'month' 

# ====================================================================
# WORKER FUNCTION (Runs on multiple cores)
# ====================================================================
def process_grid_cell(data_packet):
    """This function runs isolated on a separate CPU core"""
    oid, xmin, ymin, xmax, ymax, tiff_path, c_size, is_circle = data_packet
    target_value = 3 
    
    try:
        lower_left = arcpy.Point(xmin, ymin)
        ncols = int(round((xmax - xmin) / c_size))
        nrows = int(round((ymax - ymin) / c_size))
        
        main_array = arcpy.RasterToNumPyArray(tiff_path, lower_left, ncols, nrows)
        
        # --- CIRCULAR MASKING LOGIC ---
        # If it's the 110m plot, mask out the corners of the bounding box
        if is_circle:
            Y, X = np.ogrid[:nrows, :ncols]
            # Center of the array
            center_r, center_c = nrows / 2.0, ncols / 2.0
            # Calculate distance of each pixel from center (in meters)
            dist_from_center = np.sqrt((X - center_c)**2 + (Y - center_r)**2) * c_size
            
            # Mask out everything outside the 55m radius
            mask = dist_from_center <= 55.0
            is_bare = (main_array == target_value) & mask
            total_pixels = np.sum(mask) # Only count pixels inside the circle
            
            # For Canopy Gap, zero out the corners so they aren't counted as gaps or canopy
            main_array[~mask] = 0 
        else:
            is_bare = (main_array == target_value)
            total_pixels = main_array.size
            
        if total_pixels == 0:
            return oid, None

        # --- 1. BGR ---
        bgr_percent = (np.sum(is_bare) / total_pixels) * 100

        # --- 2. LPI ---
        structure = np.ones((3, 3), dtype=int)
        labeled_array, num_features = label(is_bare, structure=structure)
        if num_features > 0:
            counts = np.bincount(labeled_array.ravel())
            max_patch_pixels = np.max(counts[1:]) if len(counts) > 1 else 0
            lpi_percent = (max_patch_pixels / total_pixels) * 100
        else:
            lpi_percent = 0.0

        # --- 3. Mean Fetch ---
        dist_array = distance_transform_edt(is_bare) * c_size
        
        if is_circle:
            # Only average the fetch values that are inside our circular mask
            mean_fetch_exact = np.mean(dist_array[mask])
        else:
            mean_fetch_exact = np.mean(dist_array)

        # --- 4. CANOPY GAP FRACTIONS ---
        horizontal_transects = [main_array[i, :] for i in range(nrows)]
        vertical_transects = [main_array[:, j] for j in range(ncols)]
        all_transects = horizontal_transects + vertical_transects
        
        # Total transect length is adjusted based on total valid pixels
        total_transect_length_m = len(all_transects) * (ncols * c_size)
        if is_circle:
            # Approximation for circle transects: 2 * (Area / c_size)
            total_transect_length_m = 2 * (total_pixels * c_size)

        all_gap_lengths = []
        for t in all_transects:
            # We must ignore the padding '0' values created by the circle mask
            is_gap = (t == target_value)
            padded = np.concatenate(([False], is_gap, [False]))
            diffs = np.diff(padded.astype(int))
            starts = np.where(diffs == 1)[0]
            ends = np.where(diffs == -1)[0]
            all_gap_lengths.extend((ends - starts) * c_size)

        all_gap_lengths = np.array(all_gap_lengths)

        if total_transect_length_m > 0:
            f_0_24   = (np.sum(all_gap_lengths[(all_gap_lengths < 0.25)]) / total_transect_length_m) * 100
            f_25_50  = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.25) & (all_gap_lengths <= 0.50)]) / total_transect_length_m) * 100
            f_51_100 = (np.sum(all_gap_lengths[(all_gap_lengths >= 0.51) & (all_gap_lengths <= 1.00)]) / total_transect_length_m) * 100
            f_101_200= (np.sum(all_gap_lengths[(all_gap_lengths >= 1.01) & (all_gap_lengths <= 2.00)]) / total_transect_length_m) * 100
            f_gt_200 = (np.sum(all_gap_lengths[(all_gap_lengths > 2.00)]) / total_transect_length_m) * 100
        else:
            f_0_24 = f_25_50 = f_51_100 = f_101_200 = f_gt_200 = 0.0

        return oid, [bgr_percent, lpi_percent, mean_fetch_exact, f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200]
        
    except Exception as e:
        return oid, None

# ====================================================================
# EXECUTION LOGIC
# ====================================================================
def calculate_metrics_for_fc(target_fc, raw_footprints, cell_size, is_circle=False):
    print(f'\n--- Processing Metrics for {os.path.basename(target_fc)} ---')
    
    # 1. Join the 'month' attribute from the raw footprints
    print('  -> Retrieving "month" attribute via Spatial Join...')
    temp_join = r'memory\temp_join'
    arcpy.analysis.SpatialJoin(
        target_features=target_fc,
        join_features=raw_footprints,
        out_feature_class=temp_join,
        join_operation="JOIN_ONE_TO_ONE",
        match_option="INTERSECT"
    )
    
    # Delete original features and replace with the joined ones
    arcpy.management.Delete(target_fc)
    arcpy.management.CopyFeatures(temp_join, target_fc)
    arcpy.management.Delete(temp_join)

    # 2. Add New Fields
    new_fields = ['BGR_pct', 'LPI_pct', 'Fetch_m', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
    for field in new_fields:
        arcpy.management.AddField(target_fc, field, 'DOUBLE')

    # 3. Gather Tasks
    print('  -> Reading boundaries and seasonal routing...')
    tasks = []
    
    with arcpy.da.SearchCursor(target_fc, ['OID@', 'SHAPE@', month_field_name]) as cursor:
        for oid, geom, month_val in cursor:
            month_str = str(month_val).strip().lower()
            if month_str in ['5', '05', 'may']:
                active_tiff = may_tiff
            elif month_str in ['9', '09', 'sep', 'sept', 'september']:
                active_tiff = sep_tiff
            else:
                print(f"Warning: Unknown month '{month_val}' for OID {oid}. Defaulting to May.")
                active_tiff = may_tiff

            ext = geom.extent
            min_x = round(ext.XMin / 10.0) * 10.0
            min_y = round(ext.YMin / 10.0) * 10.0
            max_x = round(ext.XMax / 10.0) * 10.0
            max_y = round(ext.YMax / 10.0) * 10.0

            tasks.append((oid, min_x, min_y, max_x, max_y, active_tiff, cell_size, is_circle))

    total_cells = len(tasks)
    print(f'  -> Found {total_cells} features. Spinning up {multiprocessing.cpu_count()} CPU cores...')
    
    # 4. Multiprocessing
    results_dict = {}
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for count, result in enumerate(executor.map(process_grid_cell, tasks), 1):
            oid, metrics = result
            if metrics:
                results_dict[oid] = metrics
            if count % 1000 == 0:
                print(f'  -> Processed {count} / {total_cells} cells...')

    # 5. Write back to Geodatabase
    print('  -> Writing results back to the attribute table...')
    update_fields = ['OID@'] + new_fields
    with arcpy.da.UpdateCursor(target_fc, update_fields) as cursor:
        for row in cursor:
            oid = row[0]
            if oid in results_dict:
                row[1:9] = results_dict[oid]
                cursor.updateRow(row)

    print(f'  -> Successfully updated {os.path.basename(target_fc)}!')

if __name__ == '__main__':
    arcpy.ResetEnvironments() 
    arcpy.env.overwriteOutput = True

    # Retrieve raster cell size from the May TIFF
    print('Reading base raster properties...')
    pixel_size = float(arcpy.management.GetRasterProperties(may_tiff, 'CELLSIZEX').getOutput(0))
    print(f'Detected pixel size: {pixel_size}m')

    # Execute for all three datasets
    # target_fc, raw_footprints_for_spatial_join, pixel_size, is_circle
    calculate_metrics_for_fc(output_grid_10m, input_footprints, pixel_size, is_circle=False)
    calculate_metrics_for_fc(output_grid_30m, input_footprints, pixel_size, is_circle=False)
    calculate_metrics_for_fc(output_nri_plots, input_footprints, pixel_size, is_circle=True)

    print('\nAll Extracted Metrics Processing Complete!')
