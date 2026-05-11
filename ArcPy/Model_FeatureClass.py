import arcpy
import os
import numpy as np
from scipy.ndimage import label, distance_transform_edt
import concurrent.futures
import multiprocessing

# ====================================================================
# USER INPUTS
# ====================================================================
input_gee_shapefile = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\Model_Grid_SHP\SRER_Model_Grid_UTM.shp'
source_tiff = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\Clipped_Classified_May_2019_UTM.tif'
output_fc = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_Grid_Metrics_Final'

# ====================================================================
# WORKER FUNCTION (Runs on multiple cores simultaneously)
# ====================================================================
def process_grid_cell(data_packet):
    """This function runs isolated on a separate CPU core"""
    oid, xmin, ymin, xmax, ymax, tiff_path, c_size = data_packet
    target_value = 3 
    
    try:
        # 1. Extract data for this specific cell
        lower_left = arcpy.Point(xmin, ymin)
        ncols = int(round((xmax - xmin) / c_size))
        nrows = int(round((ymax - ymin) / c_size))
        
        main_array = arcpy.RasterToNumPyArray(tiff_path, lower_left, ncols, nrows)
        total_pixels = main_array.size
        
        if total_pixels == 0:
            return oid, None

        # Create a boolean mask of bare ground
        is_bare = (main_array == target_value)

        # --- 1. BGR ---
        bgr_percent = (np.sum(is_bare) / total_pixels) * 100

        # --- 2. LPI (Using SciPy RegionGroup equivalent) ---
        # SciPy 'label' finds contiguous clusters natively in RAM
        labeled_array, num_features = label(is_bare)
        if num_features > 0:
            counts = np.bincount(labeled_array.ravel())
            # counts[0] is the background, counts[1:] are the patches
            max_patch_pixels = np.max(counts[1:]) if len(counts) > 1 else 0
            lpi_percent = (max_patch_pixels / total_pixels) * 100
        else:
            lpi_percent = 0.0

        # --- 3. Mean Fetch (Using SciPy EucDistance equivalent) ---
        # distance_transform_edt calculates distance to the nearest ZERO value.
        # So we invert the mask (~is_bare) to make obstacles zero.
        dist_array = distance_transform_edt(~is_bare) * c_size
        valid_fetch = dist_array[is_bare] # Only keep distances on bare ground
        mean_fetch_exact = np.mean(valid_fetch) if valid_fetch.size > 0 else 0.0

        # --- 4. CANOPY GAP FRACTIONS ---
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

        # Return the OID and the list of calculated metrics
        return oid, [bgr_percent, lpi_percent, mean_fetch_exact, f_0_24, f_25_50, f_51_100, f_101_200, f_gt_200]
        
    except Exception as e:
        return oid, None

# ====================================================================
# MAIN EXECUTION BLOCK
# ====================================================================
if __name__ == '__main__':
    arcpy.ResetEnvironments() 
    arcpy.env.overwriteOutput = True

    print('Projecting and prepping data...')
    spatial_ref_nad83 = arcpy.SpatialReference(26912)
    projected_grid = 'in_memory\\Projected_Grid'
    
    arcpy.management.Project(
        in_dataset=input_gee_shapefile, 
        out_dataset=projected_grid, 
        out_coor_system=spatial_ref_nad83, 
        transform_method='WGS_1984_(ITRF00)_To_NAD_1983'
    )

    arcpy.management.CopyFeatures(projected_grid, output_fc)

    new_fields = ['BGR_pct', 'LPI_pct', 'Fetch_m', 'Gap_0_24', 'Gap_25_50', 'Gap_51_100', 'Gap_101_200', 'Gap_gt_200']
    for field in new_fields:
        arcpy.management.AddField(output_fc, field, 'DOUBLE')

    cell_size = float(arcpy.management.GetRasterProperties(source_tiff, 'CELLSIZEX').getOutput(0))

    # --- 1. Gather all tasks ---
    print("Reading grid boundaries...")
    tasks = []
    with arcpy.da.SearchCursor(output_fc, ["OID@", "SHAPE@"]) as cursor:
        for oid, geom in cursor:
            ext = geom.extent
            # Create a "packet" of data for each grid cell to send to a CPU core
            tasks.append((oid, ext.XMin, ext.YMin, ext.XMax, ext.YMax, source_tiff, cell_size))

    total_cells = len(tasks)
    print(f"Found {total_cells} cells to process.")
    
    # --- 2. MULTIPROCESSING: Distribute tasks to all available CPU cores ---
    print(f"Spinning up {multiprocessing.cpu_count()} CPU cores. Hold on tight...")
    results_dict = {}
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # executor.map automatically handles distributing the 'tasks' list to the worker function
        for count, result in enumerate(executor.map(process_grid_cell, tasks), 1):
            oid, metrics = result
            if metrics:
                results_dict[oid] = metrics
            
            if count % 1000 == 0:
                print(f"Processed {count} / {total_cells} cells...")

    # --- 3. Write results back to the geodatabase ---
    print("Writing results back to the Feature Class attribute table...")
    update_fields = ["OID@"] + new_fields
    with arcpy.da.UpdateCursor(output_fc, update_fields) as cursor:
        for row in cursor:
            oid = row[0]
            if oid in results_dict:
                # Insert the 8 metrics into indices 1 through 8
                row[1:9] = results_dict[oid]
                cursor.updateRow(row)

    print('\nProcessing Complete! Your CPU can rest now.')
