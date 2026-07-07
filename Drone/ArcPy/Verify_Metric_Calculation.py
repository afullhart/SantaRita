import arcpy
import os
import numpy as np
from arcpy.sa import *

# ====================================================================
# USER INPUT
# ====================================================================
source_tiff = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif"
raw_shapefile = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Model_Grid_Input_Shp\SRER_Grid.shp" 
grid_fc = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_Grid_Metrics_Final"

target_oid = 11622

# ====================================================================
# SETUP & DATA PREP
# ====================================================================
arcpy.ResetEnvironments() 
arcpy.env.workspace = "memory"
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

print("Importing raw shapefile to Geodatabase and enforcing projection...")

try:
    # 1. IMPORT & PROJECT
    sr_wgs84_utm12 = arcpy.SpatialReference(32612)
    arcpy.management.Project(
        in_dataset=raw_shapefile,
        out_dataset=grid_fc,
        out_coor_system=sr_wgs84_utm12
    )
    print("--> Shapefile successfully imported to Geodatabase!")

    cell_size = float(arcpy.management.GetRasterProperties(source_tiff, "CELLSIZEX").getOutput(0))
    
    print(f"\nReading exact boundaries for OID {target_oid}...")
    
    # 2. DYNAMIC EXTRACTION WITH SENTINEL-2 SNAP
    with arcpy.da.SearchCursor(grid_fc, ["SHAPE@"], f"OBJECTID = {target_oid}") as cursor:
        for row in cursor:
            ext = row[0].extent
            
            # THE SNAP: mathematically erase the GEE decimal drift
            min_x = round(ext.XMin / 10.0) * 10.0
            min_y = round(ext.YMin / 10.0) * 10.0
            max_x = round(ext.XMax / 10.0) * 10.0
            max_y = round(ext.YMax / 10.0) * 10.0

    print(f"Anchoring test raster to SNAPPED coordinates:\nX: {min_x} to {max_x}\nY: {min_y} to {max_y}")

    # Enforce exactly 200x200 extraction based on 10m grid and 5cm pixels
    ncols = 200
    nrows = 200
    lower_left = arcpy.Point(min_x, min_y)

    print(f"Slicing exact pixel matrix ({ncols} cols x {nrows} rows)...")

    # 3. Extract directly to NumPy (Injecting NoData handling)
    main_array = arcpy.RasterToNumPyArray(source_tiff, lower_left, ncols, nrows, nodata_to_value=-9999)
    
    valid_mask = main_array != -9999
    total_valid_pixels = np.sum(valid_mask)
    print(f"Total Valid Pixels: {total_valid_pixels}")
    
    if total_valid_pixels != 40000:
        print("!! WARNING: Did not extract exactly 40,000 valid pixels. Raster may not fully cover this grid cell. !!")
    else:
        print("SUCCESS: Exactly 40,000 valid pixels extracted (10m x 10m box).")
    
    # 4. Convert BACK to a temporary raster
    temp_raster = arcpy.NumPyArrayToRaster(main_array, lower_left, cell_size, cell_size)
    in_raster = Int(temp_raster)

    # ====================================================================
    # CORE METRICS (BGR, LPI, FETCH, HERB/WOODY)
    # ====================================================================
    print("\nCalculating Core Metrics...")
    target_value = 3 
    
    # BGR
    bare_ground_pixels = np.sum(main_array[valid_mask] == target_value)
    bgr_percent = (float(bare_ground_pixels) / float(total_valid_pixels)) * 100

    # LPI
    isolated_class = Con(in_raster == target_value, 1)
    grouped_patches = RegionGroup(isolated_class, "EIGHT", "WITHIN", "NO_LINK")
    patch_array = arcpy.RasterToNumPyArray(grouped_patches)
    
    patch_ids, counts = np.unique(patch_array, return_counts=True)
    valid_patch_mask = patch_ids > 0 
    max_patch_pixels = np.max(counts[valid_patch_mask]) if np.any(valid_patch_mask) else 0
    
    lpi_percent = (float(max_patch_pixels) / float(total_valid_pixels)) * 100

    # HERB-TO-WOODY RATIO & COVER PERCENTAGES
    herb_val = 1
    woody_val = 2
    herb_pixels = np.sum(main_array[valid_mask] == herb_val)
    woody_pixels = np.sum(main_array[valid_mask] == woody_val)
    
    # NEW: Calculate Herb_Pct and Woody_Pct
    herb_pct = (float(herb_pixels) / float(total_valid_pixels)) * 100
    woody_pct = (float(woody_pixels) / float(total_valid_pixels)) * 100
    
    if woody_pixels > 0:
        herb_woody_ratio = float(herb_pixels) / float(woody_pixels)
        ratio_display = f"{herb_woody_ratio:.4f}"
    else:
        herb_woody_ratio = None
        ratio_display = "Undefined (0 Woody Pixels)"

    # MODIFIED MEAN FETCH
    print("Calculating Mean Fetch (Including 0s)...")
    obstacle_class = Con(in_raster != target_value, 1)
    fetch_dist_raster = EucDistance(obstacle_class)
    raw_fetch_array = arcpy.RasterToNumPyArray(fetch_dist_raster, lower_left, ncols, nrows)
    forced_fetch_array = np.where(main_array != target_value, 0.0, raw_fetch_array)
    
    sum_fetch = np.sum(forced_fetch_array[valid_mask])
    mean_fetch_exact = sum_fetch / float(total_valid_pixels) 

    # ====================================================================
    # RAP CANOPY GAP FRACTIONS (ALIGNED TO MAIN SCRIPT LOGIC)
    # ====================================================================
    print("Calculating Exhaustive Canopy Gap Fractions...")
    
    def get_gap_lengths(transect_array, gap_val, p_size):
        is_gap = (transect_array == gap_val)
        padded = np.concatenate(([False], is_gap, [False]))
        diffs = np.diff(padded.astype(int))
        starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
        return (ends - starts) * p_size

    exhaust_gap_arrays = []
    masked_for_gaps = np.where(valid_mask, main_array, -9999)
    
    for r in range(nrows):
        exhaust_gap_arrays.append(get_gap_lengths(masked_for_gaps[r, :], target_value, cell_size))
    for c in range(ncols):
        exhaust_gap_arrays.append(get_gap_lengths(masked_for_gaps[:, c], target_value, cell_size))
        
    all_exhaust_gaps = np.concatenate(exhaust_gap_arrays)
    total_exhaust_length_m = total_valid_pixels * cell_size * 2
    
    ex_0_24 = (np.sum(all_exhaust_gaps[all_exhaust_gaps < 0.25]) / total_exhaust_length_m) * 100
    ex_25_50 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.25) & (all_exhaust_gaps <= 0.50)]) / total_exhaust_length_m) * 100
    ex_51_100 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 0.51) & (all_exhaust_gaps <= 1.00)]) / total_exhaust_length_m) * 100
    ex_101_200 = (np.sum(all_exhaust_gaps[(all_exhaust_gaps >= 1.01) & (all_exhaust_gaps <= 2.00)]) / total_exhaust_length_m) * 100
    ex_gt_200 = (np.sum(all_exhaust_gaps[all_exhaust_gaps > 2.00]) / total_exhaust_length_m) * 100

    # ====================================================================
    # FINAL RESULTS
    # ====================================================================
    print(f"\n--- Final Site Results ---")
    print(f"BGR: {bgr_percent:.4f}%")
    print(f"Herb Cover: {herb_pct:.4f}%")
    print(f"Woody Cover: {woody_pct:.4f}%")
    print(f"LPI: {lpi_percent:.4f}%")
    print(f"Mean Fetch: {mean_fetch_exact:.6f} m")
    print(f"Herb-to-Woody Ratio: {ratio_display}")
    
    print(f"\n--- Canopy Gap Fractions ---")
    print(f"Gap 0-24 cm:    {ex_0_24:.4f}%")
    print(f"Gap 25-50 cm:   {ex_25_50:.4f}%")
    print(f"Gap 51-100 cm:  {ex_51_100:.4f}%")
    print(f"Gap 101-200 cm: {ex_101_200:.4f}%")
    print(f"Gap > 200 cm:   {ex_gt_200:.4f}%")

    total_gaps = ex_0_24 + ex_25_50 + ex_51_100 + ex_101_200 + ex_gt_200
    print(f"-----------------------------")
    print(f"Total Gap Fraction: {total_gaps:.4f}%")

except Exception as e:
    print(f"An error occurred: {e}")

finally:
    arcpy.management.Delete("memory")
    arcpy.CheckInExtension("Spatial")
    
