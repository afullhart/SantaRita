import arcpy
import os
import numpy as np
from arcpy.sa import *

# ====================================================================
# USER INPUT
# ====================================================================
source_tiff = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Classified_May_2019_UTM12N_Mosaic.tif"
raw_shapefile = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\SRER_Model_Grid_Shp\SRER_Model_Grid_Export.shp" 
grid_fc = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_Grid_Metrics_Final"

target_oid = 4061 

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

    # 3. Extract directly to NumPy
    main_array = arcpy.RasterToNumPyArray(source_tiff, lower_left, ncols, nrows)
    
    total_valid_pixels = main_array.size
    print(f"Total Valid Pixels: {total_valid_pixels}")
    
    if total_valid_pixels != 40000:
        print("!! WARNING: Did not extract exactly 40,000 pixels. The raster may not fully cover this grid cell. !!")
    else:
        print("SUCCESS: Exactly 40,000 pixels extracted (10m x 10m box).")
    
    # 4. Convert BACK to a temporary raster
    temp_raster = arcpy.NumPyArrayToRaster(main_array, lower_left, cell_size, cell_size)
    in_raster = Int(temp_raster)

    # ====================================================================
    # CORE METRICS (BGR, LPI, FETCH)
    # ====================================================================
    print("\nCalculating Core Metrics...")
    target_value = 3 
    
    # BGR
    bare_ground_pixels = np.sum(main_array == target_value)
    bgr_percent = (float(bare_ground_pixels) / float(total_valid_pixels)) * 100

    # LPI
    isolated_class = Con(in_raster == target_value, 1)
    grouped_patches = RegionGroup(isolated_class, "EIGHT", "WITHIN", "NO_LINK")
    patch_array = arcpy.RasterToNumPyArray(grouped_patches)
    
    patch_ids, counts = np.unique(patch_array, return_counts=True)
    valid_patch_mask = patch_ids > 0 
    max_patch_pixels = np.max(counts[valid_patch_mask]) if np.any(valid_patch_mask) else 0
    
    lpi_percent = (float(max_patch_pixels) / float(total_valid_pixels)) * 100

    # Mean Fetch
    obstacle_class = Con(in_raster != target_value, 1)
    fetch_dist_raster = EucDistance(obstacle_class)
    fetch_bare_only = Con(in_raster == target_value, fetch_dist_raster)
    
    fetch_array = arcpy.RasterToNumPyArray(fetch_bare_only)
    valid_fetch = fetch_array[fetch_array > 0]
    mean_fetch_exact = np.mean(valid_fetch) if valid_fetch.size > 0 else 0.0

    # ====================================================================
    # RAP CANOPY GAP FRACTIONS
    # ====================================================================
    print("Calculating Canopy Gap Fractions...")
    
    transect_indices = np.arange(10, 200, 20)
    
    horizontal_transects = [main_array[i, :] for i in transect_indices]
    vertical_transects = [main_array[:, j] for j in transect_indices]

    all_transects = horizontal_transects + vertical_transects
    total_transect_length_m = len(all_transects) * (ncols * cell_size) 

    def get_gap_lengths(transect_array, gap_val, p_size):
        is_gap = (transect_array == gap_val)
        padded = np.concatenate(([False], is_gap, [False]))
        diffs = np.diff(padded.astype(int))
        starts = np.where(diffs == 1)[0]
        ends = np.where(diffs == -1)[0]
        return (ends - starts) * p_size

    all_gap_lengths = []
    for t in all_transects:
        gaps = get_gap_lengths(t, target_value, cell_size)
        all_gap_lengths.extend(gaps)

    all_gap_lengths = np.array(all_gap_lengths)

    class_0_24    = all_gap_lengths[(all_gap_lengths < 0.25)]
    class_25_50   = all_gap_lengths[(all_gap_lengths >= 0.25) & (all_gap_lengths <= 0.50)]
    class_51_100  = all_gap_lengths[(all_gap_lengths >= 0.51) & (all_gap_lengths <= 1.00)]
    class_101_200 = all_gap_lengths[(all_gap_lengths >= 1.01) & (all_gap_lengths <= 2.00)]
    class_gt_200  = all_gap_lengths[(all_gap_lengths > 2.00)]

    fraction_0_24    = (np.sum(class_0_24) / total_transect_length_m) * 100
    fraction_25_50   = (np.sum(class_25_50) / total_transect_length_m) * 100
    fraction_51_100  = (np.sum(class_51_100) / total_transect_length_m) * 100
    fraction_101_200 = (np.sum(class_101_200) / total_transect_length_m) * 100
    fraction_gt_200  = (np.sum(class_gt_200) / total_transect_length_m) * 100

    # ====================================================================
    # FINAL RESULTS
    # ====================================================================
    print(f"\n--- Final Site Results ---")
    print(f"BGR: {bgr_percent:.4f}%")
    print(f"LPI: {lpi_percent:.4f}%")
    print(f"Mean Fetch: {mean_fetch_exact:.6f} m")
    
    print(f"\n--- Canopy Gap Fractions ---")
    print(f"Gap 0-24 cm:    {fraction_0_24:.4f}%")
    print(f"Gap 25-50 cm:   {fraction_25_50:.4f}%")
    print(f"Gap 51-100 cm:  {fraction_51_100:.4f}%")
    print(f"Gap 101-200 cm: {fraction_101_200:.4f}%")
    print(f"Gap > 200 cm:   {fraction_gt_200:.4f}%")

    total_gaps = fraction_0_24 + fraction_25_50 + fraction_51_100 + fraction_101_200 + fraction_gt_200
    print(f"-----------------------------")
    print(f"Total Gap Fraction: {total_gaps:.4f}%")

except Exception as e:
    print(f"An error occurred: {e}")

finally:
    arcpy.management.Delete("memory")
    arcpy.CheckInExtension("Spatial")
