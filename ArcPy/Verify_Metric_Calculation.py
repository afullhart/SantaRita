import arcpy
import os
import numpy as np
from arcpy.sa import *

# ====================================================================
# USER INPUT
# ====================================================================
source_tiff = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\Clipped_Classified_May_2019_UTM.tif"

# Using your exact coordinates
min_x = 515370.0
min_y = 3518619.99

# ====================================================================
# SETUP
# ====================================================================
arcpy.ResetEnvironments() 
arcpy.env.workspace = "memory"
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

print("Slicing exact 200x200 pixel matrix...")

try:
    # 1. Get exact cell size to ensure our math is perfect
    cell_size = float(arcpy.management.GetRasterProperties(source_tiff, "CELLSIZEX").getOutput(0))
    
    # A 10m box at 5cm resolution is exactly 200 x 200 pixels
    ncols = 200
    nrows = 200
    lower_left = arcpy.Point(min_x, min_y)

    # 2. Extract directly to NumPy (Bypasses ExtractByRectangle bugs entirely)
    # This forces it to grab exactly 40,000 pixels, no more, no less.
    main_array = arcpy.RasterToNumPyArray(source_tiff, lower_left, ncols, nrows)
    
    total_valid_pixels = main_array.size
    print(f"Total Valid Pixels: {total_valid_pixels}")
    
    if total_valid_pixels != 40000:
        print("!! WARNING: Did not extract exactly 40,000 pixels. Check extent. !!")
    
    # 3. Convert the perfect matrix BACK to a temporary raster for Spatial Analyst tools
    temp_raster = arcpy.NumPyArrayToRaster(main_array, lower_left, cell_size, cell_size)
    in_raster = Int(temp_raster)

    # ====================================================================
    # METRICS
    # ====================================================================
    print("Calculating Metrics...")
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

    print(f"\n--- Final Site Results ---")
    print(f"BGR: {bgr_percent:.4f}%")
    print(f"LPI: {lpi_percent:.4f}%")
    print(f"Mean Fetch: {mean_fetch_exact:.6f} m")

except Exception as e:
    print(f"An error occurred: {e}")

finally:
    arcpy.management.Delete("memory")
    arcpy.CheckInExtension("Spatial")
    
