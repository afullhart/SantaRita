import arcpy
import os
import glob

arcpy.env.overwriteOutput = True

# ====================================================================
# USER INPUTS
# ====================================================================
# ====================================================================
# USER INPUTS
# ====================================================================

input_folder = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\GEE_Tiles"

# The folder where the final stitched TIFFs will be saved
output_folder = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data"

output_may_name = "SRER_Classified_May_2019_UTM12N_Mosaic.tif"
output_sep_name = "SRER_Classified_Sep_2019_UTM12N_Mosaic.tif"
# ====================================================================
# EXECUTION
# ====================================================================
print("Scanning folder for tiles...")

# Use glob to automatically grab all files matching the name patterns
may_tiles = glob.glob(os.path.join(input_folder, "*May_2019*.tif"))
sep_tiles = glob.glob(os.path.join(input_folder, "*Sep_2019*.tif"))

print(f"Found {len(may_tiles)} tiles for May and {len(sep_tiles)} tiles for September.\n")

# Enforce WGS 1984 UTM Zone 12N to prevent visual shifting
sr = arcpy.SpatialReference(32612) 

# Mosaic Parameters for classified drone data
pixel_type = "8_BIT_UNSIGNED" 
cell_size = 0.05
number_of_bands = 1

if may_tiles:
    print("Stitching May 2019 tiles. This might take a minute...")
    arcpy.management.MosaicToNewRaster(
        input_rasters=may_tiles,
        output_location=output_folder,
        raster_dataset_name_with_extension=output_may_name,
        coordinate_system_for_the_raster=sr,
        pixel_type=pixel_type,
        cellsize=cell_size,
        number_of_bands=number_of_bands,
        mosaic_method="LAST", 
        mosaic_colormap_mode="FIRST"  # <--- FIXED PARAMETER NAME
    )
    print("--> May Mosaic Complete!")

if sep_tiles:
    print("\nStitching September 2019 tiles. Hang tight...")
    arcpy.management.MosaicToNewRaster(
        input_rasters=sep_tiles,
        output_location=output_folder,
        raster_dataset_name_with_extension=output_sep_name,
        coordinate_system_for_the_raster=sr,
        pixel_type=pixel_type,
        cellsize=cell_size,
        number_of_bands=number_of_bands,
        mosaic_method="LAST",
        mosaic_colormap_mode="FIRST"  # <--- FIXED PARAMETER NAME
    )
    print("--> September Mosaic Complete!")

print("\nAll done! You now have two seamless master TIFFs.")
