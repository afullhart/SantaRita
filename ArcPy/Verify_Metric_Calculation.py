import arcpy
import os
from arcpy.sa import *

# 1. Your input coordinates in Decimal Degrees (Longitude, Latitude)
coords = [[[-110.8376236670315, 31.80316647982381],
           [-110.83762375798096, 31.803076268958048],
           [-110.83751816564339, 31.80307616372741],
           [-110.83751798374443, 31.803166390139502],
           [-110.8376236670315, 31.80316647982381]]]

# Access the default geodatabase of the current ArcGIS Pro project
aprx = arcpy.mp.ArcGISProject("CURRENT")
default_gdb = aprx.defaultGeodatabase

# Get the main folder where your ArcGIS Pro project lives
project_folder = aprx.homeFolder

# Set environment variables (must happen AFTER default_gdb is defined)
arcpy.env.workspace = default_gdb
arcpy.env.scratchWorkspace = default_gdb
arcpy.env.overwriteOutput = True

# Ensure Spatial Analyst is checked out
arcpy.CheckOutExtension("Spatial")

# ====================================================================
# PART 1: RASTER PROCESSING & AREA CALCULATION
# ====================================================================

# 1. Define inputs
input_raster_path = arcpy.management.Clip(
  in_raster='Clipped_Classified_May_2019.tif',
  rectangle='-110.837624030617 31.8030778206418 -110.837518984704 31.8031675281144 GEOGCRS["GCS_WGS_1984",DYNAMIC[FRAMEEPOCH[1990.5],MODEL["AM0-2"]],DATUM["D_WGS_1984",ELLIPSOID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],CS[ellipsoidal,2],AXIS["Latitude (lat)",north,ORDER[1]],AXIS["Longitude (lon)",east,ORDER[2]],ANGLEUNIT["Degree",0.0174532925199433],ID["EPSG","4326"]]',
  out_raster=os.path.join(default_gdb, 'Clipped_GridCell'),
  in_template_dataset='SRER_BareGroundPatch',
  nodata_value='256',
  clipping_geometry='NONE',
  maintain_clipping_extent='NO_MAINTAIN_EXTENT'
)

target_value = 3  # Target class (e.g., bare ground)

# Define the original coordinate system as WGS 1984 (EPSG 4326)
wgs84_sr = arcpy.SpatialReference(4326)
arcpy.management.DefineProjection(
    in_dataset=input_raster_path, 
    coor_system=wgs84_sr
)

# ===================

# 2. Define the target projection: NAD 1983 UTM Zone 12N (EPSG 26912)
utm12n_sr = arcpy.SpatialReference(26912)

# 3. Project the raster to UTM Zone 12N so cell size is measured in meters
projected_raster_path = os.path.join(default_gdb, "Clipped_GridCell_UTM")

# Note: Using nearest neighbor ("NEAREST") resampling is critical for classified/categorical data!
arcpy.management.ProjectRaster(
    in_raster=input_raster_path, 
    out_raster=projected_raster_path, 
    out_coor_system=utm12n_sr, 
    resampling_type="NEAREST"
)

# 4. Load the newly projected raster as a Spatial Analyst object
in_raster_proj = Raster(projected_raster_path)

# 5. Isolate the target class 
isolated_class = Con(in_raster_proj == target_value, 1)

# 6. Run Region Group
grouped_patches = RegionGroup(isolated_class, "EIGHT", "WITHIN", "NO_LINK")

# Define the path as a .tif file outside the geodatabase (B/C of file lock issue)
output_data_dir = os.path.join(project_folder, 'DATA')
if not os.path.exists(output_data_dir):
    os.makedirs(output_data_dir)
    
output_grouped_path = os.path.join(output_data_dir, 'Grouped_Patches_UTM.tif')

# Save the grouped output
grouped_patches.save(output_grouped_path)

# 7. Extract cell dimensions (now accurately in meters)
cell_size_x = float(arcpy.management.GetRasterProperties(grouped_patches, "CELLSIZEX").getOutput(0))
cell_size_y = float(arcpy.management.GetRasterProperties(grouped_patches, "CELLSIZEY").getOutput(0))
pixel_area_sq_meters = cell_size_x * cell_size_y 

# 8. Find the largest patch using a Search Cursor
max_pixel_count = 0

with arcpy.da.SearchCursor(output_grouped_path, ["Count"]) as cursor:
    for row in cursor:
        if row[0] > max_pixel_count:
            max_pixel_count = row[0]

# 9. Calculate the final area of the largest patch
largest_patch_area = max_pixel_count * pixel_area_sq_meters

# ====================================================================
# PART 2: OVERALL PERCENTAGE CALCULATION
# ====================================================================

# Ensure the projected raster has an attribute table built so we can read counts
arcpy.management.BuildRasterAttributeTable(projected_raster_path)

total_test_pixel_count = 0
bare_ground_pixel_count = 0

# Loop through the attribute table of the clipped, projected raster to get totals
with arcpy.da.SearchCursor(projected_raster_path, ["Value", "Count"]) as cursor:
    for row in cursor:
        val = row[0]
        count = row[1]
        
        total_test_pixel_count += count
        
        # Check if the class is bare ground
        if val == target_value:
            bare_ground_pixel_count = count

# Calculate the percentage safely avoiding division by zero
if total_test_pixel_count > 0:
    bare_ground_percent = (bare_ground_pixel_count / total_test_pixel_count) * 100
else:
    bare_ground_percent = 0.0

# ====================================================================
# PART 3: MEAN FETCH CALCULATION (Distance to Nearest Obstacle)
# ====================================================================

print("Calculating Base Fetch Distance...")

# 1. Isolate the "obstacles" (everything that is NOT bare ground)
obstacle_class = Con(in_raster_proj != target_value, 1)

# 2. Calculate Euclidean Distance to the nearest obstacle (distance across bare ground)
fetch_distance_raster = EucDistance(obstacle_class)
fetch_dist_path = os.path.join(output_data_dir, 'Fetch_Distance_UTM.tif')
fetch_distance_raster.save(fetch_dist_path)

# --------------------------------------------------------------------
# METHOD A: MAP ALGEBRA (Exact Mean of Non-Zero Pixels)
# --------------------------------------------------------------------
print("Calculating Exact Mean Fetch (Map Algebra)...")

# Isolate non-zero pixels using Map Algebra (sets 0 to NoData)
non_zero_fetch_raster = Con(fetch_distance_raster > 0, fetch_distance_raster)
non_zero_fetch_path = os.path.join(output_data_dir, 'NonZero_Fetch_UTM.tif')
non_zero_fetch_raster.save(non_zero_fetch_path)

# Calculate statistics and extract the exact Mean across all valid pixels
arcpy.management.CalculateStatistics(non_zero_fetch_path)
mean_fetch_exact = float(arcpy.management.GetRasterProperties(non_zero_fetch_path, "MEAN").getOutput(0))


# --------------------------------------------------------------------
# METHOD B: RANDOM SAMPLING (1000 Points)
# --------------------------------------------------------------------
print("Calculating Sampled Mean Fetch (1000 Random Points)...")

# Create 1000 random points within the boundaries of the test polygon
random_points_fc = "Fetch_RandomPoints"
arcpy.env.randomGenerator = "1234 ACM599"
arcpy.management.CreateRandomPoints(
    out_path=default_gdb,
    out_name=random_points_fc,
    constraining_feature_class="SRER_BareGroundPatch",
    number_of_points_or_field=1000
)

# Extract the fetch distance values to those random points
extracted_points_fc = os.path.join(default_gdb, "Extracted_Fetch_Points")
ExtractValuesToPoints(
    in_point_features=os.path.join(default_gdb, random_points_fc),
    in_raster=fetch_distance_raster,
    out_point_features=extracted_points_fc,
    interpolate_values="NONE",
    add_attributes="VALUE_ONLY"
)

# Calculate the average of those 1000 points
total_fetch_sampled = 0.0
valid_points_sampled = 0

with arcpy.da.SearchCursor(extracted_points_fc, ["RASTERVALU"]) as cursor:
    for row in cursor:
        val = row[0]
        # ---> THE FIX: Mathematically align with Map Algebra by ignoring 0 values
        if val is not None and val > 0:
            total_fetch_sampled += val
            valid_points_sampled += 1

mean_fetch_sampled = (total_fetch_sampled / valid_points_sampled) if valid_points_sampled > 0 else 0.0


# ====================================================================
# FINAL OUTPUTS
# ====================================================================
print(f"\n--- Analysis Results ---")
print(f"Projected Grid cell size: {cell_size_x:.2f}m x {cell_size_y:.2f}m")
print(f"Largest patch area (LPI): {largest_patch_area:.2f} square meters")
print(f"Overall Bare Ground Percentage: {bare_ground_percent:.2f}%")
print(f"Mean Fetch (Exact - Map Algebra): {mean_fetch_exact:.2f} meters")
print(f"Mean Fetch (Sampled - 1000 pts): {mean_fetch_sampled:.2f} meters")

# Return extension
arcpy.CheckInExtension("Spatial")
