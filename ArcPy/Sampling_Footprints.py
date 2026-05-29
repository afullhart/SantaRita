import arcpy
import math
import os

# --- CONFIGURATION ---
input_footprints = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\Data\srer_drone_veg_monitoring_footprints_2019\srer_drone_veg_monitoring_footprints_2019.shp'
out_gdb = r'C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb'

# Output paths
output_nri_plots = os.path.join(out_gdb, 'SRER_NRI_Plots_110m')
output_grid_10m = os.path.join(out_gdb, 'SRER_Grid_10m')
output_grid_30m = os.path.join(out_gdb, 'SRER_Grid_30m')

# NRI plot diameter is 110m (55m radius)
radius = 55.0
min_dist = 110.0

def get_greedy_centers(poly_geom, poly_idx, total_polys):
    """
    Extremity-First Greedy Algorithm with STRICT Boundaries.
    Guarantees no bleeding while maximizing edge-placement.
    """
    centers = []

    # 1. STRICT BUFFER: Buffer by full radius to prevent any boundary bleeding.
    safe_zone = poly_geom.buffer(-radius)

    # Catch completely invalid or narrow polygons
    if safe_zone is None or safe_zone.area <= 0:
        print(f'  -> [Skipped] Footprint {poly_idx + 1} of {total_polys}: Too narrow to fit a 110m plot.')
        return centers

    # 2. Path-Tracing: Get the boundary
    safe_path = safe_zone.boundary()
    
    # Generate points along the line every 2 meters.
    potential_pts = []
    for i in range(0, int(safe_path.length), 2):
        potential_pts.append(safe_path.positionAlongLine(i))

    # 3. Sort by distance from centroid (Furthest points first)
    centroid = safe_zone.centroid
    potential_pts.sort(key=lambda p: p.distanceTo(centroid), reverse=True)

    # 4. Greedy Packing: From outside-in
    for pt in potential_pts:
        if not centers or min([pt.distanceTo(c) for c in centers]) >= min_dist:
            centers.append(pt)

    print(f'  -> [Success] Footprint {poly_idx + 1} of {total_polys}: Packed {len(centers)} plots.')
    return centers


def generate_intersecting_grid(target_poly, cell_size, out_feature_class):
    """
    Generates a grid of perfect squares over the extent of the target polygon,
    then saves ONLY the squares that are completely inside the target polygon.
    """
    print(f'Generating {cell_size}m x {cell_size}m grid...')
    
    # Get Extent of the working polygons to frame the fishnet
    desc_poly = arcpy.Describe(target_poly)
    ext = desc_poly.extent
    sr = desc_poly.spatialReference
    
    origin_coord = f"{ext.XMin} {ext.YMin}"
    y_axis_coord = f"{ext.XMin} {ext.YMax}"
    corner_coord = f"{ext.XMax} {ext.YMax}"
    
    temp_fishnet = rf'memory\temp_fishnet_{cell_size}m'
    
    # Force the tool to use the correct spatial reference during creation
    with arcpy.EnvManager(outputCoordinateSystem=sr):
        arcpy.management.CreateFishnet(
            out_feature_class=temp_fishnet,
            origin_coord=origin_coord,
            y_axis_coord=y_axis_coord,
            cell_width=cell_size,
            cell_height=cell_size,
            number_rows=0,
            number_columns=0,
            corner_coord=corner_coord,
            labels="NO_LABELS",
            template=target_poly,
            geometry_type="POLYGON"
        )
    
    print('  -> Selecting only squares completely within the footprint...')
    # Select only the grid cells that are COMPLETELY WITHIN the footprint
    fishnet_lyr = arcpy.management.MakeFeatureLayer(temp_fishnet, f"lyr_{cell_size}m")
    arcpy.management.SelectLayerByLocation(
        in_layer=fishnet_lyr,
        overlap_type="COMPLETELY_WITHIN",  # <--- The crucial change
        select_features=target_poly,
        selection_type="NEW_SELECTION"
    )
    
    # Export selected grid cells to the final database
    arcpy.management.CopyFeatures(fishnet_lyr, out_feature_class)
    print(f'  -> Saved completely contained {cell_size}m grid to: {out_feature_class}')
    
    # Cleanup memory
    arcpy.management.Delete(temp_fishnet)
    arcpy.management.Delete(fishnet_lyr)

# --- MAIN EXECUTION ---
arcpy.env.overwriteOutput = True

print('Checking Spatial Reference...')
desc = arcpy.Describe(input_footprints)
working_footprints = input_footprints

# --- AUTO-PROJECT TO METERS ---
if desc.spatialReference.type == 'Geographic':
    print(f"Notice: Input is in '{desc.spatialReference.name}' (Degrees).")
    print("Auto-projecting to WGS 1984 UTM Zone 12N (Meters) in memory...")
    
    # EPSG 32612 is WGS 84 / UTM zone 12N (Standard for Southern Arizona)
    metric_sr = arcpy.SpatialReference(32612) 
    working_footprints = r'memory\projected_footprints'
    arcpy.management.Project(input_footprints, working_footprints, metric_sr)
    
    # Update our spatial reference description to the new metric one
    desc = arcpy.Describe(working_footprints)

# ---------------------------------------------------------
# PART 1: 110m CIRCULAR PLOTS (Original Logic)
# ---------------------------------------------------------
print('\n--- Starting Strict Extremity-First Greedy Packing ---')

# Prepare data in memory 
dissolved = arcpy.management.Dissolve(working_footprints, r'memory\dissolved')
singlepart = arcpy.management.MultipartToSinglepart(dissolved, r'memory\singlepart')

final_points = []
total_polys = int(arcpy.management.GetCount(singlepart).getOutput(0))
print(f'Total footprints to process: {total_polys}')

# Process footprint-by-footprint
with arcpy.da.SearchCursor(singlepart, ['SHAPE@']) as cursor:
    for i, row in enumerate(cursor):
        final_points.extend(get_greedy_centers(row[0], i, total_polys))

# Export to final feature class
if not final_points:
    print("WARNING: No valid 110m plots could be placed in any footprint.")
else:
    print('Finalizing output points...')
    sr = desc.spatialReference
    out_pts = arcpy.management.CreateFeatureclass('memory', 'greedy_pts', 'POINT', spatial_reference=sr)

    with arcpy.da.InsertCursor(out_pts, ['SHAPE@']) as icursor:
        for pt in final_points:
            icursor.insertRow([pt])

    # Buffer back to 110m circular footprints
    print('Generating final 110m plot boundaries...')
    arcpy.analysis.Buffer(out_pts, output_nri_plots, f'{radius} Meters')
    print(f'Success! Strict-packing plots saved to: {output_nri_plots}')

# ---------------------------------------------------------
# PART 2: 10m and 30m OVERLAPPING GRIDS
# ---------------------------------------------------------
print('\n--- Generating Grid Footprints ---')

# Use the dissolved footprint layer as the target to easily capture overlaps
generate_intersecting_grid(dissolved, 10, output_grid_10m)
generate_intersecting_grid(dissolved, 30, output_grid_30m)

print('\nAll processing complete!')
