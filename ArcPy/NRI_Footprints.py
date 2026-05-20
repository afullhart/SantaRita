import arcpy
import math

# --- CONFIGURATION ---
input_footprints = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_Grid_Metrics_Final"
output_nri_plots = r"C:\Users\andre\Documents\ArcGIS\Projects\MyProject1\MyProject1.gdb\SRER_NRI_Plots_Greedy_Strict"

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
    # This keeps the plot center far enough away that the plot boundary touches 
    # the footprint edge exactly, but does not cross it.
    safe_zone = poly_geom.buffer(-radius) 
    
    # Catch completely invalid or narrow polygons
    if safe_zone is None or safe_zone.area <= 0:
        print(f"  -> [Skipped] Footprint {poly_idx + 1} of {total_polys}: Too narrow.")
        return centers
    
    # 2. Path-Tracing: Densify the boundary of the safe zone to trace the edge.
    # Step size of 2 ensures we catch even small pockets.
    safe_path = safe_zone.boundary().densify("DISTANCE", 2)
    potential_pts = []
    for i in range(0, int(safe_path.length), 2):
        potential_pts.append(safe_path.positionAlongLine(i))
        
    # 3. Sort by distance from centroid (Furthest points first)
    # Forces the algorithm to start packing at the extreme "tips" of the footprint.
    centroid = safe_zone.centroid
    potential_pts.sort(key=lambda p: p.distanceTo(centroid), reverse=True)
    
    # 4. Greedy Packing: From outside-in
    for pt in potential_pts:
        if not centers or min([pt.distanceTo(c) for c in centers]) >= min_dist:
            centers.append(pt)
            
    print(f"  -> [Success] Footprint {poly_idx + 1} of {total_polys}: Packed {len(centers)} plots.")
    return centers

# --- MAIN EXECUTION ---
arcpy.env.overwriteOutput = True
print("Starting Strict Extremity-First Greedy Packing...")

# 1. Prepare data in memory
dissolved = arcpy.management.Dissolve(input_footprints, "memory/dissolved")
singlepart = arcpy.management.MultipartToSinglepart(dissolved, "memory/singlepart")

final_points = []
total_polys = int(arcpy.management.GetCount(singlepart).getOutput(0))
print(f"Total footprints to process: {total_polys}")

# 2. Process footprint-by-footprint
with arcpy.da.SearchCursor(singlepart, ["SHAPE@"]) as cursor:
    for i, row in enumerate(cursor):
        final_points.extend(get_greedy_centers(row[0], i, total_polys))

# 3. Export to final feature class
print("Finalizing output points...")
sr = arcpy.Describe(singlepart).spatialReference
out_pts = arcpy.management.CreateFeatureclass("memory", "greedy_pts", "POINT", spatial_reference=sr)

with arcpy.da.InsertCursor(out_pts, ["SHAPE@"]) as icursor:
    for pt in final_points:
        icursor.insertRow([pt])

# 4. Buffer back to 110m circular footprints
print("Generating final 110m plot boundaries...")
arcpy.analysis.Buffer(out_pts, output_nri_plots, f"{radius} Meters")

print(f"Success! Strict-packing plots saved to: {output_nri_plots}")
