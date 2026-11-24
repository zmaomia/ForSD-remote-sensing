"""
Compute tree-level forest structural diversity indices in a moving window.

For each detected tree (focal tree), a circular window with radius R_w is placed
on the tree, and only the trees inside the forest boundary are used (edge correction).

Output fields:
  Tnum       number of trees in the window (focal tree included)
  Tnum_ha    spatial tree density, SDen (trees/ha); the number of trees is divided by
             the window area inside the forest boundary
  GINI       spatial Gini index of tree height, SGINI (focal tree included)
  MLISp      spatial mingling index of tree species, SMLI (focal tree excluded)
  Deciduous  number of deciduous trees in the window (Class 3 or 4)
  Mheight    mean tree height in the window (m)
  Mcrown     mean crown area in the window (m^2)
  Edge_ratio window area inside the boundary / full window area
             (1 = not at the edge, < 1 = close to the edge)

Input point shapefile fields: "height" (tree height), "Class" (tree species),
"Shape_area" (crown area). The point and boundary files should use a projected
coordinate system in metres, e.g., ETRS-TM35FIN (EPSG:3067).
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from scipy.spatial import cKDTree

# =======================
# Settings
# =======================
shapefile_path = "./Anne_open_indices/Ann_open_1_pred.shp"  # detected trees (points)
boundary_path = "./Boundary/Ann_111.shp"                     # forest site boundary (polygons)
output_path = "./Anne_open_indices/Ann_open_1_index.shp"    # output shapefile with the indices

radius = 10  # moving-window radius R_w in metres

gdf = gpd.read_file(shapefile_path)
boundary_gdf = gpd.read_file(boundary_path)

# Use the same projected coordinate system (in metres) for points and boundary
if gdf.crs.is_geographic:
    gdf = gdf.to_crs(epsg=3067)
boundary_gdf = boundary_gdf.to_crs(gdf.crs)
boundary = boundary_gdf.unary_union  # merge multiple polygons into one boundary

coords = np.array([(point.x, point.y) for point in gdf.geometry])
tree = cKDTree(coords)

GINI_list = []
Mean_H = []
Tree_N = []
Tree_ha = []
Deciduous_N = []
crown_list = []
Edge_ratio = []  # used to correct the edge effect

# =======================
# Spatial tree density and Gini index of tree height
# =======================
for i, (x, y) in enumerate(coords):
    point_geom = Point(x, y)
    buffer_geom = point_geom.buffer(radius)

    # Intersect the window with the boundary to obtain the valid neighborhood
    valid_area = buffer_geom.intersection(boundary)
    if valid_area.is_empty:
        # The tree lies completely outside the boundary
        print('---trees are outside---')
        Tree_N.append(0)
        Tree_ha.append(0)
        Deciduous_N.append(0)
        Mean_H.append(0)
        crown_list.append(0)
        GINI_list.append(0)
        Edge_ratio.append(0)
        continue

    # Share of the window area inside the boundary
    area_ratio = valid_area.area / buffer_geom.area
    Edge_ratio.append(area_ratio)

    # Find all trees within the radius (the focal tree is included)
    indices = tree.query_ball_point([x, y], r=radius)
    # Keep only the trees inside the forest boundary
    neighbors = [idx for idx in indices if valid_area.contains(gdf.geometry.iloc[idx])]

    heights = gdf.iloc[neighbors]["height"].values
    neighbor_classes = gdf.iloc[neighbors]["Class"].astype(int).values
    crown_sizes = gdf.iloc[neighbors]["Shape_area"].values

    # Spatial tree density: divide by the window area inside the boundary,
    # instead of the full window area
    Tnum = len(heights)
    if valid_area.area > 0:
        Tnum_ha = Tnum * 10000 / valid_area.area
    else:
        Tnum_ha = 0

    # Number of deciduous trees (Class 3: birch, Class 4: aspen/other deciduous)
    Deciduous = np.sum((neighbor_classes == 3) | (neighbor_classes == 4))

    # Gini index of tree height and mean values
    # (SGINI is set to 0 if the window contains only the focal tree)
    if len(heights) > 1:
        sorted_heights = np.sort(heights)
        n = len(sorted_heights)
        mean_height = round(np.sum(heights) / n, 3)
        mean_crown = round(np.sum(crown_sizes) / n, 3)
        gini_index = round(
            (2 * np.sum((np.arange(1, n + 1) * sorted_heights)) -
             (n + 1) * np.sum(sorted_heights)) / (n * np.sum(sorted_heights)), 3)
    else:
        gini_index = 0
        mean_height = 0
        mean_crown = 0

    Tree_N.append(Tnum)
    Tree_ha.append(Tnum_ha)
    Deciduous_N.append(Deciduous)
    Mean_H.append(mean_height)
    crown_list.append(mean_crown)
    GINI_list.append(gini_index)

gdf["Tnum"] = Tree_N
gdf["Tnum_ha"] = Tree_ha
gdf["Deciduous"] = Deciduous_N
gdf["Mheight"] = Mean_H
gdf["Mcrown"] = crown_list
gdf["GINI"] = GINI_list
gdf["Edge_ratio"] = Edge_ratio   # 1 = not at the edge, < 1 = close to the edge

# =======================
# Spatial mingling index of tree species
# =======================
MLI1_list = []

for i, (x, y) in enumerate(coords):
    focal_species = gdf.iloc[i]["Class"]
    point_geom = Point(x, y)
    valid_area = point_geom.buffer(radius).intersection(boundary)

    if valid_area.is_empty:
        MLI1_list.append(0)
        continue

    # Neighbors inside the boundary, excluding the focal tree itself
    indices_m = tree.query_ball_point([x, y], r=radius)
    indices_m = [idx for idx in indices_m if idx != i and valid_area.contains(gdf.geometry.iloc[idx])]

    # Share of neighbors with a different species (SMLI is set to 0 if there are no neighbors)
    if len(indices_m) > 0:
        neighbor_species = gdf.iloc[indices_m]["Class"].values
        MLI1 = np.sum(neighbor_species != focal_species) / len(neighbor_species)
    else:
        MLI1 = 0
    MLI1_list.append(MLI1)

gdf["MLISp"] = MLI1_list

# =======================
# Output results
# =======================
gdf.to_file(output_path)
print(f"Results saved to: {output_path}")
