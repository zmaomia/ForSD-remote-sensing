# ForSD-remote-sensing

This repository provides supplementary materials, visualizations, and experimental results for our research on boreal forest structural diversity charaterrization. 
Project page with results and figures: https://zmaomia.github.io/ForSD-remote-sensing/

## Calculating the structural diversity indices

`compute_indices.py` calculates three tree-level structural diversity indices from individual tree detection (ITD) results, using a moving window with radius $R_w$ centered on each tree:

- **SDen**: spatial tree density (trees/ha)
- **SGINI**: spatial Gini index of tree height
- **SMLI**: spatial mingling index of tree species

Only the trees inside the forest boundary are used, and SDen is calculated using the area of the window inside the boundary (edge-effect correction).

### Input

- A point shapefile of the detected trees with the fields `height` (tree height in m), `Class` (tree species code), and `Shape_area` (crown area in m²).
- A polygon shapefile of the forest site boundary.

### Output

The output shapefile contains the original fields and the following new fields:

| Field | Description |
|---|---|
| `Tnum` | Number of trees in the window (focal tree included) |
| `Tnum_ha` | SDen, spatial tree density (trees/ha) |
| `GINI` | SGINI, spatial Gini index of tree height (focal tree included) |
| `MLISp` | SMLI, spatial mingling index of tree species (focal tree excluded) |
| `Deciduous` | Number of deciduous trees in the window (`Class` Birch (3) or Aspen (4)) |
| `Mheight` | Mean tree height in the window (m) |
| `Mcrown` | Mean crown area in the window (m²) |
| `Edge_ratio` | Window area inside the boundary / full window area (1 = not at the edge) |

If the window contains only the focal tree, SGINI and SMLI are set to 0.

## Data availability

The data and trained models are available from the corresponding author upon reasonable request.


