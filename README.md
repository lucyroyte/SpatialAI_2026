# SpatialAI_2026

## Measuring building volume in a Rhino model

`rhino_volume.py` adds up the volume of the objects in a Rhino `.3dm` file, per layer,
without Rhino installed.

**You need** Python 3.9 or newer, `rhino3dm` and `numpy`:

```
pip install -r requirements.txt
```

**Run it** with the path to your `.3dm` file:

```
python rhino_volume.py path/to/model.3dm
```

It prints a table of volume per layer (in cubic metres, converted from the model's
units) and writes the same table to `volumes_by_layer.csv`. Options:

| Option | What it does |
| --- | --- |
| `--csv FILE` | Write the per-layer table somewhere else. |
| `--top-level` | Add sub-layers into their top-level layer, so `Buildings::Facade` and `Buildings::Roof` show as one `Buildings` row. |
| `--roof-layer LAYER --base-layer LAYER` | Roof / base mode, for models whose buildings are loose surfaces instead of closed solids (see below). |
| `--cell SIZE` | Grid size for roof / base mode, in model units. Default 1 ft. Use a smaller one for a small model. |

The per-layer table only measures closed objects: extrusions, closed polysurfaces and
closed meshes. Open surfaces, curves and text are listed as "skipped".

### NYC 3D building model

The NYC DCP 3D model (`NYC_3DModel_MN01.3dm` and so on) stores each building as separate
roof, facade and footprint surfaces, so almost everything is "skipped" in the per-layer
table. Use roof / base mode instead:

```
python rhino_volume.py NYC_3DModel_MN01.3dm --top-level \
    --roof-layer "Buildings::RoofTop Surface" \
    --base-layer "Buildings::FootPrint Surface"
```

This cuts the plan into a 1 ft grid and adds up, cell by cell, the height from the
footprint to the highest roof above it. Facades aren't needed. Where roof surfaces of
one building overlap in plan (about 60 buildings in MN01 to MN03), only the highest
counts, and a footprint stored twice (4 in MN01) counts once. The script also prints
what adding up every roof without that correction would give, and how much roof or
footprint area has nothing above or below it.

Results (matching the separate per-building calculation in this project to within 0.01%):

| District | Building volume |
| --- | --- |
| MN01 | 66.2M m³ |
| MN02 | 39.3M m³ |
| MN03 | 26.6M m³ |

A 370 to 560 MB district file takes one to two minutes.
