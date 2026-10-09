# SpatialAI_2026

## Measuring building volume in a Rhino model

`rhino_volume.py` adds up the volume of the objects in a Rhino `.3dm` file, per layer,
without Rhino installed.

**You need** Python 3.9 or newer and `rhino3dm`:

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

This measures the space between each roof and the ground under it: every roof surface
is a prism down to z=0, every footprint is a prism down to z=0, and the volume is the
difference. Facades aren't needed. The script also prints the roof and footprint areas
seen from above; they should be close, and it warns if they differ by more than 5%.
It assumes no building has a roof overhanging another roof (true for this model, where
the two areas agree to within 2%).

A 370 MB borough district file takes about 30 seconds.
