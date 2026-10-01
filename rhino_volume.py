"""Sum the volume of closed objects in a Rhino .3dm file, per layer, without Rhino.

Usage:
    pip install rhino3dm
    python rhino_volume.py path/to/model.3dm [--csv volumes_by_layer.csv]

How each object is measured:
  * Extrusion (e.g. a building footprint pushed up): profile area x height, exact.
  * Brep / polysurface with flat faces: exact, from its edges.
  * Mesh, or Brep with curved faces: from the render mesh saved in the file.
Objects that can't be measured (open surfaces, curved Breps saved without
render meshes, curves, text...) are counted as "skipped" and reported.
"""
import argparse
import csv
import sys
from collections import defaultdict

import rhino3dm as r3

# Conversion from the model's unit to cubic metres, by rhino3dm UnitSystem name.
TO_METRES = {
    "Millimeters": 0.001, "Centimeters": 0.01, "Meters": 1.0, "Kilometers": 1000.0,
    "Inches": 0.0254, "Feet": 0.3048, "Yards": 0.9144, "Miles": 1609.344,
}


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def xyz(p):
    return (p.X, p.Y, p.Z)


def newell(pts):
    """Twice the area-weighted normal of a closed polygon (sum of p_i x p_i+1)."""
    n = [0.0, 0.0, 0.0]
    for i in range(len(pts)):
        c = cross(pts[i], pts[(i + 1) % len(pts)])
        n[0] += c[0]; n[1] += c[1]; n[2] += c[2]
    return tuple(n)


def curve_points(crv):
    """Vertices of a polyline-like curve, or None if it has curved segments."""
    pl = crv.TryGetPolyline()
    if pl is None:
        return None
    return [xyz(pl[i]) for i in range(len(pl))]


def extrusion_volume(ext):
    """Profile area x height. Returns None if a profile isn't straight-sided."""
    if not ext.IsSolid:
        return None
    area_vec = None
    for i in range(ext.ProfileCount):
        pts = curve_points(ext.Profile3d(i, 0.0))
        if pts is None:
            return None
        n = newell(pts)
        if area_vec is None:
            area_vec = n
        else:
            # Holes subtract, whichever direction they were drawn.
            s = -1 if dot(n, area_vec) > 0 else 1
            area_vec = tuple(area_vec[k] + s * n[k] for k in range(3))
    height_vec = sub(xyz(ext.PathEnd), xyz(ext.PathStart))
    return abs(dot(area_vec, height_vec)) / 2.0


def brep_volume(brep):
    """Exact volume of a closed Brep whose faces are all flat with straight edges."""
    if not brep.IsSolid:
        return None
    total = 0.0
    for fi in range(len(brep.Faces)):
        face = brep.Faces[fi]
        if not face.IsPlanar(1e-3):
            return None
        for loop in face.Loops:
            pts = []
            for trim in loop.Trims:
                if trim.EdgeIndex < 0:
                    continue
                edge_pts = curve_points(brep.Edges[trim.EdgeIndex])
                if edge_pts is None:
                    return None
                if trim.IsReversed:
                    edge_pts = edge_pts[::-1]
                pts.extend(edge_pts[:-1])
            if len(pts) < 3:
                continue
            n = newell(pts)
            if face.OrientationIsReversed:
                n = (-n[0], -n[1], -n[2])
            total += dot(pts[0], n)
    return abs(total) / 6.0


def mesh_volume(mesh):
    total = 0.0
    v = mesh.Vertices
    for i in range(len(mesh.Faces)):
        f = mesh.Faces[i]
        a, b, c = xyz(v[f[0]]), xyz(v[f[1]]), xyz(v[f[2]])
        total += dot(a, cross(b, c))
        if f[2] != f[3]:
            d = xyz(v[f[3]])
            total += dot(a, cross(c, d))
    return total / 6.0


def object_volume(geom):
    if isinstance(geom, r3.Extrusion):
        vol = extrusion_volume(geom)
        if vol is not None:
            return vol
        meshes = [geom.GetMesh(r3.MeshType.Any)]
    elif isinstance(geom, r3.Brep):
        vol = brep_volume(geom)
        if vol is not None:
            return vol
        if not geom.IsSolid:
            return None
        meshes = [geom.Faces[i].GetMesh(r3.MeshType.Any) for i in range(len(geom.Faces))]
    elif isinstance(geom, r3.Mesh):
        if not geom.IsClosed:
            return None
        meshes = [geom]
    else:
        return None
    if not meshes or any(m is None for m in meshes):
        return None
    return abs(sum(mesh_volume(m) for m in meshes))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="path to the .3dm file")
    ap.add_argument("--csv", default="volumes_by_layer.csv", help="where to write the per-layer table")
    args = ap.parse_args()

    print(f"Reading {args.model} (this can take a few minutes for a big file)...", flush=True)
    model = r3.File3dm.Read(args.model)
    if model is None:
        sys.exit("Could not read that file. Is it a .3dm file, and is the path right?")

    unit = str(model.Settings.ModelUnitSystem).split(".")[-1]
    to_m3 = TO_METRES.get(unit, 1.0) ** 3
    layers = {model.Layers[i].Index: model.Layers[i].FullPath for i in range(len(model.Layers))}

    vol = defaultdict(float)
    counted = defaultdict(int)
    skipped = defaultdict(int)
    n = len(model.Objects)
    for i, obj in enumerate(model.Objects):
        layer = layers.get(obj.Attributes.LayerIndex, "(unknown layer)")
        v = object_volume(obj.Geometry)
        if v is None:
            skipped[layer] += 1
        else:
            vol[layer] += v
            counted[layer] += 1
        if (i + 1) % 50000 == 0:
            print(f"  {i + 1:,} / {n:,} objects", flush=True)

    rows = sorted(set(vol) | set(skipped))
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["layer", f"volume_{unit.lower()}3", "volume_m3", "objects_measured", "objects_skipped"])
        for layer in rows:
            w.writerow([layer, round(vol[layer], 2), round(vol[layer] * to_m3, 2), counted[layer], skipped[layer]])

    print(f"\nModel units: {unit}")
    print(f"{'layer':45s} {'volume (m3)':>18s} {'measured':>9s} {'skipped':>8s}")
    for layer in rows:
        print(f"{layer[:45]:45s} {vol[layer] * to_m3:18,.0f} {counted[layer]:9,d} {skipped[layer]:8,d}")
    total = sum(vol.values())
    print(f"\nTOTAL: {total * to_m3:,.0f} m3  ({total:,.0f} cubic {unit.lower()})")
    print(f"Measured {sum(counted.values()):,} objects, skipped {sum(skipped.values()):,}.")
    print(f"Per-layer table written to {args.csv}")


if __name__ == "__main__":
    main()
