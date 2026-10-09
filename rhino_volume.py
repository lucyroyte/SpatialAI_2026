"""Sum the volume of buildings in a Rhino .3dm file, per layer, without Rhino.

Usage:
    pip install rhino3dm
    python rhino_volume.py path/to/model.3dm [more.3dm ...] [--csv volumes_by_layer.csv]

How each object is measured:
  * Extrusion (e.g. a building footprint pushed up): profile area x height, exact.
  * Brep / polysurface with flat faces: exact, from its edges.
  * Mesh, or Brep with curved faces: from the render mesh saved in the file.

Open-surface models (like the NYC DCP 3D building model, where every building
is loose facade, roof and footprint surfaces rather than a closed solid):
the walls are vertical, so a building's volume is the area under its roofs
minus the area under its footprints, each weighted by height:
    volume = sum over roofs of (plan area x roof height)
           - sum over footprints of (plan area x ground height)
Layers whose name contains --roof-layer / --ground-layer are measured this way.

Anything else that can't be measured (other open surfaces, curves, text...)
is counted as "skipped" and reported.
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


def curve_points(crv, samples=16):
    """Vertices of a polyline-like curve; a curved one is sampled at `samples` segments."""
    pl = crv.TryGetPolyline()
    if pl is not None:
        return [xyz(pl[i]) for i in range(len(pl))]
    t0, t1 = crv.Domain.T0, crv.Domain.T1
    return [xyz(crv.PointAt(t0 + (t1 - t0) * i / samples)) for i in range(samples + 1)]


def extrusion_volume(ext):
    """Profile area x height. Returns None if a profile isn't straight-sided."""
    if not ext.IsSolid:
        return None
    area_vec = None
    for i in range(ext.ProfileCount):
        n = newell(curve_points(ext.Profile3d(i, 0.0)))
        if area_vec is None:
            area_vec = n
        else:
            # Holes subtract, whichever direction they were drawn.
            s = -1 if dot(n, area_vec) > 0 else 1
            area_vec = tuple(area_vec[k] + s * n[k] for k in range(3))
    height_vec = sub(xyz(ext.PathEnd), xyz(ext.PathStart))
    return abs(dot(area_vec, height_vec)) / 2.0


def face_loops(brep, face):
    """Boundary loops of a Brep face as 3D point lists, outer loop first."""
    loops = []
    for loop in face.Loops:
        pts = []
        for trim in loop.Trims:
            if trim.EdgeIndex < 0:
                continue
            edge_pts = curve_points(brep.Edges[trim.EdgeIndex])
            if trim.IsReversed:
                edge_pts = edge_pts[::-1]
            pts.extend(edge_pts[:-1])
        if len(pts) >= 3:
            loops.append(pts)
    return loops


def height_integral(loops):
    """Plan area of a flat face and the integral of its height over that area.

    Holes (every loop after the first) subtract. Exact for planar faces, since
    height varies linearly across them; a fan of triangles handles any polygon.
    """
    area = z_area = 0.0
    for k, p in enumerate(loops):
        a = za = 0.0
        x0, y0, z0 = p[0]
        for i in range(1, len(p) - 1):
            (x1, y1, z1), (x2, y2, z2) = p[i], p[i + 1]
            t = ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) / 2.0
            a += t
            za += t * (z0 + z1 + z2) / 3.0
        sign = (1 if k == 0 else -1) * (1 if a >= 0 else -1)
        area += sign * a
        z_area += sign * za
    return area, z_area


def surface_height_integral(geom):
    """Sum of height_integral over every face of an open Brep (or None)."""
    if not isinstance(geom, r3.Brep):
        return None
    area = z_area = 0.0
    for fi in range(len(geom.Faces)):
        a, za = height_integral(face_loops(geom, geom.Faces[fi]))
        area += a
        z_area += za
    return area, z_area


def brep_volume(brep):
    """Exact volume of a closed Brep whose faces are all flat with straight edges."""
    if not brep.IsSolid:
        return None
    total = 0.0
    for fi in range(len(brep.Faces)):
        face = brep.Faces[fi]
        if not face.IsPlanar(1e-3):
            return None
        for pts in face_loops(brep, face):
            n = newell(pts)
            if face.OrientationIsReversed:
                n = (-n[0], -n[1], -n[2])
            total += dot(pts[0], n)
    return abs(total) / 6.0


def mesh_volume(mesh, origin):
    """Signed volume, measured from `origin` (near the mesh, so far-off map
    coordinates don't swamp the result in rounding error)."""
    total = 0.0
    v = mesh.Vertices
    for i in range(len(mesh.Faces)):
        f = mesh.Faces[i]
        a, b, c = (sub(xyz(v[f[k]]), origin) for k in range(3))
        total += dot(a, cross(b, c))
        if f[2] != f[3]:
            d = sub(xyz(v[f[3]]), origin)
            total += dot(a, cross(c, d))
    return total / 6.0


def object_volume(geom):
    # (mesh, sign): a Brep face's render mesh follows its surface, which can
    # point the opposite way to the face itself.
    if isinstance(geom, r3.Extrusion):
        vol = extrusion_volume(geom)
        if vol is not None:
            return vol
        meshes = [(geom.GetMesh(r3.MeshType.Any), 1)]
    elif isinstance(geom, r3.Brep):
        vol = brep_volume(geom)
        if vol is not None:
            return vol
        if not geom.IsSolid:
            return None
        meshes = [(geom.Faces[i].GetMesh(r3.MeshType.Any), -1 if geom.Faces[i].OrientationIsReversed else 1)
                  for i in range(len(geom.Faces))]
    elif isinstance(geom, r3.Mesh):
        if not geom.IsClosed:
            return None
        meshes = [(geom, 1)]
    else:
        return None
    if not meshes or any(m is None for m, _ in meshes):
        return None
    origin = xyz(geom.GetBoundingBox().Center)
    return abs(sum(s * mesh_volume(m, origin) for m, s in meshes))


def measure(path, roof_key, ground_key):
    """Per-layer volumes (model units) for one file, plus the model's unit name."""
    print(f"Reading {path} (this can take a few minutes for a big file)...", flush=True)
    model = r3.File3dm.Read(path)
    if model is None:
        sys.exit(f"Could not read {path}. Is it a .3dm file, and is the path right?")

    unit = str(model.Settings.ModelUnitSystem).split(".")[-1]
    layers = {model.Layers[i].Index: model.Layers[i].FullPath for i in range(len(model.Layers))}

    vol = defaultdict(float)
    counted = defaultdict(int)
    skipped = defaultdict(int)
    plan = {"roof": [0.0, 0.0], "ground": [0.0, 0.0]}  # [plan area, height x area]
    n = len(model.Objects)
    for i, obj in enumerate(model.Objects):
        layer = layers.get(obj.Attributes.LayerIndex, "(unknown layer)")
        geom = obj.Geometry
        v = object_volume(geom)
        part = "roof" if roof_key in layer else "ground" if ground_key in layer else None
        if v is None and part:
            hz = surface_height_integral(geom)
            if hz is not None:
                plan[part][0] += hz[0]
                plan[part][1] += hz[1]
                counted[layer] += 1
                continue
        if v is None:
            skipped[layer] += 1
        else:
            vol[layer] += v
            counted[layer] += 1
        if (i + 1) % 50000 == 0:
            print(f"  {i + 1:,} / {n:,} objects", flush=True)

    rows = {layer: [vol[layer], counted[layer], skipped[layer]] for layer in set(vol) | set(counted) | set(skipped)}
    (roof_area, roof_z), (ground_area, ground_z) = plan["roof"], plan["ground"]
    if roof_area:
        name = f"(open surfaces: '{roof_key}' minus '{ground_key}')"
        rows[name] = [roof_z - ground_z, 0, 0]
        print(f"  Open-surface buildings: roof plan area {roof_area:,.0f}, footprint plan area {ground_area:,.0f} "
              f"square {unit.lower()}")
        if ground_area:
            print(f"  (these two should be close)  mean height above ground: "
                  f"{(roof_z - ground_z) / ground_area:,.1f} {unit.lower()}")
    return unit, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("models", nargs="+", help="path(s) to .3dm files")
    ap.add_argument("--csv", default="volumes_by_layer.csv", help="where to write the per-layer table")
    ap.add_argument("--roof-layer", default="RoofTop Surface",
                    help="layers containing this text hold open roof surfaces (default: %(default)s)")
    ap.add_argument("--ground-layer", default="FootPrint Surface",
                    help="layers containing this text hold open footprint surfaces (default: %(default)s)")
    args = ap.parse_args()

    grand = 0.0
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "layer", "unit", "volume_model_units3", "volume_m3", "objects_measured", "objects_skipped"])
        for path in args.models:
            unit, rows = measure(path, args.roof_layer, args.ground_layer)
            to_m3 = TO_METRES.get(unit, 1.0) ** 3
            print(f"\n{path}  (model units: {unit})")
            print(f"{'layer':55s} {'volume (m3)':>16s} {'measured':>9s} {'skipped':>8s}")
            for layer in sorted(rows):
                v, c, k = rows[layer]
                w.writerow([path, layer, unit, round(v, 2), round(v * to_m3, 2), c, k])
                print(f"{layer[:55]:55s} {v * to_m3:16,.0f} {c:9,d} {k:8,d}")
            total = sum(r[0] for r in rows.values()) * to_m3
            grand += total
            print(f"TOTAL for this file: {total:,.0f} m3")
            print(f"Measured {sum(r[1] for r in rows.values()):,} objects, "
                  f"skipped {sum(r[2] for r in rows.values()):,}.")
    if len(args.models) > 1:
        print(f"\nGRAND TOTAL ({len(args.models)} files): {grand:,.0f} m3")
    print(f"Per-layer table written to {args.csv}")


if __name__ == "__main__":
    main()
