"""Sum the volume of closed objects in a Rhino .3dm file, per layer, without Rhino.

Usage:
    pip install rhino3dm
    python rhino_volume.py path/to/model.3dm [--csv volumes_by_layer.csv] [--top-level]

    # Models where buildings are loose roof / footprint surfaces (e.g. the NYC
    # DCP 3D building model), not closed solids:
    python rhino_volume.py NYC_3DModel_MN01.3dm \
        --roof-layer "Buildings::RoofTop Surface" --base-layer "Buildings::FootPrint Surface"

How each object is measured:
  * Extrusion (e.g. a building footprint pushed up): profile area x height, exact.
  * Brep / polysurface with flat faces: exact, from its edges.
  * Mesh, or Brep with curved faces: from the render mesh saved in the file.
Objects that can't be measured (open surfaces, curved Breps saved without
render meshes, curves, text...) are counted as "skipped" and reported.

Roof / base mode (--roof-layer and --base-layer): the plan is cut into a fine
grid (1 ft by default), and each cell adds (highest roof - ground) x cell area.
Walls are not needed and it does not matter which way surfaces face. Where roof
surfaces overlap in plan, only the highest one counts, and a footprint stored
twice counts once.
"""
import argparse
import csv
import sys
from collections import defaultdict

import numpy as np
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
    centre = xyz(brep.GetBoundingBox().Center)
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
            total += dot(sub(pts[0], centre), n)
    return abs(total) / 6.0


def mesh_volume(mesh, centre):
    """Signed volume, with tetrahedra measured from a point near the object.

    Measuring from the world origin instead loses all precision at map
    coordinates (NYC state plane is ~1,000,000 ft from the origin)."""
    total = 0.0
    v = mesh.Vertices
    for i in range(len(mesh.Faces)):
        f = mesh.Faces[i]
        a, b, c = (sub(xyz(v[f[k]]), centre) for k in range(3))
        total += dot(a, cross(b, c))
        if f[2] != f[3]:
            d = sub(xyz(v[f[3]]), centre)
            total += dot(a, cross(c, d))
    return total / 6.0


def prism(pts):
    """Signed (projected area, area x mean height) of a planar polygon over z=0."""
    area = az = 0.0
    p0 = pts[0]
    for i in range(1, len(pts) - 1):
        p1, p2 = pts[i], pts[i + 1]
        t = ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0
        area += t
        az += t * (p0[2] + p1[2] + p2[2]) / 3.0
    return area, az


def mesh_prism(mesh):
    area = az = 0.0
    v = mesh.Vertices
    for i in range(len(mesh.Faces)):
        f = mesh.Faces[i]
        tris = [(f[0], f[1], f[2])] + ([(f[0], f[2], f[3])] if f[2] != f[3] else [])
        for tri in tris:
            a, b = prism([xyz(v[k]) for k in tri])
            area += a; az += b
    return area, az


def face_prism(brep, fi):
    """(projected area, area x mean height) of one Brep face, facing up, or None."""
    face = brep.Faces[fi]
    area = az = 0.0
    outer = None
    for loop in face.Loops:
        pts = []
        for trim in loop.Trims:
            if trim.EdgeIndex < 0:
                continue
            edge_pts = curve_points(brep.Edges[trim.EdgeIndex])
            if edge_pts is None:
                pts = None
                break
            if trim.IsReversed:
                edge_pts = edge_pts[::-1]
            pts.extend(edge_pts[:-1])
        if pts is None or not face.IsPlanar(1e-3):
            mesh = face.GetMesh(r3.MeshType.Any)
            if mesh is None:
                return None
            area, az = mesh_prism(mesh)
            break
        if len(pts) < 3:
            continue
        a, b = prism(pts)
        if outer is None:
            outer = 1 if a > 0 else -1
        # Holes are wound the other way round from the outer loop, so they subtract.
        area += outer * a
        az += outer * b
    s = 1 if area >= 0 else -1
    return s * area, s * az


def surfaces_prism(geom):
    """(projected area, area x mean height) summed over a surface object, or None."""
    if isinstance(geom, r3.Extrusion):
        geom = geom.ToBrep(False)
    if isinstance(geom, r3.Brep):
        if geom.IsSolid:
            return None  # a closed solid isn't a roof or a base; measure it on its own layer
        total = [0.0, 0.0]
        for fi in range(len(geom.Faces)):
            r = face_prism(geom, fi)
            if r is None:
                return None
            total[0] += r[0]; total[1] += r[1]
        return tuple(total)
    if isinstance(geom, r3.Mesh):
        a, b = mesh_prism(geom)
        return (a, b) if a >= 0 else (-a, -b)
    return None


def mesh_triangles(mesh):
    v = mesh.Vertices
    out = []
    for i in range(len(mesh.Faces)):
        f = mesh.Faces[i]
        out.append([xyz(v[f[0]]), xyz(v[f[1]]), xyz(v[f[2]])])
        if f[2] != f[3]:
            out.append([xyz(v[f[0]]), xyz(v[f[2]]), xyz(v[f[3]])])
    return out


def surface_triangles(geom):
    """Triangles covering a surface object, from its render mesh, or None."""
    if isinstance(geom, r3.Extrusion):
        geom = geom.ToBrep(False)
    if isinstance(geom, r3.Mesh):
        return mesh_triangles(geom)
    if not isinstance(geom, r3.Brep) or geom.IsSolid:
        return None
    out = []
    for fi in range(len(geom.Faces)):
        face = geom.Faces[fi]
        mesh = face.GetMesh(r3.MeshType.Any)
        if mesh is not None:
            out.extend(mesh_triangles(mesh))
            continue
        # No render mesh saved: a single straight-edged convex loop can be fanned.
        loops = list(face.Loops)
        if len(loops) != 1:
            return None
        pts = []
        for trim in loops[0].Trims:
            if trim.EdgeIndex < 0:
                continue
            edge_pts = curve_points(geom.Edges[trim.EdgeIndex])
            if edge_pts is None:
                return None
            if trim.IsReversed:
                edge_pts = edge_pts[::-1]
            pts.extend(edge_pts[:-1])
        signs = {prism([pts[i - 1], pts[i], pts[(i + 1) % len(pts)]])[0] > 0 for i in range(len(pts))}
        if len(pts) < 3 or len(signs) != 1:
            return None
        out.extend([pts[0], pts[i], pts[i + 1]] for i in range(1, len(pts) - 1))
    return out


def rasterize(tris, cell, origin, take_max):
    """Height per grid cell (cell centres) over a set of triangles.

    Returns (cell keys, heights): the highest triangle over each cell if take_max,
    else the lowest."""
    keys, zs = [], []
    for t in tris:
        (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t
        d = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if abs(d) < 1e-12:
            continue  # vertical or degenerate: covers no area seen from above
        ix = np.arange(np.ceil((min(x0, x1, x2) - origin[0]) / cell - 0.5),
                       np.floor((max(x0, x1, x2) - origin[0]) / cell - 0.5) + 1)
        iy = np.arange(np.ceil((min(y0, y1, y2) - origin[1]) / cell - 0.5),
                       np.floor((max(y0, y1, y2) - origin[1]) / cell - 0.5) + 1)
        if not len(ix) or not len(iy):
            continue
        gx, gy = np.meshgrid(ix, iy)
        px = origin[0] + (gx.ravel() + 0.5) * cell
        py = origin[1] + (gy.ravel() + 0.5) * cell
        a = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) / d
        b = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) / d
        c = 1 - a - b
        inside = (a >= 0) & (b >= 0) & (c > 0)  # half-open, so shared edges count once
        keys.append(gx.ravel()[inside].astype(np.int64) * (1 << 32) + gy.ravel()[inside].astype(np.int64))
        zs.append((a * z0 + b * z1 + c * z2)[inside])
    if not keys:
        return np.zeros(0, np.int64), np.zeros(0)
    keys, zs = np.concatenate(keys), np.concatenate(zs)
    order = np.lexsort((-zs if take_max else zs, keys))
    keys, zs = keys[order], zs[order]
    first = np.r_[True, keys[1:] != keys[:-1]]
    return keys[first], zs[first]


def envelope_volume(roof_tris, base_tris, cell):
    """Volume between the highest roof and the ground, on a grid of square cells.

    Overlapping roofs (a roof surface stacked over another) and duplicated
    footprints are each counted once. Returns (volume, roof cells without a base,
    base cells without a roof) with the last two as areas."""
    allpts = np.array([p for t in roof_tris + base_tris for p in t])
    origin = allpts[:, :2].min(0)
    rk, rz = rasterize(roof_tris, cell, origin, take_max=True)
    bk, bz = rasterize(base_tris, cell, origin, take_max=False)
    both, ri, bi = np.intersect1d(rk, bk, assume_unique=True, return_indices=True)
    h = rz[ri] - bz[bi]
    area = cell * cell
    return float(h.sum() * area), (len(rk) - len(both)) * area, (len(bk) - len(both)) * area


def object_volume(geom):
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
        # A face's render mesh follows the surface, not the solid, so flip reversed faces.
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
    centre = xyz(geom.GetBoundingBox().Center)
    return abs(sum(s * mesh_volume(m, centre) for m, s in meshes))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="path to the .3dm file")
    ap.add_argument("--csv", default="volumes_by_layer.csv", help="where to write the per-layer table")
    ap.add_argument("--top-level", action="store_true",
                    help="add up sub-layers into their top-level layer (e.g. all of Buildings::*)")
    ap.add_argument("--roof-layer", help="layer holding roof surfaces (roof / base mode, see above)")
    ap.add_argument("--base-layer", help="layer holding footprint / ground surfaces (roof / base mode)")
    ap.add_argument("--cell", type=float,
                    help="grid size for roof / base mode, in model units (default: 1 ft)")
    args = ap.parse_args()
    if bool(args.roof_layer) != bool(args.base_layer):
        ap.error("--roof-layer and --base-layer go together")

    print(f"Reading {args.model} (this can take a few minutes for a big file)...", flush=True)
    model = r3.File3dm.Read(args.model)
    if model is None:
        sys.exit("Could not read that file. Is it a .3dm file, and is the path right?")

    unit = str(model.Settings.ModelUnitSystem).split(".")[-1]
    to_m3 = TO_METRES.get(unit, 1.0) ** 3
    layers = {model.Layers[i].Index: model.Layers[i].FullPath for i in range(len(model.Layers))}
    missing = [l for l in (args.roof_layer, args.base_layer) if l and l not in layers.values()]
    if missing:
        sys.exit(f"No layer called {missing[0]!r}. Layers in this file:\n  " + "\n  ".join(sorted(layers.values())))
    # [projected area, area x height, surfaces used, surfaces skipped] for roof and base
    prisms = {args.roof_layer: [0.0, 0.0, 0, 0], args.base_layer: [0.0, 0.0, 0, 0]}
    triangles = {args.roof_layer: [], args.base_layer: []}

    vol = defaultdict(float)
    counted = defaultdict(int)
    skipped = defaultdict(int)
    n = len(model.Objects)
    for i, obj in enumerate(model.Objects):
        layer = layers.get(obj.Attributes.LayerIndex, "(unknown layer)")
        if args.roof_layer and layer in prisms:
            p = surfaces_prism(obj.Geometry)
            t = surface_triangles(obj.Geometry)
            acc = prisms[layer]
            if p is None or t is None:
                acc[3] += 1
            else:
                acc[0] += p[0]; acc[1] += p[1]; acc[2] += 1
                triangles[layer].extend(t)
        if args.top_level:
            layer = layer.split("::")[0]
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

    if args.roof_layer:
        roof, base = prisms[args.roof_layer], prisms[args.base_layer]
        cell = args.cell or 0.3048 / TO_METRES.get(unit, 1.0)  # 1 ft by default
        print(f"\nWorking out roof / base volume on a {cell:g} {unit.lower()} grid...", flush=True)
        ev, roof_only, base_only = envelope_volume(triangles[args.roof_layer], triangles[args.base_layer], cell)
        pv = roof[1] - base[1]
        a2 = to_m3 ** (2 / 3)
        print(f"ROOF / BASE VOLUME: {ev * to_m3:,.0f} m3  ({ev:,.0f} cubic {unit.lower()})")
        print(f"  roofs: {roof[2]:,} surfaces, {roof[0] * a2:,.0f} m2 seen from above ({roof[3]:,} skipped)")
        print(f"  bases: {base[2]:,} surfaces, {base[0] * a2:,.0f} m2 seen from above ({base[3]:,} skipped)")
        print(f"  Adding up every roof without removing overlaps would give {pv * to_m3:,.0f} m3.")
        print(f"  {roof_only * a2:,.0f} m2 of roof has no base under it and {base_only * a2:,.0f} m2 of base"
              f" has no roof over it; neither counts.")


if __name__ == "__main__":
    main()
