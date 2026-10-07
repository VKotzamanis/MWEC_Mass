"""Verify a STEP file with gmsh/OpenCASCADE (test-only).

Reports, as JSON on stdout: the declared length unit (read from the file text), the number of
solids and faces, the number of open boundary edges, per solid the exact OpenCASCADE volume, the
volume from the divergence theorem on a surface mesh, and the bounding box. The file is imported
with Geometry.OCCTargetUnit = M, so lengths are metres whatever unit the file declares; the declared
unit is reported separately.

An open boundary edge is a non-degenerate edge used by exactly one face. A closed solid has none:
the shared edges of two faces are used twice, and the seam of a periodic face is listed twice in
that face's boundary.

Usage:
  python3 -I tests/step_check.py FILE.step [--mesh-size H] [--expect-solids N] [--expect-units METRE]
         [--expect-open-edges N] [--expect-volume V] [--volume-rtol R] [--bbox-atol A]
Exit status 1 if any given --expect-* check fails; without them the file is only reported.
--expect-volume V is compared with the sum of the exact solid volumes (relative tolerance R, default
1e-9). --mesh-size H is the maximum surface-mesh edge length (default 1/80 of the bounding-box
diagonal); the mesh volume carries the chord error of the flat triangles and is reported, not gated.
"""
import argparse
import json
import re
import sys
from collections import Counter

import gmsh
import numpy as np


def declared_length_units(path):
    """Length units named in the file, for example ['METRE'] or ['MILLIMETRE']."""
    text = open(path, errors="replace").read()
    units = []
    for m in re.finditer(r"LENGTH_UNIT\s*\(\s*\)\s*NAMED_UNIT\s*\(\s*\*\s*\)\s*SI_UNIT\s*\(\s*([^,]*?)\s*,\s*\.(\w+)\.\s*\)", text):
        prefix = m.group(1).strip().strip(".")
        prefix = "" if prefix in ("$", "") else prefix
        units.append(prefix + m.group(2))
    for m in re.finditer(r"CONVERSION_BASED_UNIT\s*\(\s*'([^']*)'", text):
        units.append(m.group(1).upper())
    return units


def surface_mesh_volume(solid_tag, h):
    """Divergence-theorem volume of the triangulated boundary of one solid.

    The triangles of a face are consistently oriented but not necessarily outward, so the sign of
    each face comes from a probe: the point of the face nearest to one triangle, moved a small step
    along that triangle's normal, is classified against the solid by OpenCASCADE (inside: the
    triangles point inward).
    """
    boundary = gmsh.model.getBoundary([(3, solid_tag)], combined=False, oriented=False)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    index = {int(t): i for i, t in enumerate(node_tags)}
    total = 0.0
    for _, face in boundary:
        _, _, nodes = gmsh.model.mesh.getElements(2, abs(face))
        face_sum = 0.0
        sign = 0.0
        for block in nodes:
            tri = np.array([[index[int(n)] for n in row] for row in block.reshape(-1, 3)])
            p0, p1, p2 = xyz[tri[:, 0]], xyz[tri[:, 1]], xyz[tri[:, 2]]
            face_sum += np.sum(np.einsum("ij,ij->i", p0, np.cross(p1, p2))) / 6.0
            if sign == 0.0:
                k = len(p0) // 2
                n = np.cross(p1[k] - p0[k], p2[k] - p0[k])
                n = n / np.linalg.norm(n)
                centroid = (p0[k] + p1[k] + p2[k]) / 3.0
                on_face = np.array(gmsh.model.getClosestPoint(2, abs(face), list(centroid))[0])
                sign = -1.0 if gmsh.model.isInside(3, solid_tag, list(on_face + 1e-3 * h * n)) else 1.0
        total += sign * face_sum
    return float(total)


def check(path, mesh_size=None):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setString("Geometry.OCCTargetUnit", "M")
    try:
        gmsh.model.add("step_check")
        gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        solids = [t for _, t in gmsh.model.getEntities(3)]
        faces = [t for _, t in gmsh.model.getEntities(2)]
        uses = Counter()
        for f in faces:
            for _, e in gmsh.model.getBoundary([(2, f)], combined=False, oriented=False):
                uses[abs(e)] += 1
        # Degenerate edges (the poles of a sphere) have zero length and are used once by design.
        open_edges = sorted(e for e, n in uses.items()
                            if n == 1 and gmsh.model.occ.getMass(1, e) > 1e-12)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
        h = mesh_size or float(np.linalg.norm([xmax - xmin, ymax - ymin, zmax - zmin])) / 80.0
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.mesh.generate(2)
        per_solid = []
        for t in solids:
            bb = gmsh.model.getBoundingBox(3, t)
            per_solid.append({"tag": t,
                              "volume_occ": gmsh.model.occ.getMass(3, t),
                              "volume_mesh": surface_mesh_volume(t, h),
                              "bbox": list(bb)})
        return {"file": path,
                "declared_length_units": declared_length_units(path),
                "n_solids": len(solids),
                "n_faces": len(faces),
                "n_open_edges": len(open_edges),
                "open_edges": open_edges,
                "mesh_size": h,
                "solids": per_solid,
                "volume_occ_total": float(sum(s["volume_occ"] for s in per_solid)),
                "bbox": [xmin, ymin, zmin, xmax, ymax, zmax]}
    finally:
        gmsh.finalize()


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--mesh-size", type=float)
    ap.add_argument("--expect-solids", type=int)
    ap.add_argument("--expect-units")
    ap.add_argument("--expect-open-edges", type=int)
    ap.add_argument("--expect-volume", type=float)
    ap.add_argument("--volume-rtol", type=float, default=1e-9)
    args = ap.parse_args(argv[1:])
    report = check(args.file, args.mesh_size)
    failures = []
    if args.expect_solids is not None and report["n_solids"] != args.expect_solids:
        failures.append("solids %d, expected %d" % (report["n_solids"], args.expect_solids))
    if args.expect_units is not None and report["declared_length_units"] != [args.expect_units]:
        failures.append("declared units %s, expected [%s]" % (report["declared_length_units"], args.expect_units))
    if args.expect_open_edges is not None and report["n_open_edges"] != args.expect_open_edges:
        failures.append("open edges %d, expected %d" % (report["n_open_edges"], args.expect_open_edges))
    if args.expect_volume is not None:
        err = abs(report["volume_occ_total"] - args.expect_volume) / abs(args.expect_volume)
        report["volume_relative_error"] = err
        if err > args.volume_rtol:
            failures.append("volume %.17g, expected %.17g (relative error %.3e > %.1e)"
                            % (report["volume_occ_total"], args.expect_volume, err, args.volume_rtol))
    report["failures"] = failures
    print(json.dumps(report, indent=1))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
