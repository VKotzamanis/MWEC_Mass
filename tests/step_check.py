"""Import a STEP file with gmsh (OpenCASCADE) and print JSON with checks on the imported model.

usage: python3 step_check.py FILE.step [--mesh-size H]

The model is imported with Geometry.OCCTargetUnit = M, so coordinates and volumes are in metres
and cubic metres whatever length unit the file declares; a file declaring MILLI therefore comes
back 1000 times smaller.

JSON keys
  declared_length_unit  prefix and unit of the file's LENGTH_UNIT, read from the file text
  n_volumes, n_surfaces entity counts after import
  bbox_occ              [xmin, ymin, zmin, xmax, ymax, zmax] from OpenCASCADE (includes its shape tolerance)
  bbox_mesh             the same box over the surface-mesh nodes
  occ_volumes           OpenCASCADE getMass of every volume
  mesh_volumes          divergence-theorem volume of the surface mesh of every volume; gmsh
                        orients each triangle by the orientation of its face in the imported shell
                        (outward on the outer shell, into the cavity on a void), so no flip is applied
  open_edges            mesh edges whose incident triangles do not pair up as one forward and one
                        reverse use, summed over the volumes and the surfaces that bound no volume
  mesh_size             element size used for the surface mesh (default: bbox diagonal / 20)
"""
import json
import re
import sys

import gmsh
import numpy as np


def declared_length_unit(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = re.search(r"LENGTH_UNIT\(\)[^;]*?SI_UNIT\(\s*(\$|\.[A-Z]+\.)\s*,\s*\.([A-Z]+)\.\s*\)", text)
    if not m:
        return None
    prefix = "" if m.group(1) == "$" else m.group(1).strip(".")
    return prefix + m.group(2)


def open_edge_count(tri):
    if len(tri) == 0:
        return 0
    a = np.concatenate([tri[:, 0], tri[:, 1], tri[:, 2]])
    b = np.concatenate([tri[:, 1], tri[:, 2], tri[:, 0]])
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    sign = np.where(a < b, 1, -1)
    keys, inv = np.unique(np.stack([lo, hi], axis=1), axis=0, return_inverse=True)
    inv = inv.ravel()
    count = np.bincount(inv, minlength=len(keys))
    net = np.bincount(inv, weights=sign, minlength=len(keys))
    return int(np.count_nonzero((count != 2) | (net != 0)))


def surface_triangles(tag, node_index):
    types, _, nodes = gmsh.model.mesh.getElements(2, tag)
    out = [np.array([node_index[n] for n in conn], dtype=np.int64).reshape(-1, 3)
           for t, conn in zip(types, nodes) if t == 2]
    return np.vstack(out) if out else np.zeros((0, 3), dtype=np.int64)


def main(argv):
    path = argv[1]
    mesh_size = None
    if "--mesh-size" in argv:
        mesh_size = float(argv[argv.index("--mesh-size") + 1])

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setString("Geometry.OCCTargetUnit", "M")
    gmsh.model.add("step_check")
    gmsh.model.occ.importShapes(path)
    gmsh.model.occ.synchronize()

    volumes = [t for _, t in gmsh.model.getEntities(3)]
    surfaces = [t for _, t in gmsh.model.getEntities(2)]
    bbox = list(gmsh.model.getBoundingBox(-1, -1))
    diag = float(np.linalg.norm(np.array(bbox[3:]) - np.array(bbox[:3])))
    if mesh_size is None:
        mesh_size = diag / 20.0
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.model.mesh.generate(2)

    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.asarray(coords).reshape(-1, 3)
    node_index = {int(t): i for i, t in enumerate(ntags)}

    occ_volumes = [gmsh.model.occ.getMass(3, v) for v in volumes]
    groups = []
    mesh_volumes = []
    bounding = set()
    for v in volumes:
        tris = []
        for _, s in gmsh.model.getBoundary([(3, v)], combined=False, oriented=True):
            tris.append(surface_triangles(abs(s), node_index))
            bounding.add(abs(s))
        tri = np.vstack(tris) if tris else np.zeros((0, 3), dtype=np.int64)
        groups.append(tri)
        p = xyz - np.array(bbox[:3])
        a, b, c = p[tri[:, 0]], p[tri[:, 1]], p[tri[:, 2]]
        mesh_volumes.append(float(np.einsum("ij,ij->", a, np.cross(b, c)) / 6.0))
    free = [surface_triangles(s, node_index) for s in surfaces if s not in bounding]
    if free:
        groups.append(np.vstack(free))

    all_tri = np.vstack([surface_triangles(s, node_index) for s in surfaces]) if surfaces else np.zeros((0, 3), dtype=np.int64)
    used = xyz[np.unique(all_tri)] if len(all_tri) else np.zeros((0, 3))
    bbox_mesh = list(used.min(axis=0)) + list(used.max(axis=0)) if len(used) else None

    result = {
        "file": path,
        "gmsh_version": gmsh.__version__,
        "declared_length_unit": declared_length_unit(path),
        "n_volumes": len(volumes),
        "n_surfaces": len(surfaces),
        "bbox_occ": bbox,
        "bbox_mesh": bbox_mesh,
        "occ_volumes": occ_volumes,
        "mesh_volumes": mesh_volumes,
        "open_edges": sum(open_edge_count(g) for g in groups),
        "mesh_size": mesh_size,
    }
    gmsh.finalize()
    print(json.dumps(result))


if __name__ == "__main__":
    main(sys.argv)
