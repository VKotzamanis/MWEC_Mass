"""Write small STEP files with gmsh/OpenCASCADE for tests/step/test_step_check.m (test-only).

Usage: python3 -I tests/step/make_step_fixtures.py OUTDIR
Writes cube.step (unit cube, declared METRE), cylinder.step (radius 0.5, height 2), void.step
(4 x 4 x 4 box with a sphere of radius 1 cut out, an inner shell), sheet.step (unit square surface,
no solid) and cube_mm.step (the cube with its unit left at gmsh's default MILLIMETRE label).
gmsh writes coordinates as they are and labels them MILLI.METRE, so the label is replaced by METRE.
"""
import os
import sys

import gmsh


def write(model_builder, path):
    gmsh.model.add(os.path.basename(path))
    model_builder()
    gmsh.model.occ.synchronize()
    gmsh.write(path)
    gmsh.model.remove()


def relabel_metre(path):
    text = open(path).read()
    assert "SI_UNIT(.MILLI.,.METRE.)" in text
    open(path, "w").write(text.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)"))


def main(out):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ

    def void():
        box = occ.addBox(0, 0, 0, 4, 4, 4)
        sphere = occ.addSphere(2, 2, 2, 1)
        occ.cut([(3, box)], [(3, sphere)])

    jobs = {"cube.step": lambda: occ.addBox(0, 0, 0, 1, 1, 1),
            "cylinder.step": lambda: occ.addCylinder(0, 0, 0, 0, 0, 2, 0.5),
            "void.step": void,
            "sheet.step": lambda: occ.addRectangle(0, 0, 0, 1, 1),
            "cube_mm.step": lambda: occ.addBox(0, 0, 0, 1, 1, 1)}
    for name, build in jobs.items():
        path = os.path.join(out, name)
        write(build, path)
        if name != "cube_mm.step":
            relabel_metre(path)
    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1])
