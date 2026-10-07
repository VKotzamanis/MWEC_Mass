"""Write tests/baseline/matlab_v1_reference.json from the MATLAB v1.0 results Output/C1_*_results.mat.

The keys match tools/baseline_run.m (summarise), so the Octave baseline and the MATLAB reference
can be compared field by field. Test-only; requires numpy and scipy.

Usage: python3 tools/make_matlab_v1_reference.py [repo_root]
"""
import json
import os
import sys

import numpy as np
from scipy.io import loadmat

STAGE3_NAMES = ["t_wall", "vertical_shift", "draft", "t_max", "t_min_active", "V_wall",
                "V_air", "M_total", "CG_z_world", "GM_realised", "T_heave_realised",
                "T_pitch_realised", "mass_balance_error_pct", "phi_star", "feasible", "exitflag"]


def row(v):
    return [float(x) for x in np.atleast_1d(np.asarray(v, dtype=float)).ravel()]


def summarise(mode, results, final_props):
    s2 = results.stage2_3d
    c = results.steel_data if mode == "thin_shell" else results.constructability
    summary = {
        "mode": mode,
        "stage1_x": row(results.stage1_2d.x_optimal),
        "stage2_x": row(s2.x_optimal),
        "stage2_exitflag": float(s2.exitflag),
        "stage2_fval": float(s2.fval),
        "stage2_iterations": float(s2.output.iterations),
        "mass_total": float(final_props.mass_total),
        "mass_buoyant_force": float(final_props.mass_buoyant_force),
        "CG_total_z": float(np.atleast_1d(final_props.CG_total)[2]),
        "GM_L": float(final_props.GM_L),
        "T_heave_coupled": float(np.atleast_1d(final_props.coupled_periods)[1]),
        "T_pitch_coupled": float(np.atleast_1d(final_props.coupled_periods)[2]),
        "vertical_shift": float(final_props.vertical_shift),
        "draft": float(final_props.draft),
        "realised_strip_density": row(final_props.realised_strip_density),
    }
    stage3 = {}
    # The wall material names the thickness and volume keys; the v1.0 field names are t_steel and
    # V_steel in both modes (they held UHPC quantities in modular precast).
    material = "steel" if mode == "thin_shell" else "uhpc"
    for name in STAGE3_NAMES:
        key = name.replace("wall", material)
        for source in (key, key.replace("uhpc", "steel")):
            if hasattr(c, source):
                stage3[key] = float(np.asarray(getattr(c, source), dtype=float).ravel()[0])
                break
    for name in ("z_ballast", "z_fill"):
        if hasattr(c, name):
            stage3["z_ballast"] = float(np.asarray(getattr(c, name), dtype=float).ravel()[0])
            break
    if hasattr(c, "t_offset_strip") and np.asarray(c.t_offset_strip).size:
        t = np.asarray(c.t_offset_strip, dtype=float).ravel()
        stage3["t_offset_strip_finite"] = [float(x) for x in t[np.isfinite(t)]]
    summary["stage3"] = stage3
    # Stage-2 design (results.Final3D), the reference that Stage 3 is compared with.
    f = results.Final3D
    summary["stage2"] = {
        "vertical_shift": float(f.vertical_shift),
        "draft": float(f.draft),
        "mass_total": float(f.mass_total),
        "CG_total_z": float(np.atleast_1d(f.CG_total)[2]),
        "GM_L": float(f.GM_L),
        "T_heave_coupled": float(np.atleast_1d(f.coupled_periods)[1]),
        "T_pitch_coupled": float(np.atleast_1d(f.coupled_periods)[2]),
    }
    return summary


def main(argv):
    root = argv[1] if len(argv) > 1 else "."
    out = {"source": "MATLAB v1.0 results, Output/C1_<mode>_results.mat"}
    for mode in ("modular_precast", "thin_shell"):
        d = loadmat(os.path.join(root, "Output", "C1_%s_results.mat" % mode),
                    squeeze_me=True, struct_as_record=False)
        out[mode] = summarise(mode, d["results"], d["final_props"])
    path = os.path.join(root, "tests", "baseline", "matlab_v1_reference.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")


if __name__ == "__main__":
    main(sys.argv)
