# The HAMS-MREL route

This suite can build its hydrodynamic-coefficient cache in two ways: from a WAMIT run set or directly from the hull geometry with the HAMS-MREL solver. This page documents the second route.

At each draft (vertical shift of the hull relative to the waterline) the route meshes the wetted
hull and the waterplane lid, writes the HAMS-MREL solver's four input files, runs the solver, and
parses the added mass, radiation damping and diffraction excitation it returns into the same
25-field `hydro_table` structure the WAMIT route also produces. The rest of the suite reads either
cache in the same way, so downstream code does not need to know which route built it.

## What the user supplies

The route calls a solver binary that is not distributed with this repository. Obtain or build the binary from the [HAMS-MREL GitHub repository](https://github.com/vaibhavraghavan/HAMS-MREL), an open-source multi-body boundary-element solver for marine renewable energies developed by Raghavan, Loukogeorgaki, Mantadakis, Metrikine, and Lavidas. Consult the upstream repository for its current license information. The solver is cited as follows.

- Raghavan, V., Loukogeorgaki, E., Mantadakis, N., Metrikine, A. V., and Lavidas, G. (2024). _HAMS-MREL, a new open-source multiple body solver for marine renewable energies: Model description, application and validation._ _Renewable Energy_, 237, 127857. https://doi.org/10.1016/j.renene.2024.121577
- Raghavan, V., Metrikine, A. V., and Lavidas, G. (2026). “Theory and validation of the new features in BIEM solver HAMS-MREL.” _Journal of Ocean Engineering and Marine Energy_, 12, 53-71. https://doi.org/10.1007/s40722-025-00431-8

The binary must be built or obtained by the user and placed at `Input/HAMS_MREL/HAMS_MREL`. This route currently supports Linux only.

Before running the binary, the route checks that the Intel oneAPI runtime is reachable, since the
binary links against it dynamically and Intel's own libraries are not normally on a distribution's
library path. The following is a shortened MATLAB excerpt of that guard (the ellipsis is
intentional; it stands for the full error-construction expression):

```matlab
if ~ispc
    ld_path = getenv('LD_LIBRARY_PATH');
    if ~contains(ld_path, 'intel/oneapi') && ~contains(ld_path, 'intel\oneapi')
        error(... 'Intel oneAPI runtime not on LD_LIBRARY_PATH. HAMS-MREL is linked against ' ...
              'Intel MKL + iomp5 and cannot load them without setvars.sh.' ...)
```

The check runs before the solver is ever invoked, so a missing runtime stops with this message
rather than the binary's own cryptic "shared object file not found" abort. The fix is to source
Intel's `setvars.sh` in the shell MATLAB is launched from, then start MATLAB from that same shell.

## Switching it on

Four settings in `WEC_User_Input.m` control the route:

- `in.bem.run_HAMS_MREL`: `true` runs this route; `false` (the default) loads the WAMIT cache instead and does not call the solver.
- `in.bem.hams_dir`: the solver's run workspace, default `fullfile('Output', 'hams_mrel')`.
- `in.bem.hams_exe`: the solver binary path, default `fullfile('Input', 'HAMS_MREL', 'HAMS_MREL')`.
- `in.bem.hams_cache_file`: the name of the cache this route writes, default
  `'C1_hams_cache.mat'`, resolved under `in.bem.hams_dir`.

## What it writes

`in.bem.hams_dir` is created at run time, not pre-populated, with:

- `Input/`: the four solver input files this route writes at each draft: hull and waterplane meshes, the control file, and the hydrostatic-matrix file.
- `Output/Hams_format/`, `Output/Hydrostar_format/`, and `Output/Wamit_format/`: the solver's output folders. The route reads the WAMIT-format `.1` files for added mass and radiation damping, and `.3` files for excitation, from `Output/Wamit_format/`.
- `Output/ErrorCheck.txt`: created empty because the solver requires the file to exist before it runs.

The finished `hydro_table` (coarse draft nodes plus adaptively refined ones) is saved under
`in.bem.hams_dir` as `in.bem.hams_cache_file`, by default
`Output/hams_mrel/C1_hams_cache.mat`. This is deliberately a different file from the one the WAMIT
route reads, `Input/C1_wamit_cache.mat`, so switching this route on cannot overwrite the WAMIT cache. To use a cache produced by this route in the main pipeline, point `in.bem.bem_cache_file` to a copy placed under `Input/`.

## Function map

`+mwecmass.bem.hams_mrel`:

- `run`: builds a `hydro_table` by running HAMS-MREL. It sets the mesh size, builds the geometry configuration, and generates the draft nodes.
- `generate_draft_nodes`: generates a coarse and adaptively refined list of draft (vertical-shift) nodes for HAMS runs.
- `run_at_draft`: runs the full HAMS pipeline at one vertical shift: mesh generation, input writing, solver execution, and parsing.
- `run_solver`: executes the HAMS-MREL solver binary and scans its output for abort markers that the binary does not always report through its exit status.
- `HamsWriter`: a static-method class that writes the four fixed-format solver input files: `control_file`, `hydrostatic_file`, `pnl_file`, and `matrix_6x6`.
- `parse_wamit_1_file`: reads a HAMS-produced WAMIT `.1`-format output file containing dimensionalised added mass and radiation damping.
- `parse_wamit_3_file`: reads a HAMS-produced WAMIT `.3`-format output file containing dimensionalised excitation force.
- `band_averaged_damping`: averages radiation damping over a period band or at the nearest frequency.
- `default_hams_params`: returns default HAMS analysis parameters, including the period grid, for a geometry configuration.
- `hams_constants`: HAMS-internal normalisation constants for water density and gravity, used only on this route.

Mesh functions this route calls, from `+mwecmass.mesh`:

- `panel_grid_counts`: selects the panel grid (`Nu`, `Nv`) from a target panel size; used only on this route.
- `generate`: builds a panel mesh with vertices, panels, and normals from the parsed hull model at a draft and grid density.
- `hull_waterline_polygon`: extracts an ordered waterline polygon from a hull mesh.
- `waterplane_mesh_unstructured`: meshes the waterplane lid by constrained Delaunay triangulation followed by greedy quadrilateral matching.

A second waterplane mesher, `waterplane_mesh_structured` (two-curve transfinite interpolation), is used only by the panel-mesh diagnostic figure and is not called by this route.

## Release status

The HAMS-MREL binary and Intel oneAPI runtime are not distributed with this release. Users who need this optional route must obtain or build HAMS-MREL separately from the [HAMS-MREL GitHub repository](https://github.com/vaibhavraghavan/HAMS-MREL) and place the executable at the documented path. The repository includes the route implementation, input writers, output parsers, and configuration required to run it.

`in.bem.run_HAMS_MREL` is `false` by default, and the supplied hydrodynamic cache for this hull was computed with WAMIT at a 74 m water depth. HAMS-MREL runs were performed before release; the binary and its runtime dependencies are simply not included in this distribution.
