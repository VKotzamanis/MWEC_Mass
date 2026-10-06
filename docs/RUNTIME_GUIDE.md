# Runtime guide

This guide is the operational reference for the WEC mass-distribution suite. It explains how to
install and run the code, defines every author-editable input and output option, and specifies the
geometry and hydrodynamic-cache contracts. For the governing equations and modelling assumptions,
see [Methods and computational engine](METHODS_ENGINE.md); for exported MAT-file fields, see
[RESULT_SCHEMA.md](RESULT_SCHEMA.md).

## Contents

- [1. Purpose and prerequisites](#1-purpose-and-prerequisites)
- [2. Choose and install an archive](#2-choose-and-install-an-archive)
- [3. First run](#3-first-run)
- [4. Execution lifecycle](#4-execution-lifecycle)
- [5. Coordinates, signs, units, and degrees of freedom](#5-coordinates-signs-units-and-degrees-of-freedom)
- [6. `WEC_User_Input.m` reference](#6-wec_user_inputm-reference)
- [7. `WEC_Output_Options.m` reference](#7-wec_output_optionsm-reference)
- [8. Required hydrodynamic-cache MAT schema](#8-required-hydrodynamic-cache-mat-schema)
- [9. MultiSurf `.ms2` geometry input](#9-multisurf-ms2-geometry-input)
- [10. Output files and directories](#10-output-files-and-directories)
- [11. Reproducible runs](#11-reproducible-runs)
- [12. Troubleshooting](#12-troubleshooting)
- [13. Safe customization examples](#13-safe-customization-examples)
- [14. Model scope and operating conditions](#14-model-scope-and-operating-conditions)

## 1. Purpose and prerequisites

The suite optimises the vertical distribution of mass in a free-floating wave-energy-converter
hull, then optionally realises that distribution as a thin steel shell with fill or as modular
precast UHPC with internal voids. One run evaluates one realisation type.

Requirements:

- MATLAB R2022b or later.
- Optimization Toolbox (`fmincon` and `optimoptions`).
- Either `WEC_MassOptimization.zip` or `Example_C1.zip`.
- Write permission in the repository, because runs create or replace files under `Output/`.

The complete C1 archive uses `Input/C1_wamit_cache.mat`; it does **not** require WAMIT or
HAMS-MREL at runtime. An external BEM solver is required only when producing coefficients for a
new hull. The HAMS-MREL executable is **not bundled** with this repository. See
[HAMS_MREL_ROUTE.md](HAMS_MREL_ROUTE.md) before enabling that route.

> **Important:** run MATLAB from the extracted source root. The root-level entry points resolve
> `Input/`, `src/`, and `Output/` relative to the repository.

## 2. Choose and install an archive

The release provides two independent ZIP files. Choose one; they are not combined.

1. Extract `WEC_MassOptimization.zip` for the code and public documentation, or extract
   `Example_C1.zip` for an immediately runnable C1 copy.
2. Start MATLAB with the extracted archive root as the current folder.
3. If using the source-only archive, supply a compatible `.ms2` deck and hydrodynamic cache under
   `Input/` before running. The complete archive already supplies `Input/C1.ms2` and
   `Input/C1_wamit_cache.mat`.

The installed tree relevant to a runtime user is:

```text
WEC_MassOptimization/                      # source-only choice
├── README.md
├── CHANGELOG.md
├── WEC_User_Input.m
├── WEC_Output_Options.m
├── Input/
│   ├── HAMS_MREL/                     # location for a user-supplied executable
│   └── WAMIT/                         # public cache-preparation helper code
├── docs/
│   ├── RUNTIME_GUIDE.md
│   ├── METHODS_ENGINE.md
│   ├── RESULT_SCHEMA.md
│   └── HAMS_MREL_ROUTE.md
├── src/
│   └── +mwecmass/
└── validation/
    └── diagnostics/                   # optional post-processing helpers

Example_C1/                                # complete choice; same tree plus:
├── Input/
│   ├── C1.ms2
│   └── C1_wamit_cache.mat
└── Output/                              # 31 verified C1 artifacts
```

The complete archive includes reference outputs; a new run with the default repeatable filenames
replaces files having the same names. The source-only archive creates `Output/` as needed.

## 3. First run

Start MATLAB with the extracted source root as the current folder. The simplest run is:

```matlab
WEC_User_Input
```

With no output argument, `WEC_User_Input` changes to the repository root, adds `src/` to the
MATLAB path, loads `WEC_Output_Options`, and runs the pipeline. The shipped default realisation is
`'thin_shell'`.

To inspect inputs without running:

```matlab
in = WEC_User_Input();
```

To run one mode programmatically:

```matlab
addpath(fullfile(pwd,'src'));
in = WEC_User_Input();
out = WEC_Output_Options();
in.materials.realisation_type = 'modular_precast';
mwecmass.driver.run(in, out);
```

To run all modes:

```matlab
addpath(fullfile(pwd,'src'));
in = WEC_User_Input();
out = WEC_Output_Options();

for mode = {'preliminary', 'thin_shell', 'modular_precast'}
    in.materials.realisation_type = mode{1};
    mwecmass.driver.run(in, out);
end
```

The three modes are:

| Value | Calculation after Stage 2 | Principal mode-specific record |
|---|---|---|
| `'preliminary'` | No material-realisation solve | Optimiser result passes through |
| `'thin_shell'` | Steel shell and solid-fill solve | `results.steel_data` |
| `'modular_precast'` | UHPC wall/ring and void solve | `results.constructability` |

## 4. Execution lifecycle

```text
WEC_User_Input + WEC_Output_Options
                 │
                 ▼
resolve Input/ and Output/ paths
                 │
                 ▼
parse MultiSurf hull and load/validate hydrodynamic cache
                 │
                 ▼
build hydrostatic, strip-geometry, and runtime tables
                 │
                 ▼
Stage 1: 2-D search / screening / calibration
                 │
                 ▼
Stage 2: full 3-D constrained optimisation
                 │
                 ▼
selected material realisation
                 │
                 ▼
figures + logs + results MAT
                 │
                 ▼
optional post-processing diagnostics
```

The default BEM branch loads and validates the WAMIT-format cache. When
`in.bem.run_HAMS_MREL = true`, the alternative branch creates/updates a cache in the configured
HAMS workspace before optimisation. Stage 1 generates a warm start; Stage 2 solves the full 3-D
problem; the selected realiser converts equivalent strip densities to physical material geometry.
Output selection occurs after the calculation and does not alter the optimised values.

## 5. Coordinates, signs, units, and degrees of freedom

### 5.1 Units

Use SI units throughout unless a field explicitly says otherwise:

- length: m; area: m²; volume: m³;
- mass: kg; density: kg/m³; mass moment of inertia: kg·m²;
- time: s; circular frequency: rad/s;
- force: N; translational stiffness: N/m; rotational stiffness: N·m/rad;
- plot font and line sizes: pt; figure sizes: cm; raster resolution: dpi;
- colours: RGB rows with each component in `[0,1]`.

### 5.2 Frames and signs

Hull geometry is expressed in a body-fixed Cartesian frame. The world frame places the still-water
surface at `z = 0`; `z` is positive upward. The axes correspond to surge `x`, sway `y`, and heave
`z`. Water depth is positive downward; `-1` selects deep water where supported.

`vertical_shift` is a signed translation, not draft:

```text
z_world = z_body + vertical_shift
z_waterline_in_body = -vertical_shift
draft = abs(hull_z_min + vertical_shift)
```

A **positive** `vertical_shift` moves the hull upward and normally makes it shallower. A negative
value moves it downward. `draft` is the positive keel depth below the world waterline. The cache
field named `drafts` stores `vertical_shift` values for historical compatibility; do not replace
them with positive physical draft values.

### 5.3 Degree-of-freedom order

All full hydrodynamic and rigid-body matrices use:

```text
[surge, sway, heave, roll, pitch, yaw] = [1, 2, 3, 4, 5, 6]
```

The optimiser uses the reduced indices `[1,3,5]`, in the order `[surge, heave, pitch]`. The model
is free-floating: surge hydrostatic stiffness is zero and `K_pto` is the zero matrix. See
[Methods and computational engine](METHODS_ENGINE.md) for the mass, added-mass, restoring, and
natural-period equations.

## 6. `WEC_User_Input.m` reference

Every leaf below exists in the shipped input function. **Default** means the shipped value;
alternatives are valid only when the stated constraints and dependencies remain satisfied.
The authoritative definitions are in `WEC_User_Input.m`; propagation and validation occur in
`src/+mwecmass/+driver/build_config.m`.

### 6.1 Physical constants and files

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.constants.rho_water` | `1025` | double scalar | kg/m³ | Positive water density used in hydrostatics and WAMIT re-dimensionalization. |
| `in.constants.g` | `9.80665` | double scalar | m/s² | Positive gravitational acceleration. |
| `in.constants.wamit_L` | `1.0` | double scalar | m | Positive WAMIT reference length used by preprocessing. |
| `in.files.ms2_file` | `'C1.ms2'` | char row | - | Filename relative to `Input/`; must identify a readable supported MultiSurf deck. |

### 6.2 Geometry and discretization

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.geometry.panel_size` | `0.1` | positive double scalar | m | Target mean BEM hull-panel edge length; chiefly used by HAMS meshing. |
| `in.geometry.wp_target_edge` | `0.1` | positive double scalar | m | Target waterplane-lid edge length; a loaded cache may restore its recorded value. |
| `in.geometry.aw_table_dz` | `0.05` | positive double scalar | m | Target vertical spacing of the waterplane-area table; smaller is finer and slower. |
| `in.geometry.num_ballast_sections` | `5` | positive integer scalar | count | Number `N` of vertical ballast strips; at least 1, and at least 3 for `modular_precast`. |
| `in.geometry.n_z_levels` | `200` | positive integer scalar | count | Vertical resolution of the 2-D surrogate span table. |
| `in.geometry.eff_w_floor` | `0.05` | positive double scalar | m | Lower bound on the 2-D surrogate effective width. |

### 6.3 Realisation type and materials

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.materials.realisation_type` | `'thin_shell'` | char row | - | Exactly `'preliminary'`, `'thin_shell'`, or `'modular_precast'`. |
| `in.materials.thin_shell.rho_shell` | `7500` | positive double scalar | kg/m³ | Structural shell density. |
| `in.materials.thin_shell.rho_void` | `1.2` | nonnegative double scalar | kg/m³ | Void/air density in the thin-shell model. |
| `in.materials.thin_shell.rho_fill` | `rho_shell` | double scalar | kg/m³ | Solid-fill density; must be `>= rho_shell` for the monotone analytic warm-start seed. |
| `in.materials.thin_shell.t_init` | `0.02` | positive double scalar | m | Initial thickness guess; it may be below `t_min` because it is a numerical seed. |
| `in.materials.thin_shell.t_min` | `0.025` | positive double scalar | m | Minimum shell thickness. |
| `in.materials.thin_shell.max_slope_factor` | `5.0` | positive double scalar | - | Limit controlling thickness taper with hull slope. |
| `in.materials.thin_shell.n_z_grid` | `300` | positive integer scalar | count | Vertical integration/solve grid. |
| `in.materials.modular_precast.rho_hull` | `2500` | positive double scalar | kg/m³ | UHPC density. |
| `in.materials.modular_precast.rho_fill` | `1.2` | nonnegative double scalar | kg/m³ | Internal void/air density; distinct from thin-shell solid fill. |
| `in.materials.modular_precast.t_min` | `0.0762` | positive double scalar | m | Minimum precast thickness (3 in by default). |
| `in.materials.modular_precast.wall_height` | `1.8` | positive double scalar | m | Height assigned to the wall region. |
| `in.materials.modular_precast.n_sub` | `100` | positive integer scalar | count | Vertical samples per strip. |
| `in.materials.modular_precast.t_init` | `thin_shell.t_init*(thin_shell.rho_shell/rho_hull)` | positive double scalar | m | Density-scaled initial thickness; recomputed from the referenced defaults. |
| `in.materials.modular_precast.max_slope_factor` | `thin_shell.max_slope_factor` | positive double scalar | - | Precast taper limit; mirrors the thin-shell default. |
| `in.materials.modular_precast.n_z_grid` | `thin_shell.n_z_grid` | positive integer scalar | count | Precast solve grid; mirrors the thin-shell default. |

### 6.4 Bounds and targets

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.bounds.ballast_density_bounds` | `[20,2500]` | double `[1x2]` | kg/m³ | Ordered `[lo,hi]` bound applied to every strip; constructability may tighten individual bounds. |
| `in.bounds.max_density_ratio` | `100.0` | double scalar | - | Must be `>=1`; limits each adjacent downward-to-upward ratio using the live guarded expression `rho_i/(rho_(i+1)+1 kg/m^3)`. See `METHODS_ENGINE.md` for the constraint definition. |
| `in.bounds.vertical_shift_bounds` | `[]` | empty or double `[1x2]` | m | `[]` auto-selects `[-hull_z_max+0.1, -hull_z_min-0.1]`; otherwise an ordered manual `[lo,hi]`. These are shifts, not drafts. |
| `in.targets.T_heave_goal` | `7.77` | positive double scalar | s | Heave natural-period design target. |
| `in.targets.T_pitch_goal` | `3.89` | positive double scalar | s | Pitch natural-period design target. |
| `in.targets.T_heave_range` | `[7.0,10.0]` | double `[1x2]` | s | Ordered target band used by the range penalty. |
| `in.targets.T_pitch_range` | `[3.0,5.0]` | double `[1x2]` | s | Ordered target band used by the range penalty. |
| `in.targets.gm_min` | `0.2` | double scalar | m | Minimum permitted longitudinal metacentric height. |
| `in.targets.gm_range` | `[0.2,0.7]` | double `[1x2]` | m | Nonnegative ordered reporting/objective band. |
| `in.targets.gm_target` | `0.5` | double scalar | m | Preferred `GM_L`; normally inside `gm_range` (outside it produces a warning). |

### 6.5 Stage-1 mode and PID tuning

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.pid.stage1_mode` | `'sweep'` | char row | - | `'sweep'`, `'skip'`, `'oneshot'`, or `'trained'`. Sweep ranks cached shifts; skip passes the initial point; oneshot solves once; trained iterates PID corrections. |
| `in.pid.k_vol_init` | `1.0` | double scalar | - | Initial volume-correction scale; unity means no correction. |
| `in.pid.k_gm_init` | `1.0` | double scalar | - | Initial GM-correction scale; unity means no correction. |
| `in.pid.mass_gains` | `[0.8,0.002,0.05]` | double `[1x3]` | - | Mass controller `[P,I,D]` gains; primarily relevant to trained Stage 1. |
| `in.pid.vol_gains` | `[0.6,0.01,0.03]` | double `[1x3]` | - | Volume controller `[P,I,D]` gains. |
| `in.pid.gm_gains` | `[0.3,0.02,0.02]` | double `[1x3]` | - | GM controller `[P,I,D]` gains. |
| `in.pid.mass_limits` | `[0.5,5.0]` | double `[1x2]` | - | Ordered saturation limits for the mass correction. |
| `in.pid.vol_limits` | `[-1.0,1.0]` | double `[1x2]` | - | Ordered saturation limits for the volume-controller output. |
| `in.pid.gm_limits` | `[-1.0,1.0]` | double `[1x2]` | - | Ordered saturation limits for the GM-controller output. |
| `in.pid.bounds_kvol` | `[0.30,3.00]` | double `[1x2]` | - | Ordered bounds on `k_vol`. |
| `in.pid.bounds_kgm` | `[0.50,3.00]` | double `[1x2]` | - | Ordered bounds on `k_gm`. |
| `in.pid.damping_vol_early` | `0.5` | double scalar | - | Volume update damping before the transition iteration. |
| `in.pid.damping_vol_late` | `0.7` | double scalar | - | Volume update damping at/after the transition iteration. |
| `in.pid.damping_gm_early` | `0.4` | double scalar | - | GM update damping before the transition iteration. |
| `in.pid.damping_gm_late` | `0.5` | double scalar | - | GM update damping at/after the transition iteration. |
| `in.pid.damping_transition_iter` | `2` | nonnegative integer scalar | count | Iteration at which late damping begins. |
| `in.pid.vol_conv_tol_pct` | `3.0` | nonnegative double scalar | % | Volume convergence tolerance. |
| `in.pid.gm_conv_tol` | `0.05` | nonnegative double scalar | m | GM convergence tolerance. |
| `in.pid.delta_kvol_stable` | `0.005` | nonnegative double scalar | - | Maximum `k_vol` change counted as stable. |
| `in.pid.delta_kgm_stable` | `0.01` | nonnegative double scalar | - | Maximum `k_gm` change counted as stable. |
| `in.pid.stable_count_needed` | `2` | positive integer scalar | count | Consecutive stable updates required. |
| `in.pid.mass_acceptable_pct` | `10` | nonnegative double scalar | % | Stage-1 acceptable mass-error band. |
| `in.pid.cg_guard_floor` | `0.01` | nonnegative double scalar | m | Avoids unstable ratios when either CG ordinate is near zero. |
| `in.pid.sat_proximity` | `0.01` | nonnegative double scalar | - | Distance from a correction bound treated as saturation. |
| `in.pid.tanh_shape_param` | `3.0` | positive double scalar | - | Shape of the initial heavy-at-keel density profile. |
| `in.pid.tanh_shape_fallback` | `1.5` | positive double scalar | - | Reduced shape used if the initial profile violates `max_density_ratio`. |

### 6.6 Objective, solver, and validation

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.objective.zone_k_amp` | `5.0` | nonnegative double scalar | - | Amplifies curvature of the range penalty outside the transition zones. |
| `in.objective.penalty_guard` | `1e4` | positive double scalar | - | Fallback objective returned by failed or nonfinite evaluations. |
| `in.solver.max_outer_iterations` | `50` | positive integer scalar | count | Maximum trained Stage-1 outer iterations. |
| `in.solver.n_sweep_refine` | `4` | nonnegative integer scalar | count | Stage-1 sweep refinement passes. |
| `in.solver.stage2_algorithm` | `'sqp'` | char row | - | Algorithm name accepted by `fmincon`; `'sqp'` is the supported production default. |
| `in.validation.autocad_Ixx` | `[]` | empty or double scalar | kg·m² | Optional external roll inertia. Applied only when the supplied `Iyy` discrepancy exceeds the threshold. |
| `in.validation.autocad_Iyy` | `[]` | empty or positive double scalar | kg·m² | Optional comparison trigger for computed pitch inertia; on excessive discrepancy, the supplied CAD values replace corresponding computed inertias. |
| `in.validation.autocad_Izz` | `[]` | empty or double scalar | kg·m² | Optional external yaw inertia, used with the `Iyy` discrepancy trigger. |
| `in.validation.autocad_discrepancy_pct` | `5` | nonnegative double scalar | % | Replacement threshold for the optional CAD inertia cross-check. |

### 6.7 BEM source and frequency grid

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.bem.run_HAMS_MREL` | `false` | logical scalar | - | `false`: load the WAMIT-format cache under `Input/`; `true`: run the user-installed HAMS-MREL route. |
| `in.bem.bem_cache_file` | `'C1_wamit_cache.mat'` | char row | - | Cache filename relative to `Input/`, used when `run_HAMS_MREL=false`. |
| `in.bem.hams_cache_file` | `'C1_hams_cache.mat'` | char row | - | Cache filename under `hams_dir`, used by the HAMS route. |
| `in.bem.water_depth` | `74` | double scalar | m | Positive downward; `-1` denotes deep water. Must match cache metadata within the runtime check. |
| `in.bem.hams_dir` | `fullfile('Output','hams_mrel')` | char path | - | Repository-relative HAMS workspace. |
| `in.bem.hams_exe` | `fullfile('Input','HAMS_MREL','HAMS_MREL')` | char path | - | Repository-relative path to a **user-supplied** executable. |
| `in.bem.T_min` | `3.0` | positive double scalar | s | Shortest BEM analysis period; HAMS/cache-generation setting. |
| `in.bem.T_max` | `20.0` | positive double scalar | s | Longest period; must exceed `T_min`. |
| `in.bem.T_step` | `0.5` | positive double scalar | s | Period-grid increment. |
| `in.bem.n_sweep_drafts` | `8` | positive integer scalar | count | Initial HAMS vertical-shift nodes. |
| `in.bem.n_adaptive_refine` | `5` | nonnegative integer scalar | count | Maximum added HAMS nodes in high-gradient intervals. |

### 6.8 Plot resolution and exported context

| Field | Default | Type / shape | Units | Allowed values, dependencies, and behavior |
|---|---:|---|---|---|
| `in.plots.viz_mesh_N` | `60` | positive integer scalar | count | Visualization mesh resolution per parametric direction. |
| `in.context.T_heave_goal` | `in.targets.T_heave_goal` | double scalar | s | Exported metadata mirror; keep synchronized with the live target. |
| `in.context.T_pitch_goal` | `in.targets.T_pitch_goal` | double scalar | s | Exported metadata mirror; keep synchronized with the live target. |
| `in.context.T_heave_range` | `in.targets.T_heave_range` | double `[1x2]` | s | Exported metadata mirror. |
| `in.context.T_pitch_range` | `in.targets.T_pitch_range` | double `[1x2]` | s | Exported metadata mirror. |
| `in.context.gm_target` | `in.targets.gm_target` | double scalar | m | Exported metadata mirror. |
| `in.context.T_surge_goal` | `8.95` | positive double scalar | s | Metadata only; does not constrain the optimiser. |
| `in.context.T_surge_range` | `[7.0,9.0]` | double `[1x2]` | s | Metadata only. |
| `in.context.WIS_station` | `''` | char row | - | Optional Wave Information Studies station identifier; metadata only. |
| `in.context.data_year` | `''` | char row | - | Optional data year or range; metadata only. |

## 7. `WEC_Output_Options.m` reference

Output flags select delivery and presentation only. Every `out.save` leaf must be a logical
scalar. A `true` selection can still be skipped when its mode/data gate is not satisfied.
The authoritative defaults are in `WEC_Output_Options.m`; runtime gating is implemented in
`src/+mwecmass/+output/dispatch.m` and the selected realisation package.

### 7.1 Save selections

| Field | Default | Output and dependency |
|---|---:|---|
| `out.save.stage1.density_2d` | `false` | Post-Stage-2 2-D equivalent-density projection, invoked after realisation with the optimised/final properties; the historical `stage1` option name does not make it a Stage-1 result. |
| `out.save.stage1.draft_landscape` | `true` | Draft-landscape figure; only available for `stage1_mode='skip'` with sweep data present. |
| `out.save.stage2.density_3d` | `false` | Optimiser equivalent-density 3-D figure. |
| `out.save.stage2.cross_section` | `true` | Final cross-section figure. |
| `out.save.stage2.summary_log` | `true` | `stage2_summary.log`. |
| `out.save.convergence` | `true` | Combined convergence figure; only when Stage 1 reports more than one iteration. |
| `out.save.stage3.steel_solve` | `true` | Steel-solve figure; thin-shell mode only. |
| `out.save.stage3.steel_solve_log` | `true` | `steel_solve.log`; thin-shell mode only. |
| `out.save.stage3.precast_midplane` | `true` | Precast XZ/midplane figure; modular-precast mode only. |
| `out.save.stage3.precast_strips` | `true` | Precast strip-plan figure; modular-precast mode only. |
| `out.save.stage3.realised_density_3d` | `false` | Realised-density 3-D figure when realised strip data exist. |
| `out.save.stage3.final_results_log` | `true` | `final_results.log`. |
| `out.save.stage3.diagnostic_panels_log` | `true` | `diagnostic_panels.log`. |
| `out.save.hydrodynamics.coefficients` | `true` | Added-mass/radiation-damping figure; current dispatcher also requires the configured HAMS executable and cache file to exist. |
| `out.save.hydrodynamics.raos` | `true` | RAO figure; same availability gate. |
| `out.save.hydrodynamics.added_mass_vs_draft` | `true` | Infinite-frequency added mass versus shift; same availability gate. |
| `out.save.hydrodynamics.mesh_diagnostic` | `false` | BEM mesh diagnostic; requires HAMS executable and readable geometry. |
| `out.save.hydrodynamics.panel_normals` | `true` | Panel-normal diagnostic; requires HAMS executable and readable geometry. |
| `out.save.hydrodynamics.cache_rewrite` | `true` | Permit trained Stage 1 to rewrite an enriched HAMS cache; it is not a figure. |
| `out.save.results_mat` | `true` | `Output/<hull>_<type>_results.mat`; must remain true for post-processing diagnostics. |
| `out.save.diagnostics.hull_at_draft` | `false` | **Opt-in.** Hull-at-draft verification image under `Output/diagnostics/`. |
| `out.save.diagnostics.uhpc_mass_balance` | `false` | **Opt-in.** UHPC mass-balance image/FIG; modular-precast only. |
| `out.save.diagnostics.stage_animations` | `false` | **Opt-in.** Stage animation GIFs; potentially slow and large. |

> **Diagnostics are post-processing, not optimisation requirements.** They run only after a
> results MAT file exists, and a diagnostic failure is converted to a warning so it does not
> invalidate the optimisation result.

### 7.2 File delivery

| Field | Default | Type / units | Allowed values and behavior |
|---|---:|---|---|
| `out.export.formats` | `{'png','fig','pdf'}` | cell array of char | Any subset/order of `'png'`, `'fig'`, and `'pdf'`; unknown extensions error. |
| `out.export.dpi` | `450` | positive scalar, dpi | Resolution used for PNG export. |
| `out.export.pdf_content` | `'vector'` | char | `exportgraphics` PDF content type, normally `'vector'`, `'image'`, or `'auto'`. |
| `out.export.timestamp_filenames` | `false` | logical scalar | If true, appends `_yyyyMMdd_HHmmss` to figure stems. |
| `out.console_echo` | `true` | logical scalar | Echo enabled reports to the MATLAB console as well as their log; false keeps the log. |
| `out.output_dir` | `'Output'` | repository-relative char path | Directory used by log writers. Standard driver figures and result MAT files use the repository `Output/` tree; retain this default for a single coherent run directory. |

### 7.3 Typography, layout, and strokes

| Field | Default | Type / units | Behavior |
|---|---:|---|---|
| `out.style.font_name` | `'Times New Roman'` | char | General figure font; the font must be installed for exact rendering. |
| `out.style.mono_font_name` | `'Courier New'` | char | Numeric-card font. |
| `out.style.font_size.title` | `14` | pt scalar | Figure/axes titles. |
| `out.style.font_size.axes` | `12` | pt scalar | Axis labels. |
| `out.style.font_size.tick_label` | `12` | pt scalar | Tick labels. |
| `out.style.font_size.colorbar` | `10.8` | pt scalar | Colorbar tick labels. |
| `out.style.font_size.legend` | `9.8` | pt scalar | Legend entries. |
| `out.style.font_size.annotation` | `11` | pt scalar | In-figure annotations. |
| `out.style.font_size.small_multiple` | `8` | pt scalar | Text inside small-multiple panels. |
| `out.style.colormap` | `'cividis'` | char | `'cividis'`, `'parula'`, or `'jet'`. |
| `out.style.line_width.main` | `2.0` | pt scalar | Legacy/main curve width. |
| `out.style.line_width.curve` | `2.0` | pt scalar | Primary data curves. |
| `out.style.line_width.boundary` | `1.4` | pt scalar | Hull/material boundaries. |
| `out.style.line_width.reference` | `1.2` | pt scalar | Targets and reference curves. |
| `out.style.line_width.axes` | `0.6` | pt scalar | Axes and ticks. |
| `out.style.line_width.grid` | `0.6` | pt scalar | Grid lines. |
| `out.style.line_width.hatch` | `0.4` | pt scalar | Void hatching. |
| `out.style.marker.size` | `6` | pt scalar | Marker size. |
| `out.style.figure_size.single_column` | `[8.5,6.65]` | cm `[width,height]` | Single-column figure class. |
| `out.style.figure_size.double_column` | `[17,6.65]` | cm `[width,height]` | Wide figure class. |
| `out.style.figure_size.tall_single` | `[12.5,10]` | cm `[width,height]` | Tall single-panel class. |
| `out.style.figure_size.tall_double_column` | `[17,12.5]` | cm `[width,height]` | Tall multi-panel class. |
| `out.style.layout.tile_spacing` | `'compact'` | char | MATLAB `tiledlayout` `TileSpacing`. |
| `out.style.layout.padding` | `'compact'` | char | MATLAB `tiledlayout` `Padding`. |
| `out.style.grid.alpha` | `0.4` | scalar `[0,1]` | Grid opacity. |
| `out.style.grid.line_style` | `':'` | char | MATLAB line style. |
| `out.style.tick_label_interpreter` | `'latex'` | char | Text interpreter for ticks/axes. |
| `out.style.hatch_spacing` | `0.06` | m scalar | Physical spacing of cross-section hatch lines. |

### 7.4 Status, material, and data colours

Every colour is an RGB row in `[0,1]`.

| Field | Default | Role |
|---|---|---|
| `out.style.status_palette.ok` | `[0.10,0.55,0.10]` | Acceptable status. |
| `out.style.status_palette.bad` | `[0.80,0.15,0.15]` | Failed status. |
| `out.style.color.cg` | `[1.00,0.00,0.00]` | Centre-of-gravity marker. |
| `out.style.color.cb` | `[0.00,0.00,1.00]` | Centre-of-buoyancy marker. |
| `out.style.color.series_a` | `[0.00,0.00,1.00]` | First comparison series. |
| `out.style.color.series_b` | `[1.00,0.00,0.00]` | Second comparison series. |
| `out.style.fill_palette.solid_material` | `[0.45,0.46,0.50]` | Solid structural material. |
| `out.style.fill_palette.jacket_material` | `[0.74,0.76,0.80]` | Jacket material. |
| `out.style.fill_palette.fill_material` | `[0.55,0.55,0.60]` | Filled-strip material. |
| `out.style.fill_palette.shell` | `[0.82,0.82,0.82]` | Shell annulus. |
| `out.style.fill_palette.void` | `[1.00,1.00,1.00]` | Void/air. |
| `out.style.fill_palette.hatch` | `[0.50,0.52,0.58]` | Void hatch. |
| `out.style.fill_palette.boundary` | `[0.10,0.10,0.10]` | Outer material boundary. |
| `out.style.fill_palette.inner_boundary` | `[0.30,0.30,0.32]` | Inner shell boundary. |
| `out.style.fill_palette.waterline` | `[0.15,0.55,0.95]` | Still-water line. |
| `out.style.fill_palette.fill_level` | `[0.95,0.55,0.10]` | Fill-level line. |
| `out.style.fill_palette.wall_boundary` | `[0.80,0.15,0.10]` | Precast wall boundary. |
| `out.style.color.dof` | `[0.12 0.47 0.71; 0.20 0.63 0.17; 0.89 0.10 0.11]` | Surge/heave/pitch primary colours by row. |
| `out.style.color.dof_secondary` | `[0.40 0.65 0.85; 0.55 0.78 0.52; 0.95 0.50 0.50]` | Surge/heave/pitch secondary colours. |
| `out.style.color.dof_pair` | `[0.58 0.40 0.74; 1.00 0.50 0.05; 0.55 0.34 0.29]` | DOF pairs `(1,3)`, `(1,5)`, `(3,5)`. |
| `out.style.color.dof_pair_secondary` | `[0.75 0.60 0.85; 1.00 0.73 0.47; 0.73 0.55 0.50]` | Secondary DOF-pair colours. |
| `out.style.color.reference` | `[0.3,0.3,0.3]` | Thresholds/asymptotes/targets. |
| `out.style.color.band` | `[0.85,0.92,1.0]` | Operational-period band. |
| `out.style.color.diagnostic` | `[0.90,0.20,0.10]` | Diagnostic highlight. |
| `out.style.color.mesh_edge` | `[0.35,0.35,0.35]` | BEM-panel edge. |
| `out.style.color.mesh_quad.fill` | `[0.65,0.82,1.0]` | Waterplane quadrilateral fill. |
| `out.style.color.mesh_quad.edge` | `[0.10,0.30,0.70]` | Waterplane quadrilateral edge. |
| `out.style.color.mesh_tri.fill` | `[0.65,0.95,0.75]` | Waterplane triangle fill. |
| `out.style.color.mesh_tri.edge` | `[0.05,0.45,0.20]` | Waterplane triangle edge. |

## 8. Required hydrodynamic-cache MAT schema

The cache file must contain one top-level variable, `hydro_table`, a scalar struct with all 25
fields below. The loader requests that variable by name; keeping it as the sole top-level variable
is the strict cache contract and avoids ambiguous auxiliary state. Validation is implemented in
`src/+mwecmass/+bem/load_hydro_cache.m`.

Let `N = numel(hydro_table.drafts)` and `M = numel(hydro_table.omega)`. Full matrices use the
6-DOF order `[surge,sway,heave,roll,pitch,yaw]`.

| Field | Required type / shape | Units and meaning |
|---|---|---|
| `drafts` | numeric vector, `N` elements | m; signed **vertical shifts**, despite the historic name. |
| `z_cg` | numeric vector, `N` elements | m, body-frame `z` up; coefficient reference CG for each node. |
| `submerged_volume` | numeric vector, `N` elements | m³, positive. |
| `displaced_mass` | numeric vector, `N` elements | kg, positive. |
| `omega` | numeric vector, `M` elements | rad/s; shared frequency grid. |
| `added_mass_inf` | cell vector, `N` cells; each `[6x6]` | Added mass at infinite frequency. |
| `radiation_damping_band_avg` | cell vector, `N` cells; each `[6x6]` | Radiation damping averaged over `period_band`. |
| `added_mass_omega` | cell vector, `N` cells; each `[6x6xM]` | Frequency-dependent added mass. |
| `radiation_damping_omega` | cell vector, `N` cells; each `[6x6xM]` | Frequency-dependent radiation damping. |
| `exciting_force_omega` | cell vector, `N` cells; each `[6xM]` | Complex excitation force/moment for one heading. |
| `period_band` | numeric `[1x2]` convention | s; band used for damping average. |
| `solver_params` | scalar struct | Solver-generation metadata; mixed units. |
| `water_depth` | numeric scalar | m, positive down; `-1` for deep water. |
| `ulen` | numeric scalar | m; BEM reference length. |
| `mesh_Nu` | numeric scalar | Hull mesh count in the first parametric direction; may be `NaN` when not applicable. |
| `mesh_Nv` | numeric scalar | Hull mesh count in the second parametric direction; may be `NaN` when not applicable. |
| `panel_size` | numeric scalar | m; generating panel-size target; may be `NaN` for a higher-order WAMIT cache. |
| `wp_target_edge` | numeric scalar | m; generating waterplane-lid edge target; may be `NaN` when not applicable. |
| `period_min` | numeric scalar | s; minimum analysis period. |
| `period_max` | numeric scalar | s; maximum analysis period. |
| `period_step` | numeric scalar | s; analysis-period increment. |
| `timestamp` | char row | Cache-build timestamp metadata. |
| `ms2_file` | char row | Geometry-deck identity metadata. |
| `ms2_date` | char row | Geometry-deck modification-date metadata. |
| `solver` | char row | Solver identity metadata. |

### 8.1 Matrix units

The physical unit of a matrix element depends on whether each DOF is translational or rotational.
For example, added-mass translation/translation terms are kg, translation/rotation terms are
kg·m, and rotation/rotation terms are kg·m². Radiation-damping and excitation units change
analogously between forces and moments. Cache values must already be dimensional SI values; the
runtime loader does not re-dimensionalize them.

### 8.2 Sorting, interpolation, and consistency

- The loader sorts `drafts` ascending and applies the same order to `z_cg`,
  `submerged_volume`, `displaced_mass`, and all five per-node cell arrays.
- Duplicate shifts are rejected because interpolation would be ambiguous.
- Every per-node array must have length `N`; every frequency-resolved item must use exactly `M`.
- The current cache contract represents one excitation heading: `[6xM]`, not `[6xMxH]`.
- The driver requires exact numeric equality between `in.bem.water_depth` and the cache's
  `water_depth`; there is no tolerance. Use the same stored scalar when configuring a run.
- Cache interpolation is performed in signed vertical shift and clamps queries to the tabulated
  shift range. A cache whose grid does not span the optimisation bounds is therefore stale for
  that design space even if it passes shape validation; regenerate or deliberately narrow the
  bounds.
- Treat changes to the hull deck, BEM reference length, water depth, period grid, mesh, solver,
  or coefficient reference convention as cache-invalidating changes. Metadata are recorded for
  provenance, but not every stale-cache condition can be inferred automatically.

A minimal inspection in MATLAB is:

```matlab
S = load(fullfile('Input','C1_wamit_cache.mat'));
assert(isfield(S,'hydro_table'));
fieldnames(S.hydro_table)
size(S.hydro_table.added_mass_omega{1})
```

## 9. MultiSurf `.ms2` geometry input

Place the deck named by `in.files.ms2_file` directly under `Input/`. It is a line-oriented ASCII
MultiSurf model bounded by `BeginModel;` and `EndModel;`; entity records terminate with `;`.
The parser reads `Units:`, `Extents:`, and `Symmetry:` headers. Coordinates are **not scaled** from
the `Units:` text, so supply geometry numerically in metres and use `Units: m kg` for an SI deck.
The deck axes become the body axes without remapping.

Supported surface definitions are:

| Surface | Runtime interpretation |
|---|---|
| `RuledSurf` | Linear interpolation between two boundary curves. |
| `RevSurf` | Profile curve revolved about a line axis over the declared angle interval. |
| `BLoftSurf` | B-spline loft through ordered section curves/snakes. |
| `DevSurf` | Ruled interpolation between a snake and a curve. |
| `MirrSurf` | Reflection of a source surface in `X=0` or `Y=0`. |

The parser also resolves supporting `FramePoint`, `MirrPoint`, `BCurve`, `Conic`, `CopyCurve`,
`Line`, `BSubCurve`, `Arc`, `PolyCurve2`, `ProjCurve`, `EdgeSnake`, `AbsBead`, `AbsRing`, and
`BSubSnake` entities. `RealList`, `Variable`, and `Pathname` are skipped; unknown keywords are not
geometry. Visible surfaces have positive visibility/layer values; hidden construction entities
may still be evaluated as dependencies.

Declared `Symmetry: x`, `Symmetry: y`, or both synthesizes reflected surfaces when no explicit
visible `MirrSurf` exists. Surface geometry should be closed or closeable and should represent a
consistent hull boundary. The authoritative vertical extent is sampled from visible surfaces;
the `Extents:` line is not trusted as the final hull `z` range.

## 10. Output files and directories

For a hull deck `<hull>.ms2` and realisation `<type>`, the standard layout is:

```text
Output/
├── <hull>_<type>_results.mat
├── preliminary/
│   ├── *.png / *.fig / *.pdf
│   └── *.log
├── thin_shell/
│   ├── *.png / *.fig / *.pdf
│   └── *.log
├── modular_precast/
│   ├── *.png / *.fig / *.pdf
│   └── *.log
├── diagnostics/                         # created only for enabled diagnostics
│   ├── WAMIT_GeomVerify_vs*.png
│   ├── UHPC_Mass_Balance.png/.fig
│   ├── Stage1_Draft_Screening.gif
│   └── Stage2_Mass_Distribution*.gif
└── hams_mrel/                           # HAMS route only
```

The result MAT file contains exactly the exported variables `results` and `final_props`. It is a
data-only export: parser objects are removed, a hull-deck SHA-256 digest is recorded when
available, and mode-specific absent fields are written as typed empties. The full field contract
is in [RESULT_SCHEMA.md](RESULT_SCHEMA.md).

Common figure stems are:

| Selection | Stem |
|---|---|
| Final cross-section | `WEC_Final_3D_CrossSection`, `_Steel`, or `_UHPC` |
| Stage-1 density | `WEC_Equivalent_Density_2D` |
| Stage-2 / realised density | `WEC_Equivalent_Density_3D` / `WEC_Realised_Density_3D` |
| Convergence | `WEC_Complete_Convergence` |
| Thin-shell solve | `Steel_Solve` |
| Precast views | `WEC_Constructability_XZ`, `WEC_Constructability_Strips` |
| Hydrodynamics | `WEC_HydroCoeffs`, `WEC_RAO`, `WEC_Ainf_vs_Draft` |
| Mesh / normals | `WEC_Mesh_Diagnostic`, `WEC_PanelNormals` |
| Draft landscape | `WEC_DraftLandscape` |

Logs are UTF-8 and are rewritten on each run:

| Log | Gate |
|---|---|
| `stage2_summary.log` | `out.save.stage2.summary_log` |
| `final_results.log` | `out.save.stage3.final_results_log` |
| `diagnostic_panels.log` | `out.save.stage3.diagnostic_panels_log` |
| `steel_solve.log` | Thin-shell mode and `out.save.stage3.steel_solve_log` |

## 11. Reproducible runs

For a traceable comparison:

1. Begin from a clean extraction of exactly one release ZIP. The complete archive is immediately
   reproducible with the supplied C1 deck/cache; the source-only archive requires a compatible deck
   and cache supplied by the user.
2. Record MATLAB release, operating system, and Optimization Toolbox version.
3. Preserve the exact `WEC_User_Input.m`, `WEC_Output_Options.m`, `.ms2` deck, and cache bytes.
4. Keep `out.export.timestamp_filenames=false` for stable filenames.
5. Run from the source root using an explicit `in` and `out` struct.
6. Retain the result MAT file and logs; the MAT export records the schema version, realisation
   type, context, runtime configuration, and hull-deck digest.
7. Compare numeric result fields, not raster pixels. Font availability, renderer, operating
   system, and MATLAB release can change figure appearance without changing the calculation.

If preserving the complete archive's reference outputs, copy the entire `Output/` directory to a
separate location before running, because standard filenames are deterministic and are overwritten.

## 12. Troubleshooting

| Symptom | Likely cause | Action |
|---|---|---|
| `Input/...ms2` not found | Source-only archive has no supplied deck, or filename changed | Place a compatible deck directly under `Input/` and match `in.files.ms2_file`. |
| Cache missing `hydro_table` or a named field | Wrong MAT file or obsolete cache schema | Rebuild/convert the cache to the 25-field schema in §8. |
| `DuplicateDraft` | Repeated signed shift nodes | Deduplicate and regenerate the corresponding per-node arrays. |
| `SizeMismatch` | `N` or `M` inconsistent across cells | Check every shape in §8; do not use multi-heading excitation arrays. |
| Water-depth mismatch | `in.bem.water_depth` differs from cache metadata | Use the cache's physical depth or regenerate coefficients for the intended site. |
| Query clamps at a cache endpoint | Optimisation bounds exceed the cache shift grid | Regenerate a wider cache or narrow `vertical_shift_bounds`. |
| Invalid realisation type | Typo or unsupported value | Use exactly one of the three strings in §3. |
| Thin-shell seed assertion | `rho_fill < rho_shell` | Restore the required density ordering or use a compatible material model. |
| Modular-precast strip assertion | Fewer than three strips | Set `num_ballast_sections >= 3`. |
| No hydrodynamic plots despite `true` flags | Dispatcher availability gate is unmet | Supply the configured HAMS executable where appropriate, or treat the cached coefficients through the exported MAT data. |
| No diagnostic files | Diagnostics are false by default, results MAT disabled, or mode gate unmet | Enable only the desired `out.save.diagnostics.*` flag and keep `out.save.results_mat=true`. |
| Figures use fallback fonts | Requested fonts are absent | Install Times New Roman/Courier New or choose installed fonts in output options. |
| Output files appear to be replaced | Stable filenames are the default | Enable timestamped figure names or archive `Output/` before the run. |
| HAMS launch fails | Executable absent, incompatible, or lacks permission/dependencies | Install HAMS-MREL separately and follow `HAMS_MREL_ROUTE.md`; it is not bundled. |
| Optimisation exits without satisfying all checks | Infeasible targets/bounds or poor conditioning | Inspect `stage2_summary.log`, `diagnostic_panels.log`, exit flag, constraint violation, and mass/GM margins before changing tolerances. |

## 13. Safe customization examples

Keep project defaults in the two root functions, but make exploratory changes on returned structs
so the exact variation is visible in the calling script.

The examples below assume MATLAB's current folder is the source root and `src/` has already been
added with `addpath(fullfile(pwd,'src'))`.

Run a preliminary calculation with only the essential record and logs:

```matlab
in = WEC_User_Input();
out = WEC_Output_Options();

in.materials.realisation_type = 'preliminary';
out.export.formats = {'png'};
out.save.stage2.cross_section = false;
out.save.convergence = false;
out.save.hydrodynamics.coefficients = false;
out.save.hydrodynamics.raos = false;
out.save.hydrodynamics.added_mass_vs_draft = false;
out.save.hydrodynamics.panel_normals = false;

mwecmass.driver.run(in, out);
```

Enable one post-processing diagnostic explicitly:

```matlab
in = WEC_User_Input();
out = WEC_Output_Options();
in.materials.realisation_type = 'modular_precast';
out.save.results_mat = true;
out.save.diagnostics.uhpc_mass_balance = true;
mwecmass.driver.run(in, out);
```

Change target bands while keeping exported context synchronized:

```matlab
in = WEC_User_Input();
out = WEC_Output_Options();
in.targets.T_heave_goal = 8.0;
in.targets.T_heave_range = [7.5, 9.5];
in.context.T_heave_goal = in.targets.T_heave_goal;
in.context.T_heave_range = in.targets.T_heave_range;
mwecmass.driver.run(in, out);
```

Manually restrict vertical shift only after confirming cache coverage:

```matlab
in = WEC_User_Input();
out = WEC_Output_Options();
in.bounds.vertical_shift_bounds = [-0.4, 2.5];  % signed shifts [m], not drafts
mwecmass.driver.run(in, out);
```

## 14. Model scope and operating conditions

- Results apply to the model, geometry, coefficient cache, targets, and numerical settings used. They are intended for preliminary design and do not constitute a structural-certification or seaworthiness assessment.
- The optimisation degrees of freedom are surge, heave, and pitch. Full cached matrices retain sway, roll, and yaw terms.
- The natural-period calculation uses a zero PTO-stiffness matrix and does not include mooring or PTO restoring models.
- Hydrodynamic data are interpolated across cached signed shifts and clamped outside the cache grid. Configure optimisation bounds within the available cache range.
- The cache schema stores one excitation heading.
- The MultiSurf parser reads the entity subset listed in §9. Supply geometry in SI units because the `Units:` header is recorded but does not rescale coordinates.
- The HAMS-MREL route uses a separately installed executable. It is optional for the supplied C1 example, which uses the WAMIT cache.
- Post-processing diagnostics are optional and disabled by default. They operate on explicitly supplied result data.
