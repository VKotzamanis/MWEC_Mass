# Exact geometry, Stage-3 rebuild and STEP export — implementation plan

> Read `AGENTS.md` first. Its ground rules (§1) bind every task in this plan. Every decision in
> AGENTS.md §8 is resolved.

**Goal.** Remove every shape assumption from the geometry, mass and realisation code; rebuild
Stage 3 (modular precast and thin shell) on one exact geometry kernel; make Stage 3 follow the
owner's process; and write STEP files of the realised design.

**Architecture.** A new package `src/+mwecmass/+solid/` holds the exact geometry kernel:
outer-surface sampling and normals from `MS2Parser`, normal offset with fold trimming, spline
surface fitting, body assembly, slicing, and mass-property integration. Stage-2 floors,
both Stage-3 realisations, the stored contours, the figures and the STEP writer
(`src/+mwecmass/+output/+step/`) all call this kernel. Nothing else computes shell geometry.

**Tech stack.** MATLAB code (must also run in GNU Octave for testing), Octave test scripts in
`tests/`, Python `gmsh` (OpenCASCADE) as a test-only STEP checker.

---

## Execution model (subagent-driven)

- **Implementers:** Sonnet 5.5, high effort for kernel, solver and STEP tasks; Sonnet 5.5, medium
  effort for cleanup and docs tasks. Opus 5.5, high effort is also allowed.
- **Grading:** every task is graded by an Opus 5.5, xhigh effort grader. The grader accepts only
  with no rule violation, every acceptance item reproduced by the grader itself, all tests
  passing, and a score ≥ 9 (the minimum of correctness, rules, evidence, quality, scope). Up to
  four grading rounds, with fixer subagents between rounds. Only accepted branches merge. A
  second grader checks each merge: the diff equals the accepted branches, the full test suite
  reruns, the push is verified.
- **Worktrees:** tasks run in git worktrees on task branches. The orchestrator merges accepted
  branches into `claude/lucid-cray-7o9442` and pushes after the full test suite passes. Work
  happens in the cloud container only.
- **Waves.** Tasks in one wave run in parallel; a wave starts when the previous one is merged.
  - A: T0 harness, T1 kernel rows and normals, T9 STEP writer.
  - B: T0c parser entity lookups, T0d geometry cache, T2 offset, fold trimming and spline fit.
  - C: T0b renames, T3 bodies and exact properties. **Owner checkpoint after T3.**
  - D: T4a, then T4b (both edit `build_config.m` and the Stage-2 files).
  - E: T5 then T6 (UHPC Stage 3), in parallel with T7 (thin shell). **Owner checkpoint after T6.**
  - F: T8 figures, T10 Stage-3 STEP exports. **Owner checkpoint after T10.**
  - G: T11 cleanup, then T12 docs, then T13 final review.
- **Baseline:** the Octave regression baseline (T0, `tests/baseline/`) may change only in a task
  that intends to change results (T4a, T4b, T5, T6, T7). That task updates the baseline in the
  same commit and prints the changed quantities. Renames (T0b), parser lookups (T0c) and the
  geometry cache (T0d) leave every number identical. Full-pipeline tests run only with
  `MWEC_REGRESSION=1` (a C1 run takes 24–94 min in Octave). Owner (2026-10-07): no
  whole-pipeline run before Waves B–E are implemented; components are tested on their own.
- Checkpoints with the owner report measured numbers, not claims.
- Deletions happen in the same commit as the code that replaces them. The pipeline must not sit
  in a broken state between tasks.
- Tests may assert exact-by-construction properties (symmetry, closure, shared edges, identity of
  two code paths). Approximation errors are measured and reported; a numeric pass/fail gate on
  an approximation needs owner approval (AGENTS.md §1 rule 5, OD5).

---

## Deletion register (what goes, and what replaces it)

| ID | Delete | Replaced by | Task |
|---|---|---|---|
| D1 | `hydrostatics/compute_perpendicular_shell_volume.m` | kernel shell volume at t_min | T4b |
| D2 | `geometry/compute_rmin_at_z.m` and every call | exact `t_max` from the kernel | T5, T7 |
| D3 | `build_config.m` Passes 1–3 (r_min profile, r′, s_max; dead) | nothing | T4b |
| D4 | `thin_shell/hull_slope_cos_at_z.m`; `t/cos α` and cap in `inner_properties_at_z.m` | kernel slices of the offset surface | T7 |
| D5 | `extract_strip_geometry.m` L_k block and `R_eq` fallback (138–168, 218) | kernel slices | T5 |
| D6 | `max_slope_factor` inputs, config fields and arguments | nothing | T5, T7 |
| D7 | `(A_in/A_out)^2` moment fallback in `inner_properties_at_z.m` | exact moments | T7 |
| D8 | void squeezes and shell lines in `plot_modular_precast.m`, `compute_inner_profile.m`, `plot_inner_spline.m`, `plot_steel_solve.m` offset, `stage_animations.m` squeeze | true sections of the realised solid | T8 |
| D9 | `strip_scale_factor`, `strip_r_min`, `feasibility.s_max`, `feasibility.rho_min_achievable` | nothing | T5, T11 |
| D10 | dead `config.shell` branches, `compute_strip_equivalent_density.m` | nothing | T8 |
| D11 | `internal/offset_polygon.m` if orphaned | nothing | T11 |
| D12 | doc sentences describing D1–D10 and D14–D17 | new method text | T12 |
| D13 | `build_silhouette_profile.m` max-x mirror, `extract_midplane_profile.m` 5 cm band (figure inputs) | exact y = 0 slice | T8 |
| D14 | `c_mono` in `optim/stage2_constraints.m` (lines 35, 39, 51; `rho_max` at line 15, used only by it; its share of the fallback size at line 64; the monotonicity wording in the comments at lines 3–4 and 10–13) and the Stage-1 copy `c_monotonic` in `optim/solve_2d_surrogate.m` (lines 224–238, its entry in the concatenation at line 248, its share of the fallback size at lines 261–264 and 269, and the comment at lines 189–192) | nothing: stability by GM ≥ `gm_min` (Stage 2) and GM = GM_Stage2 (Stage 3) | T4a |
| D15 | `c_mass_min` in both files of D14 (`stage2_constraints.m:42–49`, `solve_2d_surrogate.m:240–246`) and its fallback size (`stage2_constraints.m:60–63`, `solve_2d_surrogate.m:265–268`); `config.m_min_constructability` (`build_config.m:609, 647`) once no reader is left. The m_min sum stays for the feasibility check at `build_config.m:620–645` | nothing: implied by the density bounds | T4a |
| D16 | UHPC Stage-3 feasibility pre-check, `modular_precast/solve_and_extract.m:64–106` (0.95 and 1.05 factors, `MassTooLight` and `MassTooHeavy` errors) | the closest-fail rule; the Stage-2 floors from the kernel | T5 |
| D17 | Thin-shell `t_max` numbers in `thin_shell/solve.m`: the 10 %–90 % probe range (lines 76–81), the 0.95 · t_max bound (lines 123, 237), the 0.9 · t_max start clamp (line 110), the 1 % `t_min_active` flag (line 260, used at lines 296 and 416) | `t_max` = half the neck thickness from the exact geometry | T7 |

Kept as exact geometry: `MS2Parser` evaluation of `RevSurf` (revolved spline), arcs and conics;
`isocurve_revsurf.m`; `compute_wetted_surface_area.m`; mesh panel normals.

---

## Tasks

### T0 — Toolchain, test harness, baseline (Sonnet medium)

Already in place: `tools/install_toolchain.sh` (verified: Octave 8.4.0, gmsh 4.15.2) and
`tests/octave_shims/startsWith.m`, `endsWith.m` (Octave 8 single-space-pattern bug).

Files to create: `tests/run_tests.m`; further shims in `tests/octave_shims/` — `contains`,
`discretize`, `datetime`, `optimoptions` (returns a struct), and `fmincon` implemented on
Octave's core `sqp` (map `c ≤ 0` to `h = −c ≥ 0`, accept infeasible starts, return MATLAB-style
`exitflag`/`output`); `tests/reference/c1_reference.py`
(independent Python evaluation of `Input/C1.ms2` — RevSurf as revolved spline, RuledSurf as
defined — producing reference sections, profile, normals and offset clearances);
`tests/README.md`. Shims never go on the production path.

Acceptance: install script reproduces the toolchain in a fresh container; `octave --no-gui
tests/run_tests.m` runs; `MS2Parser` parses C1 under Octave; a baseline file records the v1.0 C1
numbers quoted in AGENTS.md §2 for later comparison. Optionally add a SessionStart hook that runs
the install script in future cloud sessions.

### T0b — Rename the ballast and density variables (Sonnet medium)

Wave C (after T0d, which also edits `build_config.m`), so every Stage-2 and Stage-3 task uses
the new names. Pure renaming: no change to any
formula or value; every number in the baseline stays identical. Covers `src/`, `WEC_User_Input.m`, `WEC_Output_Options.m`, `validation/`,
`Input/WAMIT/` if affected, `docs/`, `README.md`, `AGENTS.md`.

| Old name | New name | Meaning |
|---|---|---|
| `z_fill` (inputs, config, opts, local variables, `results.constructability.z_fill`, `results.steel_data.z_fill`, docs) | `z_ballast` | top of the solid ballast region |
| `in.materials.modular_precast.rho_fill`, `config.constructability_rho_fill`, `constructability.rho_fill` | `rho_air`, `config.constructability_rho_air`, `constructability.rho_air` | air in the precast voids |
| `in.materials.thin_shell.rho_void` | `rho_air` | air inside the steel shell (`config.rho_air` already) |
| `in.materials.thin_shell.rho_fill`, `config.rho_fill`, `steel_data.rho_fill` | `rho_ballast` | solid steel ballast density |
| thin-shell `V_fill`, `M_fill`, `z_cg_fill`, `strip_V_fill` and similar | `V_ballast`, `M_ballast`, `z_cg_ballast`, `strip_V_ballast` | ballast region quantities |
| precast `rho_steel` (holds the UHPC density), `V_steel`, `M_steel`, `t_steel` in the UHPC path | `rho_uhpc`, `V_uhpc`, `M_uhpc`, `t_uhpc` | UHPC quantities |
| `steel_data.rho_steel` (holds the fill density) | removed; use `rho_ballast` | — |

Acceptance: `grep` finds none of the old names in the live code or docs (old names may appear
only in a schema note that maps old `.mat` fields to new ones); all `.m` files still parse in
Octave; the T0 smoke tests pass; `RESULT_SCHEMA.md`, `export_schema.m` and
`check_export_schema.m` list the new field names.

### T0c — Parser: resolve each entity once (Sonnet high)

Wave B; owner request (2026-10-07). Files: `src/+mwecmass/+geometry/MS2Parser.m`, new tests and
fixtures.

- Every point evaluation looks up its curves and surfaces by name in a `containers.Map`. Octave
  profile of `build_config` (C1, thin shell): 2,067 s in total, of which `containers.Map`
  subsref, isKey and key encoding take 1,083 s over 14.5 million lookups, and `feval` dispatch
  54 s. Resolve each entity's references and its evaluator once (at parse time or on first use)
  and evaluate through them.
- No change to any formula, operation or operation order: every output stays bit-identical. The
  public API and every caller stay unchanged.

Acceptance: before the change, save the outputs of the unchanged parser to `tests/fixtures/`
(every visible surface and every curve on a dense parameter grid; `extract_isocurve_at_z` at 50
or more heights); after the change the same calls return identical arrays (`isequal`). The T0
regression baseline reproduces exactly. The Octave time of `build_config` and the top profile
entries are printed before and after.

### T0d — Save and reload the geometry products of `build_config` (Sonnet high)

Wave B; owner request (2026-10-07). Files: `src/+mwecmass/+driver/build_config.m`, one new helper
under `src/+mwecmass/+internal/`, `WEC_User_Input.m`, `.gitignore`, `tests/run_tests.m`.

- Save the products of the expensive geometry steps of `build_config` (boundary cache;
  waterplane, V_sub/CB_z and S_wet tables; strip geometry; Y-span table; per-strip density
  bounds; hull volume) to a `.mat` file and reload them when the fingerprint matches; otherwise
  recompute and overwrite. Print which of the two happened. Store plain arrays and structs only
  (Octave cannot save classdef objects).
- Fingerprint: a hash of the `.ms2` file bytes, of every input field these steps read, and of the
  source of every `.m` file they can call (whole folders are acceptable; a missing file is not).
  Hash with `java.security.MessageDigest`; Octave uses the T0 shim.
- Cache folder from a new input `in.files.geometry_cache_dir` (default `Output/cache/`, ignored
  by git); tests use a temporary folder. Later tasks that add expensive geometry steps to
  `build_config` (T4b) put them inside the cached block.

Acceptance: a reloaded config equals a freshly built one (`isequal`); changing the `.ms2` file,
any keyed input, or any keyed source file triggers a rebuild (one test each); the grader checks,
by reading the cached steps, that every input they read is in the key; the T0 regression
baseline reproduces exactly; build and load times are printed.

### T1 — Kernel A: outer surface rows and normals (Sonnet high)

Files: `src/+mwecmass/+solid/outer_rows.m`, `surface_normals.m`, tests.

- Sample each `.ms2` patch so that rows are horizontal sections at requested heights (module
  edges, fill level, waterline, quadrature nodes) and columns follow the patch parameter.
- Normals from S_u × S_v using `MS2Parser.eval_surface` (analytic derivative where the parser
  provides it, central differences otherwise), oriented outward; handle mirror patches, seams
  and the undefined normal at the keel point and the neck top.

Acceptance: sections equal `extract_isocurve_at_z` and the Python reference; normals are unit,
outward, continuous across patch seams; exact symmetry under the deck's mirror planes.

### T2 — Kernel B: normal offset, fold trimming, spline surfaces (Sonnet high)

Files: `src/+mwecmass/+solid/offset_surface.m`, `trim_fold.m`, `fit_bspline_surface.m`,
`eval_bspline_surface.m`, `slice_bspline_surface.m`, tests.

- Adaptive, error-bounded fitting as specified in AGENTS.md §5 item 9: initial nodes dense where
  curvature is high or the void is narrow; offset by t + ε/2 along the exact normal; trim the
  fold; split faces along creases; fit cubic B-splines; check M1–M3 on points between the nodes;
  insert knots only in failing spans; remove unneeded knots at the end; stop with an error report
  at the iteration cap.
- Slice the fitted surface at any height into an ordered closed contour.

Acceptance: M1–M3 pass on a dense check grid not used for fitting (report min / max t_local, the
number of refinement passes and knots per face); creases present where the C1 shoulder fold is
trimmed; C1 at t = 0.0762 m compared with the Python reference half-widths per height and with
the independent erosion result for module 4 (void 3.548 m³ with the v1.0 module edges and
`z_ballast`); a slender-section case (thin-shell neck, t = 0.025 m) passes M3.

### T3 — Kernel C: bodies and exact properties (Sonnet high)

Files: `src/+mwecmass/+solid/build_body.m`, `body_properties.m`, `hydrostatics_at_draft.m`,
tests.

- Body = outer surface between two heights, minus cavities bounded by per-module inner surfaces
  (t_i), the ballast level and module planes. Mass properties per region (UHPC / fill / air):
  volume, first moments, second moments about the origin and the CG, by integrating the exact
  sections with Gauss–Legendre quadrature in z on each smooth interval (breakpoints at module
  edges, fill level and fold-trim limits).
- Hydrostatics at a given draft from exact sections: V_sub, CB, A_w, I_wp (with the
  parallel-axis term about the centre of flotation).

Acceptance: volume closure to machine precision (AGENTS.md §1 rule 11) for every module and
for the whole hull, with the ballast level placed inside a module, exactly at a module edge, and
spilled into the next module; identity checks (outer volume equals `compute_hull`/divergence
result; symmetric CG x = y = 0 to machine precision for C1); convergence of every quantity with quadrature order,
reported; agreement with an independent `gmsh` mesh integration of the same body, reported.
**Owner checkpoint after T3.**

### T4a — Stage-2 changes that do not need the kernel (Sonnet high)

Files: `src/+mwecmass/+optim/stage2_constraints.m`, `solve_2d_surrogate.m`, `stage2_bounds.m`,
`run.m`, `src/+mwecmass/+driver/build_config.m`, `WEC_User_Input.m`; delete D14, D15.

- Delete `c_mono` and `c_mass_min` (D14, D15; AGENTS.md §3 items 26 and 33). Stability stays
  enforced by GM ≥ `gm_min`.
- Bottom-filled start (AGENTS.md §3 item 26): Stage 2 also runs from a start at the Stage-1 draft
  with every module at its lower density bound (its floor once T4b is merged), then modules
  filled to the mode's solid density from the keel up until mass = ρ_w · V_sub; the UHPC wall
  module stays pinned. Stage 2 keeps the start with the lower objective and logs both.
- Upper density bound per mode from the material inputs (OD9): `rho_ballast` (7500) for thin
  shell, `rho_uhpc` (2500) for modular precast.
- Set `in.materials.thin_shell.t_min` = `t_init` = 0.0254 m (owner decision, OD7) and update
  `RUNTIME_GUIDE.md` defaults.
- Assert that every UHPC path reads `rho_air` (1.2) and the thin-shell paths read `rho_air` and
  `rho_ballast` as named after T0b.

Acceptance: no reference to `c_mono`, `c_monotonic` or `c_mass_min` remains; Stage 2 runs under
the Octave shim for both modes and logs both starts; the baseline is updated in the same commit
and the changed quantities are printed.

### T4b — Stage-2 floors from the kernel, both modes (Sonnet high)

Files: `src/+mwecmass/+driver/build_config.m`, `src/+mwecmass/+optim/stage2_bounds.m`; delete D1,
D3.

- Modular precast: ρ_min,i = [ρ_UHPC·V_shell,i(t_min) + ρ_air·(V_i − V_shell,i)] / V_i with
  V_shell from the kernel; `m_min` from the same values, used by the feasibility check in
  `build_config.m`. Wall module pinned solid as today.
- Thin shell (new): the same floor with ρ_shell and the user-set minimum shell thickness, no wall
  module.
- `in.materials.modular_precast.t_init` is derived from the steel `t_init` today; the new UHPC
  Stage 3 starts from the Stage-2 split instead, so remove that derived input in T5.

Acceptance: floors reported for both modes (Python estimates on C1 — precast ≈ 608/284/211/610
kg/m³; thin shell at 25.4 mm ≈ 722/268/451/2009/2217); no remaining caller of D1; Stage 2 runs
under the Octave shim for both modes; the baseline is updated in the same commit and the changed
quantities are printed.

### T5 — UHPC Stage 3, part a: split, build, check, store (Sonnet high) — needs OD2

Files: `src/+mwecmass/+realise/+modular_precast/` (`solve_and_extract.m`, new
`split_from_stage2.m`, `realise_modules.m`, `check_against_stage2.m`, rewrite of
`extract_strip_geometry.m`), `build_realised_properties.m`; delete D5, D9 (precast path),
D2 calls in `modular_precast/solve.m`, D6 (precast fields), D16 (the pre-check in
`solve_and_extract.m`).

1. Read the whole Stage-2 solution (`Final3D`, `x_opt`).
2. V_UHPC,i = V_i(ρ_i − ρ_air)/(ρ_UHPC − ρ_air) per module; modules with ρ_i = ρ_UHPC are solid;
   the wall module stays solid.
3. Ballast module k* = lowest module with ρ_i < ρ_UHPC: the shell above the ballast starts at
   t_min (t_k* ≥ t_min is a variable in T6, like the t_i above k*), remaining UHPC as ballast
   from the module bottom (stays inside k*). Modules above k*: solve t_i ≥ t_min so the shell
   holds V_UHPC,i. The wall module and the ballast zone are full solid sections with no t_min
   shell added (AGENTS.md §1 rule 11). A module whose V_UHPC,i is below its t_min shell volume is
   flagged (cannot occur once the T4b floors are in place).
4. Restore equilibrium at the Stage-2 draft: adjust the ballast level in k* so that mass equals
   the Stage-2 displaced mass (differences come only from exact versus table module volumes).
   Evaluate the realised body with T3 at the Stage-2 draft: mass, Z_CG, Iyy, GM, coupled
   periods (added mass at the realised CG). Compare with Stage 2 using `mass_acceptable_pct`:
   |Z_CG,3 − Z_CG,2| / |Z_CG,2| with Z_CG = `CG_total(3)` (world frame), and the same relative
   form for GM, coupled T_heave and coupled T_pitch. Flotation is an equality held to the
   solver's constraint tolerance (OD10).
5. Store the realised body description, exact contours, per-module volumes and masses, the check
   report and a status flag; `final_props` always describes this realised design.

Acceptance: C1 run under Octave prints the split, the realised module geometry and the check
table; `final_props` never contains Stage-2 values for this mode.

### T6 — UHPC Stage 3, part b: optimisation, spill, closest fail (Sonnet high)

Files: `src/+mwecmass/+realise/+modular_precast/solve.m` (rewrite), new `stage3_report.m`.

- Runs only when T5's check fails. Start point: T5's split. Variables: draft, ballast level in
  k*, t_k* and the t_i of the hollow modules above k*. Equalities: flotation and
  GM = GM_Stage2. Objective (AGENTS.md §3 item 27): minimise Σ ((X3 − X2)/X2)² over Z_CG =
  `CG_total(3)`, coupled T_heave and coupled T_pitch (the coupled periods fix I20). Bounds:
  t_min ≤ t_k*, t_i ≤ t_max,i (largest t for which module i keeps a void, from the kernel);
  ballast level within k*.
- Escalation order (owner decision; equalities held to the solver's constraint tolerance, OD10):
  1. Draft fixed at the Stage-2 value; ballast within k*.
  2. Draft fixed; ballast may enter k*+1: the same ballast-level variable, its upper bound
     widened from the top of k* to the top of k*+1 (k* is then solid; OD6, option ii).
  3. Draft free (last resort): only when mass balance cannot be met at the Stage-2 draft after
     step 2. A failed `mass_acceptable_pct` check at the Stage-2 draft does not release it.
  4. Closest fail (below).
- After every solve the `mass_acceptable_pct` check on Z_CG, GM, T_heave and T_pitch decides
  accepted or failed.
- Closest fail: if the equalities hold but the check fails, the optimum itself is the closest
  design; if the equalities cannot be met, the iterate with the smallest equality violation is.
  Set status `failed`, and report each metric (value, Stage-2 value, deviation, limit,
  pass/fail) and the active reason. Plot it, store it in `final_props` and the `.mat`, export its
  STEP files.

Acceptance: C1 run under the Octave shim; report printed; a forced-infeasible test case (e.g. an
unreachable GM) produces a stored, plotted, flagged closest-fail design. **Owner checkpoint
after T6.**

### T7 — Thin-shell rebuild on the kernel (Sonnet high)

Also in T7: evaluate coupled periods in the objective (I20), and replace the fallback to Stage 2
with the closest-fail rule (OD4): flagged status, per-metric report, plotted, stored in
`final_props` and the `.mat`, STEP files exported. Escalation order: (t, z_ballast) at the
Stage-2 draft first; free the draft only if mass balance cannot be met there (OD8, OD10). A
failed `mass_acceptable_pct` check at the Stage-2 draft does not release the draft. Closest
fail as in T6.

Files: `src/+mwecmass/+realise/+thin_shell/` (`build_geometry_grid.m`, `inner_properties_at_z.m`,
`solve.m`, `evaluate_design_point.m`, `strip_partition_volumes.m`); delete D2 (last callers),
D4, D6, D7, D17.

- Keep the thin-shell formulation: one uniform thickness t for the whole hull and `z_ballast`
  free to pass module edges without penalty (plus the draft, last resort); objective: minimise
  Σ ((X3 − X2)/X2)² over Z_CG = `CG_total(3)`, coupled T_heave and coupled T_pitch (AGENTS.md
  §3 item 27); flotation equality; **GM = GM_Stage2** (OD11, replaces GM ≥ `gm_min`); ballast
  model (full section below `z_ballast`). At fixed draft the two equalities fix t and
  `z_ballast` (AGENTS.md §5, "What the two equalities imply").
- Replace the inner geometry with the kernel's normal offset at the single t. Build it with care
  in the slender neck (C1 half-width 0.10 m): the void there is 0.20 − 2t wide and closes at
  t = 0.10 m; the rounded top (radius 0.10 m) needs no fold trimming for t < 0.10 m. Remove the
  forced `cos α = 0.1` near the top.
- `t_max` = half the thickness of the slender wall (the neck), where the two offset shells meet,
  computed on the exact geometry (AGENTS.md §3 item 30): 0.10 m for C1. It replaces
  0.5 × `compute_rmin_at_z` (0.5025 m in the C1 run). The 10 %–90 % probe range, the
  0.95 · t_max bound, the 0.9 · t_max start clamp and the 1 % `t_min_active` flag are deleted
  (D17).

Acceptance: C1 thin-shell run under the Octave shim; realised module densities reported next to
the T4b floors; no reference to D2/D4/D6/D7/D17 remains; the modular-precast path does not call
thin-shell functions.

### T8 — Figures from the realised solid (Sonnet high)

Files: `src/+mwecmass/+output/+figures/` (new `realised_section_data.m` returning the y = 0 slice
and plan sections of the realised body; `plot_modular_precast.m`, `plot_steel_solve.m`,
`plot_optimised_cross_section.m`), `validation/diagnostics/stage_animations.m`; delete D8, D10,
D13.

- Figures draw true sections of the realised (or closest-fail) solid, show the ballast level,
  and label the status. Data functions run in Octave; drawing stays MATLAB code.

Acceptance: data functions tested in Octave against T3 sections; the owner confirms the figures
in MATLAB.

### T9 — STEP writer (Sonnet high) — wave A (separate files)

Files: `src/+mwecmass/+output/+step/write_step.m` and helpers, tests.

- Outer faces: exact NURBS conversions of the `.ms2` entities where the type allows (C1: B-spline
  curves, arcs, revolution, ruled surface); inner faces: the adaptive fits from T2.
- AP214 text writer: `B_SPLINE_SURFACE_WITH_KNOTS` (and the rational form for exact
  revolutions) lateral faces whose boundary rows lie in the bounding planes (no trimming curves), `PLANE` caps, shared `EDGE_CURVE`s, `CLOSED_SHELL`,
  `MANIFOLD_SOLID_BREP` and `BREP_WITH_VOIDS`, `OPEN_SHELL` / `SHELL_BASED_SURFACE_MODEL` for
  sheets, SI units in METRE, product and solid names, presentation layers per body (lessons from
  the earlier STEP pipeline: correct unit declaration, no free construction points, every entity
  on a named layer).

Acceptance: `tests/step_check.py` imports a unit cube, a cylinder-free test solid built from
spline faces, and a solid with a void; reports solid count, zero open edges, and volume against
the analytic value.

### T10 — Stage-3 STEP exports (Sonnet high)

Files: realisation `run.m` files, `WEC_Output_Options.m` (new `out.save.stage3.step_*`
switches), `src/+mwecmass/+output/+step/` builders.

- UHPC: one STEP per module; one STEP with all modules as one connected solid (built directly as
  one B-rep with the cavity as a void shell).
- Steel: ballast solid = full outer section below `z_ballast`; shell = the exterior parametric
  surface (exact NURBS) from `z_ballast` to the deck only (no double counting of the plate);
  combined file with both bodies sharing the identical junction curve at `z_ballast`. For C1 the
  split at `z_ballast` is exact (z depends only on the profile parameter of the RevSurf and the
  RuledSurf, so the cut is an iso-parameter line found by knot insertion).
- Files go to `Output/<type>/step/`; paths stored in the `.mat`.

Acceptance: C1 files pass `tests/step_check.py`; imported volumes equal the kernel volumes
(reported). **Owner checkpoint after T10.**

### T11 — Remaining cleanup (Sonnet medium)

`offset_polygon.m` if orphaned (D11); remaining D9 fields in `export_schema.m` and `Report.m`;
check that no shape-assumption pattern remains (`grep` for `sqrt(.*/pi)`, `pi\s*\*\s*r`,
`compute_rmin`, `R_eq`, `max_slope`, `cos_alpha`, `L_k`, `homothetic`).

### T12 — Documentation (Sonnet medium)

Rewrite the methods, runtime and schema sections listed under D12; document the kernel, the new
Stage 3, the status/report fields and the STEP outputs; update `CHANGELOG.md` and `AGENTS.md`.

### T13 — Final review (Opus high)

Whole-branch review against AGENTS.md §1 and §5; full Octave test run; a short list of what the
owner must confirm in MATLAB (real `fmincon` runs, figures).
