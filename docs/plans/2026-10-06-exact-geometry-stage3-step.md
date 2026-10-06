# Exact geometry, Stage-3 rebuild and STEP export — implementation plan

> Read `AGENTS.md` first. Its ground rules (§1) bind every task in this plan. Open decisions
> OD1–OD6 (AGENTS.md §8) must be answered by the owner before the tasks that depend on them.

**Goal.** Remove every shape assumption from the geometry, mass and realisation code; rebuild
Stage 3 (modular precast and thin shell) on one exact geometry kernel; make Stage 3 follow the
owner's process; and write STEP files of the realised design.

**Architecture.** A new package `src/+mwecmass/+solid/` holds the exact geometry kernel:
outer-surface sampling and normals from `MS2Parser`, normal offset with fold trimming, spline
surface fitting, body assembly, slicing, and mass-property integration. Stage-2 floors,
both Stage-3 realisations, the stored contours, the figures and the STEP writer
(`src/+mwecmass/+output/+step/`) all call this kernel. Nothing else computes wall geometry.

**Tech stack.** MATLAB code (must also run in GNU Octave for testing), Octave test scripts in
`tests/`, Python `gmsh` (OpenCASCADE) as a test-only STEP checker.

---

## Execution model (subagent-driven)

- One implementation subagent per task, run one at a time unless a task is marked parallel.
  Implementers: **Sonnet 5.5, high effort** for kernel, solver and STEP tasks; **Sonnet 5.5,
  medium effort** for cleanup and docs tasks.
- After each task a **review subagent (Opus 5.5, high effort)** checks the diff against this plan
  and AGENTS.md §1, reruns the tests, and lists Critical / Important / Minor issues. Critical and
  Important issues are fixed by a follow-up subagent before the next task starts.
- Every task ends with its tests passing in Octave, a commit, and a push to
  `claude/lucid-cray-7o9442`. Work happens in the cloud container only.
- Checkpoints with the owner: after T3 (kernel validated), after T6 (UHPC Stage 3), after T10
  (STEP files). Report measured numbers, not claims.
- Deletions happen in the same commit as the code that replaces them. The pipeline must not sit
  in a broken state between tasks.
- Tests may assert exact-by-construction properties (symmetry, closure, shared edges, identity of
  two code paths). Approximation errors are measured and reported; a numeric pass/fail gate on
  an approximation needs owner approval (AGENTS.md §1 rule 5, OD5).

---

## Deletion register (what goes, and what replaces it)

| ID | Delete | Replaced by | Task |
|---|---|---|---|
| D1 | `hydrostatics/compute_perpendicular_shell_volume.m` | kernel wall volume at t_min | T4 |
| D2 | `geometry/compute_rmin_at_z.m` and every call | exact `t_max` from the kernel | T5, T7 |
| D3 | `build_config.m` Passes 1–3 (r_min profile, r′, s_max; dead) | nothing | T4 |
| D4 | `thin_shell/hull_slope_cos_at_z.m`; `t/cos α` and cap in `inner_properties_at_z.m` | kernel slices of the offset surface | T7 |
| D5 | `extract_strip_geometry.m` L_k block and `R_eq` fallback (138–168, 218) | kernel slices | T5 |
| D6 | `max_slope_factor` inputs, config fields and arguments | nothing | T5, T7 |
| D7 | `(A_in/A_out)^2` moment fallback in `inner_properties_at_z.m` | exact moments | T7 |
| D8 | void squeezes and wall lines in `plot_modular_precast.m`, `compute_inner_profile.m`, `plot_inner_spline.m`, `plot_steel_solve.m` offset, `stage_animations.m` squeeze | true sections of the realised solid | T8 |
| D9 | `strip_scale_factor`, `strip_r_min`, `feasibility.s_max`, `feasibility.rho_min_achievable` | nothing | T5, T11 |
| D10 | dead `config.shell` branches, `compute_strip_equivalent_density.m` | nothing | T8 |
| D11 | `internal/offset_polygon.m` if orphaned | nothing | T11 |
| D12 | doc sentences describing D1–D10 | new method text | T12 |
| D13 | `build_silhouette_profile.m` max-x mirror, `extract_midplane_profile.m` 5 cm band (figure inputs) | exact y = 0 slice | T8 |

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
`exitflag`/`output`); `tests/step_check.py` (gmsh import: solid count, open boundary edges,
mesh-divergence volume, bounding box, declared units); `tests/reference/c1_reference.py`
(independent Python evaluation of `Input/C1.ms2` — RevSurf as revolved spline, RuledSurf as
defined — producing reference sections, profile, normals and offset clearances);
`tests/README.md`. Shims never go on the production path.

Acceptance: install script reproduces the toolchain in a fresh container; `octave --no-gui
tests/run_tests.m` runs; `MS2Parser` parses C1 under Octave; a baseline file records the v1.0 C1
numbers quoted in AGENTS.md §2 for later comparison. Optionally add a SessionStart hook that runs
the install script in future cloud sessions.

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

- Inner nodes P_in = S − t·n on the T1 grid; remove nodes whose distance to the outer surface is
  below t (the fold, e.g. the C1 shoulder corner); fit a tensor-product B-spline surface through
  the remaining nodes (interpolating), per patch, per module range.
- Slice the fitted surface at any height into an ordered closed contour.

Acceptance: measured normal clearance of the fitted inner surface to the outer surface, reported
as min / max / distribution over a dense check grid not used for fitting; no self-intersection;
C1 at t = 0.0762 m compared with the Python reference half-widths per height and with the
independent erosion result for module 4 (void 3.548 m³ with the v1.0 module edges and z_fill).

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

Acceptance: identity checks (outer volume equals `compute_hull`/divergence result; symmetric CG
x = y = 0 to machine precision for C1); convergence of every quantity with quadrature order,
reported; agreement with an independent `gmsh` mesh integration of the same body, reported.
**Owner checkpoint after T3.**

### T4 — Stage-2 floors from the kernel; ρ_air naming (Sonnet high)

Files: `src/+mwecmass/+driver/build_config.m`, `WEC_User_Input.m`, delete D1, D3.

- ρ_min,i = [ρ_UHPC·V_wall,i(t_min) + ρ_air·(V_i − V_wall,i)] / V_i with V_wall from the kernel;
  `m_min_constructability` from the same values.
- Rename the modular-precast air density input to `rho_air` (keep the value 1.2); assert that
  every UHPC path reads that field; keep thin-shell `rho_fill` (solid fill) unchanged.

Acceptance: C1 floors reported (expected near 608/284/211/610 kg/m³ from the Python estimate);
no remaining caller of D1; Stage 2 still runs under the Octave shim.

### T5 — UHPC Stage 3, part a: split, build, check, store (Sonnet high) — needs OD2

Files: `src/+mwecmass/+realise/+modular_precast/` (`solve_and_extract.m`, new
`split_from_stage2.m`, `realise_modules.m`, `check_against_stage2.m`, rewrite of
`extract_strip_geometry.m`), `build_realised_properties.m`; delete D5, D9 (precast path),
D2 calls in `modular_precast/solve.m`, D6 (precast fields).

1. Read the whole Stage-2 solution (`Final3D`, `x_opt`).
2. V_UHPC,i = V_i(ρ_i − ρ_air)/(ρ_UHPC − ρ_air) per module; modules with ρ_i = ρ_UHPC are solid;
   the wall module stays solid.
3. Ballast module k* = lowest module with ρ_i < ρ_UHPC: walls at t_min, remaining UHPC as
   ballast from the module bottom (stays inside k*). Modules above k*: solve t_i ≥ t_min so the
   wall holds V_UHPC,i. A module whose V_UHPC,i is below its t_min wall volume is flagged
   (cannot occur once T4 floors are in place).
4. Evaluate the realised body with T3 at the Stage-2 draft: mass, z_CG, Iyy, GM, coupled periods
   (added mass at the realised CG). Compare with Stage 2 using `mass_acceptable_pct`.
5. Store the realised body description, exact contours, per-module volumes and masses, the check
   report and a status flag; `final_props` always describes this realised design.

Acceptance: C1 run under Octave prints the split, the realised module geometry and the check
table; `final_props` never contains Stage-2 values for this mode.

### T6 — UHPC Stage 3, part b: optimisation, spill, closest fail (Sonnet high) — needs OD1, OD6

Files: `src/+mwecmass/+realise/+modular_precast/solve.m` (rewrite), new `stage3_report.m`.

- Runs only when T5's check fails. Start point: T5's split. Variables: draft, fill level in k*,
  t_i of the hollow modules above k*. Equalities: flotation and GM = GM_Stage2. Objective:
  closeness to the Stage-2 solution (OD1). Bounds: t_min ≤ t_i ≤ t_max,i (largest t for which
  module i keeps a void, from the kernel); fill level within k*.
- If no feasible point exists, allow the fill to enter k*+1 (k* becomes solid; fill measured
  from the module bottom, OD6) and solve again.
- If still infeasible: keep the iterate with the smallest constraint violation, set status
  `failed`, and report each metric (value, Stage-2 value, deviation, limit, pass/fail) and the
  active reason. Plot it, store it in `final_props` and the `.mat`, export its STEP files.

Acceptance: C1 run under the Octave shim; report printed; a forced-infeasible test case (e.g. an
unreachable GM) produces a stored, plotted, flagged closest-fail design. **Owner checkpoint
after T6.**

### T7 — Thin-shell rebuild on the kernel (Sonnet high) — needs OD4

Files: `src/+mwecmass/+realise/+thin_shell/` (`build_geometry_grid.m`, `inner_properties_at_z.m`,
`solve.m`, `evaluate_design_point.m`, `strip_partition_volumes.m`); delete D2 (last callers),
D4, D6, D7.

- Keep the thin-shell formulation: variables `[vs, t, z_fill]`, period-penalty objective,
  flotation equality, GM ≥ `gm_min`, fill model (full section below `z_fill`).
- Replace the inner geometry with the kernel's normal offset at one global t, and `t_max` with
  the largest t for which the hull keeps a void over the middle 80 % of its height (same intent
  as today, exact geometry).

Acceptance: C1 thin-shell run under the Octave shim; no reference to D2/D4/D6/D7 remains; the
modular-precast path does not call thin-shell functions.

### T8 — Figures from the realised solid (Sonnet high)

Files: `src/+mwecmass/+output/+figures/` (new `realised_section_data.m` returning the y = 0 slice
and plan sections of the realised body; `plot_modular_precast.m`, `plot_steel_solve.m`,
`plot_optimised_cross_section.m`), `validation/diagnostics/stage_animations.m`; delete D8, D10,
D13.

- Figures draw true sections of the realised (or closest-fail) solid, show the ballast level,
  and label the status. Data functions run in Octave; drawing stays MATLAB code.

Acceptance: data functions tested in Octave against T3 sections; the owner confirms the figures
in MATLAB.

### T9 — STEP writer (Sonnet high) — parallel with T4–T7 (separate files)

Files: `src/+mwecmass/+output/+step/write_step.m` and helpers, tests.

- AP214 text writer: `B_SPLINE_SURFACE_WITH_KNOTS` lateral faces whose boundary rows lie in the
  bounding planes (no trimming curves), `PLANE` caps, shared `EDGE_CURVE`s, `CLOSED_SHELL`,
  `MANIFOLD_SOLID_BREP` and `BREP_WITH_VOIDS`, `OPEN_SHELL` / `SHELL_BASED_SURFACE_MODEL` for
  sheets, SI units in METRE, product and solid names, presentation layers per body (lessons from
  the earlier STEP pipeline: correct unit declaration, no free construction points, every entity
  on a named layer).

Acceptance: `tests/step_check.py` imports a unit cube, a cylinder-free test solid built from
spline faces, and a solid with a void; reports solid count, zero open edges, and volume against
the analytic value.

### T10 — Stage-3 STEP exports (Sonnet high) — needs OD3

Files: realisation `run.m` files, `WEC_Output_Options.m` (new `out.save.stage3.step_*`
switches), `src/+mwecmass/+output/+step/` builders.

- UHPC: one STEP per module; one STEP with all modules as one connected solid (built directly as
  one B-rep with the cavity as a void shell).
- Steel: ballast solid; shell as a 2D surface (OD3); combined file (OD3).
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
