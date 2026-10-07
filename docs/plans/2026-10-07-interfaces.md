# Interface contract for the exact-geometry, Stage-3 and STEP tasks

> Binding for every task of `2026-10-06-exact-geometry-stage3-step.md` from the lanes on. Read
> `AGENTS.md` first. Built against the accepted T1 (`+solid/outer_rows.m`, `surface_normals.m`)
> and T9 (`+output/+step/write_step.m`, `validate_brep.m`, `tests/step_check.py`). Where a task's
> file list in the plan and this contract name different paths, the contract wins (listed in §3).
> Items are cited as S1–S8 (structures), F1–F14 (functions), I1–I9 (invariants).

## 0. Conventions

- **Units:** SI, lengths in metres, everywhere (structs, tables, STEP). Densities kg/m³.
- **Body frame:** the `.ms2` coordinates (C1: keel z = −3.25, top z = 1.10). **World frame:** the
  same x, y; z_world = z_body + vs, vs = `vertical_shift` = x(1); z_world = 0 is the still-water
  line, so the waterline is z_body = −vs and draft = −(z_min + vs). Every field name or comment
  states its frame; geometry structs are body frame, `CG_total`, `CB`, `KM` are world frame.
- **Names (rule 10):** regions are `uhpc`, `air` (modular precast) and `ballast`, `shell`, `air`
  (thin shell); densities `rho_uhpc`, `rho_air`, `rho_ballast`, `rho_shell`; ballast top
  `z_ballast`; never `z_fill`, `rho_fill`, `rho_void`, `rho_steel` or `V_steel` for UHPC.
- **Errors:** ids `mwecmass:solid:<Name>` (kernel), `mwecmass:realise:<Name>`, `mwecmass:step:<Name>`.
- **Tests:** `tests/<area>/test_<name>.m` (T0 conventions); exact-by-construction assertions only,
  approximation errors printed (rule 5).

## 1. Data structures

**S1 outer patch** (`geo.outer(k)`, one per entry of `model.visible_surfs`, same order; C1: 8)
- `name`, `source` (ultimate source entity), `type` (`RevSurf`, `RuledSurf`, …), `flips` (cellstr of
  mirror planes, from `cache.mirrors(m).effective_flips`).
- `surf`: exactly T9's `'bspline'` surface: `type='bspline'`, `degree [du dv]`, `ctrl [nu x nv x 3]`
  (m, body), `knots {ku, kv}` (clamped; **not renormalised** after a split, so a piece keeps its
  parent's parameter), `weights` (`[]` or `[nu x nv]`). It can be put into `brep.surfaces` as is.
- `outward`: true when S_u × S_v points out of the hull (mirror patches reverse it).
- `exact`: true when `surf` is the `.ms2` entity converted exactly (C1: all 8).
- `z_of_u`: true when z depends on u only (§5 I3); `u_range [u0 u1]`, `z_range [z(u0) z(u1)]`.
- `pole [1x2]` logical: the row at u0 / u1 collapses to one point (C1 RevSurf: keel and neck top).
- `seam_v0`, `seam_v1`: `[patch side]` of the neighbour sharing that v-boundary (side 0 or 1).

**S1b hull geometry** `geo`: `hull_name` (deck file stem, `C1`), `outer` (S1 array), `z_range
[z_min z_max]` (body), `analytic` (`[]` for a real deck; set only by stand-in fixtures, §3).
Built once per deck; stored as `config.hull_solid` by T4b (inside T0d's cached block).

**S2 inner set** (`inner`, one per distinct thickness t): `t` [m], `eps` = 0.01·t_min [m] (rule 5),
`z_range` (body, the slab it was fitted for), `patches` (S1 array; `outward` true when the normal
points **into the void**, i.e. out of the solid; `exact` false; crease pieces are separate patches),
`report` (S2r).
**S2r fit report:** per inner patch `n_nodes`, `n_knots [u v]`, `n_passes`, `n_removed`,
`n_check`, `t_local_min`, `t_local_max` [m], `M1`, `M2`, `M3` (logical), `M3_reason` (char); totals
`ok`, `cap_reached`. Check points: the nodes of the quadrature rule of F6 in every knot span plus
span midpoints (the fit is verified where the mass integrals sample it).

**S3 design** (input of F5): `mode` (`modular_precast` | `thin_shell`), `edges [N+1 x 1]` module
edges (body; `config.strip_edges`), `vs` [m], `t [N x 1]` [m] (thin shell: all equal; NaN for a
module with no void), `z_ballast` [m, body] (≤ z_min: no ballast; one value, may cross module edges),
`solid_modules` (indices solid regardless of `z_ballast`: the UHPC wall module; `[]` thin shell).
Module i layout (rule 11): full outer section solid below `z_ballast`; above it shell t(i) + air;
modules in `solid_modules` solid throughout. No horizontal slab at module joints.

**S4 body** (output of F5)
- `design` (S3), `inner_t` (the t of each S2 set used), `planes` (sorted split heights: module edges
  and `z_ballast` inside (z_min, z_max)).
- `brep`: exactly the T9 writer struct (`vertices`, `curves`, `edges`, `surfaces`, `faces`,
  `bodies`, `uncertainty` unset = 1e-7 m). Extra face fields (ignored by `write_step` and
  `validate_brep`, which read only the named fields): `role` (`outer` | `inner` | `cap` |
  `ballast_top` | `step`), `module` (index), `inside`, `outside` (region name on the side opposite
  to / along the face normal; `exterior` outside the hull). Face normal = surface normal × (−1 if
  not `same_sense`). Lateral faces are pieces of S1/S2 patches cut at `planes` by knot insertion;
  caps are `plane` surfaces at a plane strictly inside (z_min, z_max) (annuli carry a hole loop);
  `step` faces are the annuli between two inner sections of different t at a joint; a pole is a
  vertex, not an edge (the face loop has three edges).
- `brep.bodies`: modular precast `<hull>_UHPC_module_<i>` (solid, one closed shell each; a hollow
  module is a ring); thin shell `<hull>_STEEL_ballast` (solid, full outer section below
  `z_ballast`) and `<hull>_STEEL_shell` (sheet: outer faces above `z_ballast`).
- `shells` (signed face lists for F11 and tests): precast `all_outer`, `all_voids{c}` (the fused
  `<hull>_UHPC_all` solid with its cavities); thin shell `layer`, `void` (closed, mass regions only).

**S5 section loop** (F4, F6b): `z` [m, body], `pieces(k)` with `patch` (index), `u` (iso-u value),
`curve` (T9 curve struct: the exact iso-u row, running in v), `dir` (±1); `pts [n x 3]` CCW samples
on the exact curves (no repeated first point); `area`, `centroid [1x2]`, `I [Ixx Iyy Ixy]` about
the origin [m⁴] by Green's theorem on the curves (Gauss on knot spans); `simple` (logical).
Body section (F6b): `z`, `module`, `outer` (S5), `inner` (S5 or empty), `solid` (true below
`z_ballast` or in a solid module).

**S6 body properties** (F6): for each region r and module i: `V(i)` [m³], `S(i,:)` = ∫x dV [m⁴],
`J(:,:,i)` = ∫x xᵀ dV [m⁵], body frame about the origin; `regions.<r>` holds these; `modules(i)`:
`V`, `V_<r>`, `mass`, `rho_eff` = mass/V, `CG_body [1x3]`; `total`: `mass`, `CG_body [1x3]`,
`I_origin [3x3]`, `I_cg [3x3]` (I = ∫ρ(|x|²E − x xᵀ)dV); `quad` (Gauss order used).

**S7 hydrostatics** (F7): `vs`, `draft`, `V_sub` [m³], `CB_body`, `CB` (world) [1x3], `Aw` [m²],
`xF`, `yF` (centre of flotation), `I_wp_xx`, `I_wp_yy` (about the CF; fixes I12), `I_wp_yy_origin`,
`BM_L` = I_wp_yy/V_sub, `KM` (world z) = CB(3) + BM_L, `waterline` (S5 at z_body = −vs).

**S8 realised design** (`results.stage3`, both modes; replaces `results.constructability` and
`results.steel_data`)
- `mode`, `hull_name`, `status` (`accepted` | `failed`), `reason` (char, '' when accepted),
  `step` (`split` | `fixed_draft` | `spill` | `draft_free`: the last escalation step run), `vs`,
  `draft` [m] (world).
- `stage2`: `vs`, `rho [N x 1]`, `mass`, `Z_CG` (= `Final3D.CG_total(3)`, world), `GM`
  (`Final3D.GM_L`), `T_heave`, `T_pitch` (`Final3D.periods.heave/.pitch`, coupled).
- `rho` (densities used, by region), `design` (S3), `k_star` (precast: ballast module), `V_uhpc_target
  [N x 1]` (precast split, AGENTS §3 item 4.1).
- `modules(i)`: `z_lo`, `z_hi` [m, body], `t` [m] (NaN: no void), `h_ballast` [m] (ballast height
  measured from the module bottom, OD6 ii), `V`, `V_<region>`, `mass`, `rho_eff`, `rho_stage2`,
  `rho_floor`, `CG_world [1x3]`.
- `props`: F9 output (the `final_props` fields). `check`: F10 output. `solver`: per step run
  `exitflag`, `iterations`, `fval`, `max_eq_violation`. `fit`: S2r of every inner set. `body`: S4.
- `step_files`: struct array `name`, `path`, `bodies` (filled by F11; empty before).
`final_props` = `props` plus `stage3_status` (= `status`) and `stage3_check` (= `check`).

## 2. Functions

All take and return the structs of §1; none reads globals or writes files except F11.

| | Signature | Contract |
|---|---|---|
| F1 | `geo = mwecmass.solid.outer_nurbs(model)` | S1b from the `.ms2` entity tree: BCurve directly, BSubCurve by knot insertion, Arc as rational quadratic, PolyCurve2 joined with C0 knots, ProjCurve by projecting control points, RevSurf as profile × rational arc, RuledSurf as degree 1 in v between its two curves when both get the same reparameterisation (one knot vector, equal weights; C1: `curve7` is the projection of `Edge_For_Dev`), mirrors by flipping control points. Seam rows bitwise equal. Any other entity: error `NotExact` (Open item 1). |
| F2 | `[inner, rep] = mwecmass.solid.offset_surface(model, cache, geo, t, z_range, opts)` | S2 at thickness t over `z_range` (body). Nodes from `MS2Parser` + T1 `surface_normals`, offset by t + ε/2, folds trimmed (`trim_fold.m`), faces split at creases, cubic fit (`fit_bspline_surface.m`), knot insertion where M1–M3 fail, knot removal while they hold (AGENTS §5 item 9). A `z_of_u` outer patch gives `z_of_u` inner patches with monotonic z (the offset u-curves are fitted and the patch's v-construction is kept). Neighbouring inner patches share boundary curves bitwise. `opts.t_min` (ε), `opts.max_passes`; cap reached → error `FitNotConverged` listing the failing patches and metrics. |
| F2b | `t_max = mwecmass.solid.void_closing_thickness(model, cache, geo, z_range)` | Smallest t at which offset layers from opposite sides meet inside `z_range` (not a fold of one layer). C1 neck: 0.10 m. |
| F3 | `[S, Su, Sv] = mwecmass.solid.eval_bspline_surface(surf, u, v)`; `[C, Cs] = …eval_bspline_curve(curve, s)` | Rational or not; u, v column vectors; [n x 3]. |
| F3b | `[lo, hi, u_star] = mwecmass.solid.split_bspline_surface(patch, z)` | `z_of_u` patch, z strictly inside its z range: u* from z(u*) = z, knot insertion to multiplicity du; `lo`, `hi` keep the parent parameter and share the cut row bitwise. Else error `ZNotOneParameter` / `ZOutside`. |
| F4 | `loop = mwecmass.solid.slice_bspline_surface(patches, z)` | S5 from the iso-u rows; one closed simple loop or error `SectionNotClosed` (as T1). |
| F5 | `body = mwecmass.solid.build_body(geo, design, inner)` | S4. `inner` holds one S2 set per distinct t whose `z_range` covers its modules (callers cache them, §7). Builds every face once; edges shared by index. |
| F6 | `bp = mwecmass.solid.body_properties(body, rho, opts)` | S6 by the divergence theorem with fields (f,0,0): V = ∮x n_x, ∫x = ∮x²/2 n_x, ∫y = ∮xy n_x, ∫z = ∮xz n_x, ∫x² = ∮x³/3 n_x, ∫y² = ∮xy² n_x, ∫z² = ∮xz² n_x, ∫xy = ∮x²y/2 n_x, ∫xz = ∮x²z/2 n_x, ∫yz = ∮xyz n_x. Horizontal caps contribute exactly 0, so only lateral faces are integrated; a region's integral sums faces with `inside` = r minus faces with `outside` = r. Gauss–Legendre of order `opts.n_gauss` on every knot span (exact for polynomial faces, convergent for rational ones). |
| F6b | `sec = mwecmass.solid.body_section(body, z)` | Body section (S5 block). |
| F7 | `hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, opts)` | S7 on the exact outer patches cut at z_body = −vs (F3b, F6 on the outer pieces, F4 for the waterplane). |
| F8 | `fl = mwecmass.driver.density_floors(model, cache, geo, edges, t_min, rho_solid, rho_air, solid_modules)` | Per module ρ_min = [ρ_solid·V_solid + ρ_air·V_air]/V with a t_min shell, no ballast (F2, F5, F6); solid modules give ρ_solid. Returns `rho_min`, `V`, `V_solid`, `V_air` [N x 1], `fit` (S2r). |
| F9 | `props = mwecmass.realise.evaluate_realised(bp, hs, vs, config)` | `final_props` fields of `export_schema.m` from the exact body: `mass_total`, `CG_total` = [0 0 CG_body(3)+vs] (declared symmetric model; kernel x, y kept in `bp`), `Inertia_Tensor` = `I_cg`, `Ixx/Iyy/Izz`, S7 fields, `GM_L` = KM − CG_total(3), `K_hydro` (K33 = ρ_w g A_w, K55 = M g GM), A from `mwecmass.bem.interpolate_at_draft(vs, config, CG_total(3))`, coupled periods from `coupled_periods_by_share` (as `build_realised_properties.m` today), `mass_discrepancy`. |
| F10 | `check = mwecmass.realise.check_against_stage2(props, stage2, pct, tol_eq, rho_water)` | `metrics(k)` for `Z_CG`, `GM`, `T_heave`, `T_pitch`: `value`, `stage2`, `rel_dev` = abs(X3−X2)/abs(X2), `limit` = pct/100, `pass`; `equalities(k)`: `flotation` (M/(ρ_w V_sub) − 1, as Stage 2) and `GM` (GM/GM2 − 1), `residual`, `tol` = `tol_eq` (the Stage-3 fmincon ConstraintTolerance, 1e-6 today), `pass`; `pass`, `failed` (cellstr), `reason`. |
| F11 | `files = mwecmass.output.step.export_stage3(realised, out_dir)` | Writes with T9 `write_step` into `Output/<type>/step/`: precast `<hull>_UHPC_module_<i>.step`, `<hull>_UHPC_all.step`; thin shell `<hull>_STEEL_ballast.step`, `<hull>_STEEL_shell.step`, `<hull>_STEEL_all.step` (both bodies, junction curve shared). Returns `step_files`. |
| F12 | `data = mwecmass.output.figures.realised_section_data(realised, z_plan, n_z)` | y = 0 elevation (points where F6b loops cross y = 0, found on the exact curves), plan loops at `z_plan`, `z_ballast`, module edges, waterline, `status`. Body frame plus `vs`. |
| F13 | `plot_modular_precast(realised, config)`, `plot_steel_solve(realised, config)`, `plot_optimised_cross_section(final_props, config, x_opt, realised)` | Draw F12 data only. |
| F14 | `[results, final_props] = mwecmass.realise.<type>.run(config, x_opt, opt_results)` | Signature unchanged. Sets `results.stage3` (S8) and returns `final_props`; calls F13 and F11 (`config.output.save.stage3.step`) for every status, closest fail included. |

## 3. Ownership and stand-ins

| Item | Implements | Consumers | Stand-in before the join |
|---|---|---|---|
| S1, S1b, F1, F3, F3b, F4 | T2 | T3, T4b, T5, T7, T8, T10 | SK fixtures (below) |
| S2, S2r, F2, F2b | T2 | T3, T4b, T5, T6, T7 | SK |
| S3–S6, F5, F6, F6b | T3 | T4b, T5, T6, T7, T8, T10 | SK |
| S7, F7 | T3 | T5, T6, T7 | SK |
| F8 | T4b | build_config, `stage2_bounds.m` | — |
| S8, F9, F10, schema of S8 | T5 | T6, T7, T8, T10 | SK `sti_realised`, `sti_stage2` |
| F11 | T10 | F14 | none: call gated off in tests |
| F12, F13 | T8 | F14, `dispatch.m` | none: figure switches off in tests |
| F14 | T5 (precast; T6 adds the escalation inside `solve.m`), T7 (thin shell) | driver | — |
| T1 `outer_rows`, `surface_normals` | T1 (merged) | T2, T3 tests | — |
| T9 `write_step`, `validate_brep`, `step_check.py` | T9 (merged) | T3 tests, T10 | — |

Paths that differ from the plan's file lists: F10 lives in `+realise/check_against_stage2.m` (not
`+modular_precast/`), because both modes use it (§5 item 11); F9 is the new
`+realise/evaluate_realised.m` and replaces both `build_realised_properties.m`; the F14 calls of
F11 and F13 are written by T5 and T7, so T8 and T10 do not edit `run.m`.

**Stand-in kit SK** (one graded step at the start of lane K, merged before any consumer starts):
- `tests/standins/fixtures/` (permanent): `cylinder.ms2` (RevSurf of three straight lines: bottom
  disk, wall, top disk; sharp rims) and `box.ms2` (RuledSurf faces); `sti_closed_form.m` (closed-form
  V, S, J, sections, hydrostatics, `t_max` for any S3 design on these: inner cylinder of radius
  R − t between z0 + t and z1 − t, inner box shifted by t); `sti_config.m` (minimal config: `RHO_WATER`,
  `G`, one-draft hydro table, `strip_edges`, `pid.mass_acceptable_pct`); `sti_stage2.m` (Final3D-like
  struct consistent with the closed form); `sti_realised.m` (S8 for a box or cylinder design).
- `tests/standins/+mwecmass/+solid/` and `+realise/`: closed-form versions of F1–F7, F9, F10 that
  set or require `geo.analytic` (error `mwecmass:standin:NotAnalytic` otherwise).
- Consumer tests call `addpath(fullfile(root,'tests','standins'),'-end')`: a real function in `src/`
  wins (Octave 8.4 and MATLAB merge package folders, earlier path entries first; verified in Octave).
  The join that merges a producer deletes its stand-in functions; the consumer tests then run
  unchanged on the real code (the join test). Box results are exact with the real kernel; cylinder
  results (rational faces) are compared and printed.
- C1 pieces known by construction (oracles, one line in each test): for z ∈ [−0.5, 1.0] the section
  is a stadium of radius 0.10 about (0, ±1) with flat sides x = ±0.10 (curve1's first span has
  collinear control points), area 0.4 + 0.01π m², inner offset radius 0.10 − t; z ∈ [1.0, 1.1] is a
  half cylinder (r 0.10, length 2) with two quarter spheres, V = 0.01π + (2/3)π·0.001 m³; neck
  `t_max` = 0.10 m.

## 4. File ownership map (files more than one task edits; order of edits)

A task edits a shared file only after every earlier task in its row is merged (or, inside the J2
group, accepted) and it has rebased on it; such edits are the task's last commit.

| File | Order |
|---|---|
| `src/+mwecmass/+driver/build_config.m` | T0b → T0d → T4a → T4b → T5 (precast D6 forwards) → T7 (thin-shell D6 forwards) |
| `WEC_User_Input.m` | T0b → T0d → T4a → T5 (precast `t_init`, `max_slope_factor`, `n_z_grid`) → T7 (thin-shell ones) |
| `WEC_Output_Options.m` | T0b → T10 (`out.save.stage3.step`) |
| `src/+mwecmass/+driver/run.m` | nobody (F11 uses `mwecmass.output.output_dir`) |
| `optim/stage2_bounds.m` | T4a → T4b |
| `optim/run.m`, `optim/report_assemble_results.m` | T4a (`run.m`) → T5 (`stage3` replaces `constructability`/`steel_data`) |
| `output/export_schema.m`, `check_export_schema.m`, `export_results.m`, `load_results.m`, `Report.m` | T0b → T5 (S8 for both modes) → T11 (D9 remnants) |
| `docs/RESULT_SCHEMA.md` | T0b → T5 → T12 |
| `docs/RUNTIME_GUIDE.md` | T0b → T4a → T12 |
| `docs/METHODS_ENGINE.md`, `README.md`, `AGENTS.md`, `_graph/CODE_MAP.md` | spec → T0b → T12 |
| `+realise/+modular_precast/*` | T0b → T5 → T6 |
| `+realise/+thin_shell/*`, `+realise/empty_realised_properties.m` | T0b → T7 |
| `+realise/build_realised_properties.m` | T0b → T5 (stops calling it) → T7 (deletes it) |
| `+output/+figures/plot_modular_precast.m`, `plot_steel_solve.m`, `validation/diagnostics/stage_animations.m` | T0b → T8 |
| `output/dispatch.m` | T8 (F13 call) |
| `validation/diagnostics/uhpc_mass_balance.m` | T0b → T11 |
| `tests/run_tests.m`, `.gitignore` | T0 → T0d |
| `tests/baseline/*.json` | T0 → J2 run step |
| `MS2Parser.m` | T0c only. `+solid/*`: T2 and T3 own disjoint files. `+output/+step/*`: T9 → T10 (new files; `write_step.m` changes only if a C1 body fails `step_check`, by T3). |

## 5. Invariants every implementation tests

- **I1 Volume closure (rule 11)**, relative bound of a stated multiple of eps justified in the test:
  per module V_uhpc + V_air = V_module (precast); V_ballast + V_shell + V_air = V_module and = V_hull
  (thin shell); Σ_i V_module,i = V_hull, with V_hull integrated on the unsplit outer patches using
  the cut parameters u* as extra breakpoints (split pieces keep the parent parameter, so both sample
  the same (u, v) points). With `z_ballast` inside a module, at a module edge, and spilled.
- **I2 Shared edges:** every closed shell passes `validate_brep`; the boundary rows of the two
  faces at an edge equal its curve bitwise (control points, knots, weights).
- **I3 Exact z-splits** where z depends on u only. RevSurf about an axis parallel to z: rotation
  keeps z, so every control row has one z and the weights are separable (w_ij = a_i b_j), hence
  z(u,v) = Σ N_i a_i z_i / Σ N_i a_i. C1 RuledSurf: `curve7` = `ProjCurve` of `Edge_For_Dev` onto
  y = 0 keeps z, so both ends of every ruling have the same z. A plane z = z_k then meets the patch
  in the iso-u row u*, and knot insertion cuts it without moving the surface (test: parent and
  pieces agree at the same (u, v) to rounding; cap vertices within T9's 1e-7 m). The offsets keep
  it: the normal of a surface of revolution lies in the meridian plane, and C1's rulings are all
  parallel to y, so the normal is constant along a ruling.
- **I4 Units:** METRE in every struct; `step_check` reports `declared_length_unit = METRE`.
- **I5 Names:** S8, S6 and new code use only the names of §0; `grep` for the old names in new
  files is empty.
- **I6 Frames:** CG_total(3) = CG_body(3) + vs; KM, CB world; tests check both frames.
- **I7 Orientation:** outer normals out of the hull, inner normals into the void; every region
  volume positive.
- **I8 C1 symmetry:** CG x = y = 0 and products of inertia 0 to machine precision.
- **I9 Exactness of F1:** F4 sections of the NURBS and T1 `outer_rows` at the same z agree to
  rounding (both exact).

## 6. Lanes and joins

```
J0   merge T0, T1, T9, this contract
K    SK ─► T2 ─► T3 ─────────────────────► J1 (T2+T3)  ── owner checkpoint after T3
N    T0b ─► merged at once (all later lanes start from it)        T0c (MS2Parser.m) in parallel
P    [T0b, T0c merged] ─► T0d ─► T4a ─► T4b (stand-ins; merges after J1)
U    [T0b merged] ─► T5 ─► T6
S    [T0b merged] ─► T7
O    [T0b merged] ─► T8 ∥ T10
J2   after J1 and lane P: the group T5, T6, T7, T8, T10 merges as one chain in that order
     ─► first whole-pipeline runs ─► owner checkpoints after T6 and after T10
G    T11 ─► T12 ─► T13
```
- Tasks are graded on component tests (stand-ins, C1 pieces, real kernel once J1 is in).
  Acceptance items that need a whole-pipeline run are checked only after J2 (owner: no
  whole-pipeline run before T5–T7 are implemented): (1) the T0 regression baseline on the
  integration commit after T0d (identity of T0b, T0c, T0d); (2) C1 runs of both modes on the J2
  head for T4a, T4b, T5, T6, T7, T10; then the baseline is regenerated with
  `tools/write_octave_baseline.m` and the changed quantities are printed per task. A failure there
  is fixed on a fix branch before anything else merges.
- Inside J2 the integration branch is pushed only after the last merge and the full suite.

## 7. Stage-3 evaluation per solver iteration, and what must be cached

Per evaluation of (t_k*, t_i, `z_ballast`, vs) or (t, `z_ballast`, vs): S2 sets for each distinct
t (F2, the expensive step: parser normals, fit, M1–M3), body assembly (F5: knot insertion at
`z_ballast`), F6 on the faces that changed, F7 (depends on vs only), F9 (BEM interpolation and the
3 × 3 eigenproblem), F10. Needs:
1. S1b, the outer pieces cut at the module edges and their F6 integrals: once per deck and edges.
2. S2 sets keyed by the exact t value and fit range: a finite-difference step in t must not refit an
   unchanged t, and each changed t is fitted once.
3. F7 keyed by vs: one call while the draft is fixed (steps `split`, `fixed_draft`, `spill`).
4. Smoothness: the properties must be differentiable in t and `z_ballast` within one solve, so the
   knot vectors of each S2 set stay fixed during a solve (control points refitted for a new t);
   the adaptive fit with M1–M3 runs on the final t, and a changed knot structure restarts the solve
   from that point. `z_ballast` enters only through knot insertion, which is smooth.

## 8. Decisions taken here (simplest option that meets the rules)

- One B-rep (S4) for mass properties, sections, figures and STEP: one geometry (§5 item 1, I6).
- Volume integrals use fields with no z component, so caps add 0 and need no area computation.
- Inner fit keeps the outer patch's v-construction for `z_of_u` patches: keeps exact z-splits (I3)
  and the T9 writer needs no trimming curves.
- One S2 set per distinct t, fitted over the whole slab of its modules and cut at `z_ballast`:
  equal t on both sides of a joint shares the edge; `z_ballast` moves without refitting.
- `results.stage3` for both modes with one schema list, owned by T5: one struct, one owner.
- F9 and F10 shared by both modes: the same physics and the same check (§5 item 11, OD11); the
  thin-shell process stays its own (rule 7).
- Combined steel file named `<hull>_STEEL_all.step`, parallel to `<hull>_UHPC_all.step`.

## 9. Open for the owner

1. Hulls with patches whose z depends on both parameters, or entity types without an exact NURBS
   form (not C1): exact z-splits then need trimming curves in the STEP writer, or body faces fitted
   with M1–M3 (outer faces then not exact). Until decided the kernel stops with `ZNotOneParameter`
   or `NotExact`.
2. Values of the fit pass cap (`opts.max_passes`) and the Gauss order (`opts.n_gauss`): limits,
   not gates; T2 and T3 propose them with measured pass counts and convergence tables at the T3
   checkpoint.
