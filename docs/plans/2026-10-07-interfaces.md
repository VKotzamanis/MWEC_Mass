# Interface contract for the exact-geometry, Stage-3 and STEP tasks

> Binding for every task of `2026-10-06-exact-geometry-stage3-step.md` from the lanes on; read
> `AGENTS.md` first. Built on the merged T1 (`+solid/outer_rows.m`, `surface_normals.m`; used by
> T2, T3 tests) and T9 (`+step/write_step.m`, `validate_brep.m`, `tests/step_check.py`; used by T3
> tests, T10). The contract wins where a plan task section differs (§3). Cited S1–S8, F1–F14, I1–I9.

## 0. Conventions

- **Units:** SI, lengths in metres, everywhere (structs, tables, STEP). Densities kg/m³.
- **Body frame:** the `.ms2` coordinates (C1: keel z = −3.25, top z = 1.10). **World frame:** the
  same x, y; z_world = z_body + vs, vs = `vertical_shift` = x(1); z_world = 0 is the still-water
  line, so the waterline is z_body = −vs and draft = −(z_min + vs). Every field name or comment
  states its frame; geometry structs are body frame, `CG_total`, `CB`, `KM` are world frame.
- **Names (rule 10):** regions are `uhpc`, `air` (modular precast) and `ballast`, `shell`, `air`
  (thin shell); densities `rho_uhpc`, `rho_air`, `rho_ballast`, `rho_shell`; ballast top
  `z_ballast`; never `z_fill`, `rho_fill`, `rho_void`, `rho_steel` or `V_steel` for UHPC.
- **Offset distance:** a shell of design thickness t is built at d = t + eps_fit/2, eps_fit =
  0.01·t_min (AGENTS §5 item 9.2); every closed form, oracle and stand-in below uses d.
- **Errors:** `mwecmass:solid|realise|step:<Name>`. **Tests:** `tests/<area>/test_<name>.m`;
  exact-by-construction assertions only, approximation errors printed (rule 5).

## 1. Data structures

**S1 outer patch** (`geo.outer(k)`, one per entry of `model.visible_surfs`, same order; C1: 8):
`name`, `source` (ultimate source entity), `type`, `flips` (`cache.mirrors(m).effective_flips`);
`surf`: T9's `'bspline'` surface as is (`degree [du dv]`, `ctrl [nu x nv x 3]` (m, body), `knots
{ku, kv}` clamped and **not renormalised** after a split (a piece keeps its parent's parameter),
`weights` `[]` or `[nu x nv]`); `outward`: S_u × S_v points out of the hull (mirrors reverse it);
`exact`: `surf` is the `.ms2` entity converted exactly; `z_of_u`: every control row has one z
(bitwise), weights separable (§5 I3); `u_range [u0 u1]`; `z_range [z(u0) z(u1)]` (ordered by u,
not by z); `offset_kind` (`rev_z` | `ruled_parallel` | `''`; F1, used by F2); `pole [1x2]`: the row
at u0 / u1 has all control points bitwise equal; `c0_u`, `c0_v`: interior knot values of full
multiplicity (du, dv), set by F1 (F2 for S2); `seam_u0`, `seam_u1`, `seam_v0`, `seam_v1`: `[patch
boundary]` of the neighbour on that boundary, empty for a pole; **boundaries are numbered as the
parser's EdgeSnake** (`MS2Parser.eval_edge_snake`): 1 = v0, 2 = u1, 3 = v1, 4 = u0. Bitwise
equality: element-wise IEEE `==` (`isequal`), so −0 = +0. C1 values: §3, C1 oracles.

**S1b hull geometry** `geo`: `hull_name` (deck stem), `outer` (S1 array), `z_range [z_min z_max]`
(body), `analytic` (`[]` for a real deck; stand-ins only, §3). Once per deck; `config.hull_solid`
(T4b, inside T0d's cached block).

**S2 inner set** (`inner`, one per distinct t): `t`, `d` = t + eps_fit/2, `eps_fit` [m], `z_range`
(body: heights of the **inner** surface the set covers, the full hollow range, §7), `z_lo` (lowest
point of the inner surface), `refit` (true when refitted on fixed knots, F2 `opts.knots_from`),
`patches` (S1 array; `outward` true when the normal points **into the void**, out of the solid;
`exact` false; crease pieces are separate patches, `seam_*` filled), `report` (S2r). **S2r fit
report:** per inner patch `n_nodes`, `n_knots [u v]`, `n_passes`, `n_removed`, `n_check`,
`t_local_min`, `t_local_max` [m], `M1`, `M2`, `M3`, `M3_reason`; totals `ok`, `cap_reached`. Check
points: the F6 quadrature nodes and the midpoint of every knot span (where the integrals sample).

**S3 design** (input of F5): `mode` (`modular_precast` | `thin_shell`), `edges [N+1 x 1]` module
edges (body; `config.strip_edges`), `vs` [m], `t [N x 1]` [m] (thin shell: all equal; NaN for a
module with no void), `z_ballast` [m, body] (≤ z_min: no ballast; one value, may cross module
edges), `solid_modules` (indices solid regardless of `z_ballast`: the UHPC wall module; `[]` thin
shell). Module i layout (rule 11): full outer section solid below `z_ballast`; above it shell t(i)
+ air; `solid_modules` solid throughout; no slab at module joints. If `z_ballast` ≤ the inner
`z_lo`, the inner set is not cut: the solid (precast uhpc; thin shell ballast below `z_ballast`,
shell above) fills the section up to the inner patches, which alone close the void (C1: keel pole).

**S4 body** (output of F5)
- `design` (S3), `inner_t` (t of each S2 set used), `planes` (edges and `z_ballast` in (z_min, z_max)).
- `brep`: exactly the T9 writer struct (`vertices`, `curves`, `edges`, `surfaces`, `faces`,
  `bodies`; `uncertainty` unset = 1e-7 m). Extra face fields, which `write_step` and `validate_brep`
  ignore: `role` (`outer` | `inner` | `cap` | `ballast_top` | `joint_step`), `module`, `inside`,
  `outside` (region against / along the face normal = surface normal × (−1 if not `same_sense`);
  `exterior` outside the hull). Lateral faces (outer and inner): S1/S2 pieces cut at `planes` by
  knot insertion and at their `c0_u`, `c0_v` rows (an existing knot: no control point moves; a
  `c0_u` row of a `z_of_u` patch has one z, so its edge is shared under I2). Reason: OpenCASCADE
  `getMass` is wrong on a face with an interior full-multiplicity knot row (T9 `step_check.py`: up
  to +2.2 %; split faces −2.5e-10). Caps: `plane` faces inside (z_min, z_max), annuli with a
  hole loop. A collapsed row is omitted from the loop (three edges, or two when both rows collapse).
  The hull's own bottom and top are outer patches (C1: poles; box: horizontal faces), never caps.
- Modular precast faces: `cap` at every module joint (the full section where both sides are solid,
  else the annulus between the outer section and the inner section of smaller t, plus a uhpc/air
  disk where only one side is hollow, e.g. the void top under the solid wall module); `joint_step`
  (annulus between the two inner sections where t differs, uhpc/air); `ballast_top` is only the
  void-bottom disk (uhpc/air) at `z_ballast`, also when `z_ballast` is at a joint. No face where
  uhpc meets uhpc inside a module, so none at `z_ballast` ≤ `z_lo` (the inner patches close the void).
- Thin shell faces: no joint caps (lateral faces are only split at the module planes);
  `ballast_top` is the full section at `z_ballast` (ballast/shell), split into an annulus
  (ballast/shell) and a disk (ballast/air) when the inner set is cut.
- `brep.bodies`: modular precast `<hull>_UHPC_module_<i>` (solid, one closed shell each; a hollow
  module is a ring); thin shell `<hull>_STEEL_ballast` (solid, full outer section below
  `z_ballast`; absent when `z_ballast` ≤ z_min) and `<hull>_STEEL_shell` (sheet: outer faces
  above `z_ballast`).
- `shells` (signed face lists for F11 and tests): precast `all_outer`, `all_voids{c}` (the fused
  `<hull>_UHPC_all` solid with its cavities); thin shell `layer`, `void` (closed, mass regions only).

**S5 section loop** (F4, F6b): `z` [m, body], `pieces(k)`: `patch`, `u` (iso-u), `curve` (T9 curve
struct, the exact iso-u row in v), `dir` (±1); `pts [n x 3]` CCW on the exact curves (no repeated
point); `area`, `centroid [1x2]`, `I [Ixx Iyy Ixy]` about the origin [m⁴] by Green's theorem (Gauss
on knot spans); `simple`. Body section (F6b): `z`, `module`, `outer`, `inner` (S5 or empty), `solid`.

**S6 body properties** (F6): for each region r and module i: `V(i)` [m³], `S(i,:)` = ∫x dV [m⁴],
`J(:,:,i)` = ∫x xᵀ dV [m⁵], body frame about the origin; `regions.<r>` holds these; `modules(i)`:
`V`, `V_<r>`, `mass`, `rho_eff` = mass/V, `CG_body [1x3]`; `total`: `mass`, `CG_body [1x3]`,
`I_origin [3x3]`, `I_cg [3x3]` (I = ∫ρ(|x|²E − x xᵀ)dV); `quad` (Gauss order used).

**S7 hydrostatics** (F7): `vs`, `draft`, `submersion` (`partial` | `full` | `none`), `V_sub` [m³],
`CB_body`, `CB` (world) [1x3], `S_wet` [m²] (∫|S_u × S_v| over the outer pieces below the
waterline, same Gauss rule), `Aw` [m²], `xF`, `yF` (centre of flotation), `I_wp_xx`, `I_wp_yy`
(about the CF; fixes I12), `I_wp_yy_origin`, `BM_L` = I_wp_yy/V_sub, `KM` (world z) = CB(3) + BM_L,
`waterline` (S5 at z_body = −vs). Waterline at or above z_max: `full`, V_sub = V_hull, CB = hull
centroid, `Aw`, `I_wp_*`, `BM_L` = 0, KM = CB(3), `waterline` empty. At or below z_min: `none`,
V_sub = `S_wet` = `Aw` = 0, CB, `xF`, `yF`, `BM_L`, KM = NaN, `waterline` empty.

**S8 realised design** (`results.stage3`, both modes; `[]` for `preliminary`; replaces
`results.constructability` and `results.steel_data`)
- `mode`, `hull_name`, `status` (`accepted` | `failed`), `reason` ('' when accepted), `escalation`
  (`split` | `fixed_draft` | `spill` | `draft_free`: the last step run), `vs`, `draft` [m] (world).
- `stage2`: `vs`, `rho [N x 1]`, `mass`, `Z_CG` (= `Final3D.CG_total(3)`, world), `GM`
  (`Final3D.GM_L`), `T_heave`, `T_pitch` (`Final3D.periods.heave/.pitch`, coupled).
- `rho` (densities used, by region), `design` (S3), `k_star` (precast: ballast module),
  `V_uhpc_target [N x 1]` (precast split, AGENTS §3 item 4.1). `modules(i)`: `z_lo`, `z_hi` [m,
  body], `t` [m] (NaN: no void), `h_ballast` [m] (from the module bottom, OD6 ii), `V`,
  `V_<region>`, `mass`, `rho_eff`, `rho_stage2`, `rho_floor`, `CG_world [1x3]`.
- `props`: F9 output; `check`: F10 output; `solver`: per escalation step run `exitflag`,
  `iterations`, `fval`, `max_eq_violation`; `fit`: S2r of every inner set; `body`: S4; `step_files`:
  `name`, `path`, `bodies` (F11). `final_props` = `props` + `stage3_status`, `stage3_check`.

## 2. Functions

All take and return the structs of §1; none reads globals or writes files except F11 and F13.

| | Signature | Contract |
|---|---|---|
| F1 | `geo = mwecmass.solid.outer_nurbs(model)` | S1b from the `.ms2` entity tree, every entity converted exactly. Points: FramePoint, MirrPoint, AbsBead (on the converted curve, at the parser's parameter map). Curves: BCurve as is, Line degree 1, Arc rational quadratic, BSubCurve by knot insertion, PolyCurve2 joined with C0 knots after exact degree elevation of its pieces to the highest degree, ProjCurve by projecting control points (kept coordinates copied bitwise), EdgeSnake as the boundary row or column of its converted parent surface, taken from the same array (bitwise seam; edge numbers are the parser's, S1: 1 = v0, 2 = u1, 3 = v1, 4 = u0). Surfaces: RevSurf as profile × rational arc (multiples of 90° use exact 0, ±1), RuledSurf as degree 1 in v between its two curves after one common reparameterisation (one knot vector, equal weights), mirrors by flipping control points. `offset_kind` = `rev_z` when the axis end points have equal x and y within T1's 16-ulp bound; the revolution then uses the (x, y) of the axis Line's second end point, every profile control point that is an axis end point by entity (the profile ends at a point or bead that defines the axis) takes that (x, y) bitwise, and control points are axis point + radial vector (a zero radius gives the axis point bitwise); `ruled_parallel` when the cross product of every ruling with the first is bitwise zero. Any other entity: error `NotExact` (§9 item 1). C1 results: §3, C1 oracles. |
| F2 | `[inner, rep] = mwecmass.solid.offset_surface(model, cache, geo, t, z_range, opts)` | S2 at thickness t over `z_range` (body). Nodes from `MS2Parser` + T1 `surface_normals`, offset by d, folds trimmed (`trim_fold.m`), faces split at creases into pieces that meet along u-boundaries, cubic fit (`fit_bspline_surface.m`), knot insertion where M1–M3 fail, knot removal while they hold (AGENTS §5 item 9). Keeps the outer v-construction for `rev_z` (offset profile fitted in u, revolved with the same rational arc) and `ruled_parallel` (both boundary curves offset at the same u nodes, one knot vector, ruled in v), as I3 and §8; `''`: error `OffsetNotStructured`. Inner patches of a `z_of_u` patch are `z_of_u` with monotone z. Source points are taken wherever needed so the inner patches cover `z_range`; past an open end (precast: the bottom of the wall module) F5 cuts them with F3b. v-seams must be tangent-continuous; a seam node is offset once, so neighbours share boundary curves bitwise (I2); a C0 v-seam is out of scope (§9 item 1) and fails loudly through M1–M3. `opts.t_min` (eps_fit), `opts.max_passes` (§9 item 2); cap reached: error `FitNotConverged` listing the failing patches and metrics. `opts.knots_from` (an S2 set): same pieces and knot vectors, only the control points refitted at the new t (§7 item 4); M1–M3 reported in S2r, not refined. Tests: identity at the same t; properties smooth in t. |
| F2b | `d_close = mwecmass.solid.void_closing_distance(model, cache, geo, z_range)` | Smallest offset distance at which offset layers from opposite sides meet inside `z_range` (not a fold of one layer). C1 neck: 0.10 m. The design bound is t_max = d_close − eps_fit/2: at t = t_max the void closes (d = d_close); an evaluation with d ≥ d_close errors `VoidClosed` (§8). |
| F3 | `[S, Su, Sv] = mwecmass.solid.eval_bspline_surface(surf, u, v)`; `[C, Cs] = …eval_bspline_curve(curve, s)` | Rational or not; u, v column vectors; [n x 3]. |
| F3b | `[lo, hi, u_star] = mwecmass.solid.split_bspline_surface(patch, z)` | `z_of_u` patch, z strictly between z(u0) and z(u1), in either order: u* from z(u*) = z, knot insertion to multiplicity du; `lo`, `hi` keep the parent parameter and share the cut row bitwise. Errors: `ZNotOneParameter` (not `z_of_u`), `ZOutside`, `ZNotMonotonic` (z(u) = z has more than one root). |
| F4 | `loop = mwecmass.solid.slice_bspline_surface(patches, z)` | S5 from the iso-u rows; one closed simple loop or error `SectionNotClosed` (as T1). Constant-z patches (z_range(1) = z_range(2)) are skipped. |
| F5 | `body = mwecmass.solid.build_body(geo, design, inner)` | S4. `inner` holds one S2 set per distinct t whose `z_range` covers its modules (§7). Builds every face once; edges shared by index. Layout of S3, including `z_ballast` ≤ `z_lo`; faces per mode as in S4. A patch is cut only at planes strictly inside its z range (constant-z patches never). At every `joint_step` the loop of the larger t must lie inside the loop of the smaller t, tested on the exact loops without tolerance (no intersection of the two loops, and one point of the larger-t loop inside the smaller-t loop); otherwise error `JointNotNested`, which Stage 3 reports as a failed evaluation. |
| F6 | `bp = mwecmass.solid.body_properties(body, rho, opts)` | S6 by the divergence theorem with fields (f,0,0): V = ∮x n_x, ∫x = ∮x²/2 n_x, ∫y = ∮xy n_x, ∫z = ∮xz n_x, ∫x² = ∮x³/3 n_x, ∫y² = ∮xy² n_x, ∫z² = ∮xz² n_x, ∫xy = ∮x²y/2 n_x, ∫xz = ∮x²z/2 n_x, ∫yz = ∮xyz n_x. Horizontal faces (caps and constant-z outer patches) contribute exactly 0 (n_x = 0), so only lateral faces are integrated; a region's integral sums faces with `inside` = r minus faces with `outside` = r. Gauss–Legendre of order `opts.n_gauss` on every knot span (exact for polynomial faces, convergent for rational ones). |
| F6b | `sec = mwecmass.solid.body_section(body, z)` | Body section (S5 block). |
| F7 | `hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, opts)` | S7 on the exact outer patches cut at z_body = −vs (F3b, F6 on the outer pieces, F4 for the waterplane); `full` and `none` as in S7, without a cut. |
| F8 | `fl = mwecmass.driver.density_floors(model, cache, geo, edges, t_min, rho_solid, rho_air, solid_modules)` | Per module ρ_min = [ρ_solid·V_solid + ρ_air·V_air]/V with a t_min shell, no ballast (F2, F5, F6); solid modules give ρ_solid. Returns `rho_min`, `V`, `V_solid`, `V_air` [N x 1], `fit` (S2r). |
| F9 | `props = mwecmass.realise.evaluate_realised(bp, hs, design, config)` | `final_props` fields of `export_schema.m` from the exact body. New: `vertical_shift` = `design.vs`, `mass_total`, `CG_total` = [0 0 CG_body(3)+vs] (declared symmetric model; x, y kept in `bp`), `Inertia_Tensor` = `I_cg`, `Ixx/Iyy/Izz`, S7 fields, `A_sub` = `hs.S_wet`, `GM_L` = KM − CG_total(3), `periods.heave_uncoupled` = 2π√((M+A33)/K33) and `periods.pitch_uncoupled` = 2π√((I_cg,yy+A55)/K55), Inf when K ≤ 1e-6 (`properties_3d.m` §11; today from `steel_data`), `densities_at_nodes` = `realised_strip_density` = `[bp.modules.rho_eff]'` (layout of x(2:end)), `cross_section` = `config.profile` (as `optim/run.m`), `realised_strip_edges` = `design.edges`, `components(i)` (`density` = `rho_eff`, `z_level` = module mid-height + vs, world, as Final3D), `density_profile_source` = `'realised_partition'`. As `build_realised_properties.m` today: `K_hydro` (K33 = ρ_w g A_w; K55 = M g GM if GM > 0, else 0, `properties_3d.m` §8), A = `mwecmass.bem.interpolate_at_draft(vs, config, CG_total(3))`, `periods.heave`, `.pitch`, `.surge` (Inf, K11 = 0), `coupled_periods`, `participation_factors`, `coupled_modes`, `surge_per_pitch` (`coupled_periods_by_share`), `mass_discrepancy`, `mass_buoyant_force` = ρ_w V_sub, `K_pto` = zeros(3), `K_total` = `K_hydro`, `MassMatrix_CG`/`_Origin` (`build_mass_matrix`), `fill_method` (renamed only if T0b renames it). `realised_strips` leaves the schema (T5): only the realise files T5 and T7 rewrite read it; its data is S8 `modules`. |
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

Where the plan's task sections differ, this contract wins:
- Paths: F10 is `+realise/check_against_stage2.m` (both modes, §5 item 11); F9 is the new
  `+realise/evaluate_realised.m`, replacing `build_realised_properties.m`; T5 and T7 write the F14
  calls of F11 and F13, so T8 and T10 do not edit `run.m`.
- Methods: F6 by the divergence theorem on B-rep faces and F7 on the outer patches cut at the
  waterline, not quadrature of sections in z (plan T3); sections (F4, F6b) serve figures and tests.
- Bounds: t_max,i (T6) and `t_max` (T7) are d_close − eps_fit/2 (F2b). Baseline: "updated in the
  same commit" (T4a, T4b) is superseded by §6 (regenerated once, after J2). Switches: one flag,
  `out.save.stage3.step` (F14, §4), not the `out.save.stage3.step_*` switches of the T10 Files line.

**Stand-in kit SK** (one graded step at the start of lane K, merged before any consumer starts):
- `tests/standins/fixtures/` (permanent): `cylinder.ms2` (quarter RevSurf about the z axis of a
  PolyCurve2 of three Lines, bottom disk, side, top disk, mirrored in x and y like C1; sharp rims),
  `box.ms2` (side faces RuledSurf between two vertical Lines: horizontal parallel rulings; bottom
  and top faces RuledSurf between two horizontal Lines of constant z, which F3b never cuts, F4 skips
  and which add 0 to F6: the outer patches close the hull, and F5 keeps one rule, caps only
  strictly inside the z range); `sti_closed_form.m` (V, S, J, sections, hydrostatics for any S3
  design on these, at d: inner cylinder of radius R − d between z0 + d and z1 − d, inner box shifted
  inward by d; `d_close` = min(R, H/2), or half the smallest box dimension); `sti_inner_box.m` (the
  exact S2 set of the box: its six planar ruled faces shifted inward by d, so the void is closed);
  `sti_config.m` (`config` fields `RHO_WATER`, `G`, the one-draft hydro table `hydro_drafts`,
  `hydro_z_cg`, `added_mass_diagonal`, `added_mass_full`, `radiation_damping_full`, `profile`,
  `strip_edges`, `mass_acceptable_pct`; tests put the stand-in S1b into `config.hull_solid`);
  `sti_stage2.m` (Final3D-like, consistent with the closed form); `sti_realised.m` (S8 of a design).
- `tests/standins/+mwecmass/+solid/` and `+realise/`: closed-form versions of F1–F7 including F6b,
  F9, F10 that set or require `geo.analytic` (error `mwecmass:standin:NotAnalytic` otherwise).
- Consumer tests call `addpath(fullfile(root,'tests','standins'),'-end')`: a real function in `src/`
  wins (Octave 8.4 and MATLAB merge package folders, earlier path entries first; verified in
  Octave). The join that merges a producer deletes its stand-in functions; the consumer tests then
  run unchanged on the real code (the join test).
- Fixtures per join test: the cylinder goes through the real F1, F2, F5–F7 (its rims are u-creases
  of convex radius 0 < t, so it also tests fold trimming and crease splitting, plan T2; rational
  faces, so results are compared with the closed form and printed); the join tests of T4b, T5, T6
  and T7 use it, since their code calls F2. The box serves F5–F7, F11, F12 with `sti_inner_box`
  (its vertical corners are C0 v-seams, outside F2's scope); results polynomial, asserted exact.
- C1 oracles (by construction): `BeadBottom` evaluates to x = 3.9e-16 m, `BeadTop` exactly to
  (0, 1); both ends of `curve6` lie at these beads, so `surface1` (`rev_z`, `z_range` [1.1 −3.25])
  has poles at its keel and neck-top rows. `Edge_For_Dev` = EdgeSnake edge 3 of `surface1` (v1,
  φ = 360°, the profile in y = 1); `curve7`, its projection on y = 0, copies x and z, so `surface2`
  is `ruled_parallel`, and its ridge row u0 (z = 1.1) and keel row u1 have x = 0, shared with its
  X-mirror's. `c0_u` of both: the arc-to-`curve1` joint, z = 1.0; `c0_v` empty; all 8 `exact`;
  v-seams tangent-continuous (`surface1`–`surface2` by construction, mirror seams in symmetry planes).
  For z ∈ [−0.5, 1.0] the section is a stadium of radius 0.10 about (0, ±1), flat sides x = ±0.10
  (curve1's first span has collinear control points), area 0.4 + 0.01π m², inner radius and
  half-width 0.10 − d; z ∈ [1.0, 1.1] is a half cylinder (r 0.10, length 2) with two quarter
  spheres, V = 0.01π + (2/3)π·0.001 m³; neck `d_close` = 0.10 m. Smallest convex principal radius
  0.100 m (top arc, neck hoop), so the offset folds only from t = 0.100 m on (AGENTS rule 3).
- **SK acceptance** (reproduced by the grader): every stand-in S4 (split at `c0_u` rows) passes
  `validate_brep`; the box and cylinder bodies written by `write_step` pass `step_check.py` (closed
  solids, solid count), `occ_volumes` compared with `sti_closed_form` (printed); V, S, J of the
  closed form equal an independent Gauss–Legendre integration in z of the closed-form sections
  (piecewise polynomial, so exact; asserted to rounding); region volumes sum to V_module
  (asserted); every stand-in errors `NotAnalytic` on the C1 `geo`.

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
| `output/dispatch.m` | T4a (only if the `Monotonic density` line goes with `c_mono`) → T8 (F13 call) → T10 (`validate_save_flags`: `out.save.stage3.step`) |
| `tools/baseline_run.m` | T0 → T5 (summary reads `results.stage3`) |
| `validation/diagnostics/uhpc_mass_balance.m` | T0b → T11 |
| `tests/run_tests.m`, `.gitignore` | T0 → T0d |
| `tests/baseline/*.json` | T0 → J2 run step |
| `tests/standins/*` | SK → J1 (deletes the stand-ins of F1–F7, F6b) → J2 (deletes the F9, F10 stand-ins); consumers do not edit it |
| `MS2Parser.m` | T0c only. `+solid/*`: T2 and T3 own disjoint files. `+output/+step/*`: T9 → T10 (new files; `write_step.m` changes only if a C1 body fails `step_check`, by T3). |

## 5. Invariants every implementation tests

- **I1 Volume closure (rule 11)**, to a stated multiple of machine epsilon justified in the test:
  per module V_uhpc + V_air = V_module (precast), V_ballast + V_shell + V_air = V_module (thin
  shell); Σ_i V_module,i = V_hull, V_hull integrated on the unsplit outer patches with the cut
  parameters u* as extra breakpoints (same (u, v) samples). With `z_ballast` inside a module, at a
  module edge, spilled, and below the inner `z_lo`.
- **I2 Shared edges:** every closed shell passes `validate_brep`; the boundary rows of the two faces
  at an edge equal its curve bitwise (`isequal`: control points, knots, weights), on u- and v-seams,
  including C1's keel seam of `surface2` and its X-mirror (§3, C1 oracles).
- **I3 Exact z-splits** where z depends on u only. RevSurf about a vertical axis: rotation keeps z,
  so each control row has one z and w_ij = a_i b_j, hence z(u,v) = Σ N_i a_i z_i / Σ N_i a_i. C1
  RuledSurf: `curve7` keeps z at both ends of every ruling. A plane z = z_k meets the patch in the
  iso-u row u*; knot insertion cuts it without moving the surface (test: parent and pieces agree at
  the same (u, v) to rounding; cap vertices within T9's 1e-7 m). Offsets keep it for `rev_z` and
  `ruled_parallel` only (normal in the meridian plane; normal constant along a ruling).
- **I4 Units:** METRE everywhere; `step_check` reports `declared_length_unit = METRE`. **I5 Names:**
  new code uses only the names of §0 (`grep` for the old names in new files is empty).
- **I6 Frames:** CG_total(3) = CG_body(3) + vs; KM, CB world. **I7 Orientation:** outer normals
  out of the hull, inner normals into the void; every region volume positive.
- **I8 C1 symmetry:** CG x = y = 0 and products of inertia 0 to machine precision.
- **I9 Exactness of F1:** F4 sections of the NURBS and T1 `outer_rows` at the same z: the
  difference is printed per height and asserted within T1's documented in-plane accuracy at that
  height (about sqrt(16 ulp/|z_ss|), 5.5e-8 m at the C1 shoulder z = −1), not to rounding.

## 6. Lanes and joins

```
J0   merge T0, T1, T9, this contract
K    SK ─► T2 ─► T3 ─────────────────────► J1 (T2+T3)  ── owner checkpoint after T3
N    T0b ─► merged at once (all later lanes start from it)        T0c (MS2Parser.m) in parallel
P    [T0b, SK, T0c merged] ─► T0d ─► T4a ─► T4b (stand-ins; merges after J1)
U    [T0b, SK merged] ─► T5 ─► T6
S    [T0b, SK merged] ─► T7
O    [T0b, SK merged] ─► T8 ∥ T10
J2   after J1 and lane P: the group T5, T6, T7, T8, T10 merges as one chain in that order
     ─► first whole-pipeline runs ─► owner checkpoints after T6 and after T10
G    after J2: T11 ─► T12 ─► T13
```
- Tasks are graded on component tests (stand-ins, C1 pieces, real kernel after J1). Whole-pipeline
  items are checked after J2 (owner): the T0 baseline on the commit after T0d (identity of T0b, T0c,
  T0d); C1 runs of both modes on the J2 head for T4a, T4b, T5, T6, T7, T10; then the baseline is
  regenerated once (`tools/write_octave_baseline.m`), changes printed per task. Failures are fixed
  on a fix branch; inside J2 the integration branch is pushed after the last merge and suite.

## 7. Stage-3 evaluation per solver iteration, and what must be cached

Per evaluation of (t_k*, t_i, `z_ballast`, vs) or (t, `z_ballast`, vs): S2 sets per distinct t (F2,
the expensive step), F5 (knot insertion at `z_ballast`), F6 on the changed faces, F7 (vs only), F9
(BEM interpolation, 3 × 3 eigenproblem), F10. Needs:
1. S1b, the outer pieces cut at the module edges and their F6 integrals: once per deck and edges.
2. S2 sets keyed by the bitwise value of t only, each fitted once over the full hollow range
   (precast: z_min to the bottom of the solid wall module; thin shell: the whole hull): one set
   serves every module with that t (equal t shares joint edges; `z_ballast` moves without
   refitting), and a finite-difference step on one t never refits another.
3. F7 keyed by vs: one call while the draft is fixed (`split`, `fixed_draft`, `spill`).
4. Smoothness: properties differentiable in t and `z_ballast` within one solve, so each S2 set keeps
   its knot vectors during a solve (F2 `opts.knots_from`); the adaptive fit runs on the final t, and
   a changed knot structure restarts the solve there. `z_ballast` enters only by knot insertion.

## 8. Decisions taken here (simplest option that meets the rules)

- One B-rep (S4) for mass, sections, figures and STEP: one geometry (§5 item 1); caps add 0 to F6.
- Inner fit keeps the outer v-construction only where it is the exact normal offset (I3): exact
  z-splits, no trimming curves in the T9 writer. `rev_z` within T1's 16-ulp bound and the F1 axis
  snapping: no new tolerance, and C1's poles and keel seam come out bitwise.
- One S2 set per distinct t (§7 item 2). `A_sub` = F7 `S_wet` on the exact outer pieces (one kernel).
- No margin below t_max (agent's decision, owner informed): solvers use t_max = d_close − eps_fit/2;
  an evaluation with d ≥ d_close returns error `VoidClosed` and counts as a failed evaluation.
- `results.stage3`, F9, F10 shared by both modes, schema owned by T5 (one struct, one owner;
  OD11); the thin-shell process stays its own (rule 7). `<hull>_STEEL_all.step` mirrors `_UHPC_all`.

## 9. Open for the owner, and decisions due

1. Open: hulls (not C1) with z depending on both parameters or not monotone in u, `offset_kind`
   `''`, C0 v-seams or entities without an exact NURBS form need trimming curves in the STEP writer,
   or faces fitted in both directions with M1–M3 (outer faces then not exact). Until decided the
   kernel stops with `ZNotOneParameter`, `ZNotMonotonic`, `OffsetNotStructured` or `NotExact`.
2. Due at the T3 checkpoint: `opts.max_passes` and `opts.n_gauss` (limits, not gates), set from the
   pass counts and convergence tables T2 and T3 measure and present there.
