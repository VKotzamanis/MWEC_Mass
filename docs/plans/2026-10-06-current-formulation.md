# MWEC_Mass: what the code optimises today (Stages 1–3)

Read-only investigation. No repository file was changed. Scratch scripts are in `/tmp/claude-0/-home-user-MWEC-Mass/4eed5a18-bc85-5f1a-b911-b09cefd00ef6/scratchpad/formulation/` (`peek.py`, `dump.py`, `numbers.py`, `s3.py`).

**Path convention:** paths starting with `+` are relative to `/home/user/MWEC_Mass/src/+mwecmass/`. `WEC_User_Input.m` and `Output/...` are relative to the repository root.

**Basis:** every formulation statement below comes from code I opened or grepped. I read parts of docs/ earlier in the session, before your correction. Nothing below depends on them. Docs appear only in the last section, marked "docs claim" against "code fact".

## 0. Source files read

| File | Lines read |
|---|---|
| WEC_User_Input.m | 1-156 |
| +driver/run.m | 1-113 |
| +driver/build_config.m | 1-1234 |
| +driver/build_strip_geometry_tables.m | 1-213 |
| +optim/run.m | 1-446 |
| +optim/stage1_sweep.m | 1-148 |
| +optim/stage1_screen_draft.m | 1-53 |
| +optim/stage1_oneshot.m | 1-66 |
| +optim/stage1_trained.m | 1-301 |
| +optim/solve_2d_surrogate.m | 1-478 |
| +optim/PID_Controller.m | 1-135 |
| +optim/check_3d_convergence.m | 1-65 |
| +optim/report_assemble_results.m | 1-67 |
| +optim/hams_enrichment_action.m | 1-11 |
| +optim/range_penalty.m | 1-15 |
| +optim/stage2_bounds.m | 1-22 |
| +optim/stage2_constraints.m | 1-67 |
| +optim/stage2_objective.m | 1-35 |
| +hydrostatics/properties_2d.m | 1-458 |
| +hydrostatics/properties_3d.m | 1-317 |
| +hydrostatics/coupled_periods_by_share.m | 1-86 |
| +hydrostatics/compute_strip.m | grep only: 5-101, 237-252 |
| +bem/interpolate_at_draft.m | 1-111 |
| +bem/transform_to_cg.m | 1-55 |
| +bem/retransform_at_cg.m | 1-36 |
| +bem/interpolate_matrix.m | 1-32 |
| +bem/rebuild_config_hydro.m | 1-87 |
| +bem/load_hydro_cache.m | 1-78 |
| +bem/+hams_mrel/run_at_draft.m | grep only: 1, 29-30, 86, 110, 128, 159 |
| +realise/+preliminary/run.m | 1-8 |
| +realise/+thin_shell/run.m | 1-63 |
| +realise/+thin_shell/solve.m | 1-524 |
| +realise/+thin_shell/evaluate_design_point.m | 1-219 |
| +realise/+thin_shell/build_geometry_grid.m | 1-99 |
| +realise/+thin_shell/integrate_split.m | 1-99 |
| +realise/+thin_shell/inner_properties_at_z.m | 1-64 |
| +realise/+thin_shell/strip_partition_volumes.m | 1-20 |
| +realise/+modular_precast/run.m | 1-38 |
| +realise/+modular_precast/build_realised_properties.m | 1-15 |
| +realise/+modular_precast/solve_and_extract.m | 1-138 |
| +realise/+modular_precast/solve.m | 1-478 |
| +realise/+modular_precast/evaluate_design_point.m | 1-140 |
| +realise/+modular_precast/build_geometry_grid.m | 1-115 |
| +realise/+modular_precast/extract_strip_geometry.m | grep only: 1-8, 64-74, 94, 309-430 |
| +realise/build_realised_properties.m | 1-311 |
| +output/export_results.m | grep only: 1-56 |
| Output/thin_shell/steel_solve.log | line 17 |

## 1. Pipeline and the modes the code can run

- **Order of calls.**
  - `+driver/run.m:72` builds `config = build_config(...)`.
  - `+driver/run.m:79` runs `[opt_results, x_opt, iteration] = mwecmass.optim.run(config)` (Stages 1 and 2).
  - `+driver/run.m:82-83` runs `realise_fn = str2func(['mwecmass.realise.' in.materials.realisation_type '.run'])` (Stage 3).
- **Stage-1 modes.** The switch in `+optim/run.m:99-132` has exactly four cases:
  - `'sweep'` (line 101)
  - `'skip'` (line 110)
  - `'oneshot'` (line 120)
  - `'trained'` (line 125)
  - anything else errors at line 131.
- **`screen_draft` is not a mode.** It is Tier 1 inside `sweep` (`+optim/stage1_sweep.m:42`).
- **Shipped default.** `WEC_User_Input.m:74` sets `in.pid.stage1_mode = 'sweep'`. All three C1 files have `config.stage1_mode = 'sweep'`.
- **Realisation flags** (`+driver/build_config.m:35-45`):
  - `thin_shell` gives `enable_steel_solve = true`.
  - `modular_precast` gives `enable_constructability = true`.
  - `preliminary` sets both false.

Every optimisation problem in the code:

| # | Problem | Solver call |
|---|---|---|
| 1a | Stage 1 `sweep`, Tier 1 | single evaluations, no optimiser |
| 1b | Stage 1 `sweep`, Tier 2 | fmincon at fixed vs, up to K draft nodes |
| 1c | Stage 1 `skip` | no optimisation |
| 1d | Stage 1 `oneshot` | one 2-D surrogate fmincon |
| 1e | Stage 1 `trained` | PID fixed-point loop around the 2-D surrogate fmincon, plus a pre-calibration |
| 1f | 2-D surrogate problem, shared by 1d and 1e | `+optim/solve_2d_surrogate.m` |
| 2A | Stage 2 Phase A, conditional | fmincon with vs frozen |
| 2 | Stage 2 main | fmincon |
| 3p | Stage 3 `preliminary` | no optimisation |
| 3t | Stage 3 `thin_shell` | fmincon |
| 3m | Stage 3 `modular_precast` | fmincon |

**range_penalty**, used by every objective below (`+optim/range_penalty.m:6-14`):

```
if r < -1:  delta = -r - 1;  phi = 1 + 2*delta + k_amp*delta^2
elseif r > 1: delta = r - 1; phi = 1 + 2*delta + k_amp*delta^2
else: phi = r^2
```

The 2-D path uses an identical local copy, `range_penalty_2d` (`+optim/solve_2d_surrogate.m:465-478`).

**Common inputs:**
- `k_amp = config.zone_k_amp = in.objective.zone_k_amp = 5.0` (`WEC_User_Input.m:102`; `+driver/build_config.m:1010`).
- `config.penalty_guard = in.objective.penalty_guard = 1e4` (`WEC_User_Input.m:103`; `build_config.m:1011`).

## 2. Shared physics models

### 2.1 properties_3d (used by Stage 1 sweep, Stage 2 and Phase A)

**Hydrostatics: table interpolation in the body frame.**
- `z_wl = -props.vertical_shift` (`+hydrostatics/properties_3d.m:24`).
- `z_sub_top = min(z_wl, hull_z_max)` (line 37).
- `Aw`, `I_wp_yy`, `V_sub`, `CB_z_body` come from `interp1(config.Aw_table_z, config.<table>, z_sub_top, ...)` (lines 41-45).
- `props.CB = [0,0,CB_z_body + vertical_shift]` (line 46).
- If tables are absent it falls back to `compute_submerged` (lines 58-83). Tables are always present in the driver path (`build_config.m:665-672`).

**Stability:**
- `BM_L = props.I_wp_yy / props.V_sub` (line 100).
- `props.KM = props.CB(3) + BM_L` (line 106).
- `props.mass_buoyant_force = props.V_sub * config.RHO_WATER` (line 105).
- `props.GM_L = props.KM - props.CG_total(3)` (line 180).

**Mass, CG and inertia from precomputed strip integrals with uniform density per strip.**
- Density per strip: `rho_i = max(lb,min(ub,densities_at_nodes(i)))` (lines 126-127).
- Mass and moments:
  - `strip_mass = V_strip*rho_i` (line 129)
  - `cg_z_num += strip_mass*config.strip_CB_z(i)` (line 132)
  - `Iyy_total += rho_i*config.strip_Iyy(i)` (line 133)
- CG: `props.CG_total = [0,0,cg_z_body + props.vertical_shift]` (line 146).
- Iyy about the CG: `Iyy_about_cg = max(0, Iyy_total - mass_total*cg_z_body^2)` (lines 151-152).
- **Where the strip arrays come from.** `strip_V`, `strip_CB_z` and `strip_Iyy` are computed once at configuration time by `compute_strip(parser, z_lo_i, z_hi_i, strip_quad_opts)` on the Aw/I_wp tables (`+driver/build_strip_geometry_tables.m:181-201`). From the grep of `+hydrostatics/compute_strip.m:99`, `strip.Iyy = int_x2 + int_z2`, i.e. about the body origin.
- **Consequence: mass, body CG and Iyy do not depend on vs** (they depend only on rho_i).
- **Strip layout:**
  - Non-constructability: `density_nodes_z = linspace(hull_z_min, hull_z_max, N)` (`build_config.m:245-246`). Edges are the midpoints between nodes, clamped to the hull ends (`build_strip_geometry_tables.m:14-28`). So the first and last strips are half height. C1 edges: [-3.25, -2.70625, -1.61875, -0.53125, 0.55625, 1.1].
  - Constructability: N-1 platform strips plus 1 wall strip (`build_config.m:187-233`).

**Stiffness:**
- `K33_hydro = RHO_WATER*G*Aw` (line 186).
- `K55_hydro = mass_total*G*GM_L` if `GM_L > 0`, else 0 (lines 189-193).
- `K11 = 0` (line 209).

**Added mass:** `[A11,A33,A55,A_full,B_full] = interpolate_at_draft(props.vertical_shift, config, props.CG_total(3))` (lines 202-204). See 2.3.

**Periods:**
- The uncoupled values `2*pi*sqrt(M_virtual/K)` are kept as `heave_uncoupled` / `pitch_uncoupled` (lines 218-238 and 250-251).
- **The objective reads coupled periods.**
  - `M_total = diag([m, m, Iyy_about_cg]) + props.A_full` (lines 246-247).
  - Then `coupled_periods_by_share(M_total, K_total_full)` (lines 253-257).
  - That function computes `eig(M_total \ K_full)` (`+hydrostatics/coupled_periods_by_share.m:17`).
  - Each mode is named by its largest kinetic-energy share (lines 27-37 and 42-51). It errors if a DOF with positive stiffness has no mode with at least 50% share (lines 58-69).
- **Closed form I verified numerically.**
  - In the C1 data `A_full(1,2) = A_full(2,3) = 0` and K is diagonal. So coupled heave equals uncoupled heave.
  - Coupled pitch equals `2*pi*sqrt((Iyy + A55 - A15^2/(M + A11))/K55)`.
  - Reproduced to machine precision for C1 `Final3D` (3.47112 s) and for the thin-shell `final_props` (3.81328 s).

### 2.2 properties_2d: the 2-D surrogate (used only by Stage 1 oneshot/trained)

**Geometry.**
- Midplane profile shifted by vs: `shifted_profile = config.profile + [0, props.vertical_shift]` (`+hydrostatics/properties_2d.m:18`).
- Effective width `eff_w = interp1(y_span_z_levels, y_span_table, -vs)`, floored at `eff_w_floor` (lines 23-30). The y-span table is area-corrected (`build_config.m:715-785`).
- `V_sub = sub_area*eff_w*k_vol` (line 47).
- `Aw = wp_width*eff_w` (line 48).
- `CB` is the 2-D centroid of the clipped polygon in the world frame (lines 38-51 and 445-447).
- `KM = CB(3) + BM` with `BM = I_wp_yy/V_sub` (lines 54-78). Here `I_wp_yy = eff_w*w^3/12 + Aw*d_x^2`.

**Strip integration over only `num_strips = config.n_density_strips` strips.**
- `n_density_strips = num_ballast_sections = 5` (`build_config.m:715`; `properties_2d.m:84-89`).
- `strip_vol = strip_area*strip_eff_w*k_vol` (line 128).
- **V_sub is a step function of vs.** Only whole strips with centre `z_cur <= 0` count (lines 131-133), and this overwrites the continuous V_sub (lines 175-178). The C1 data confirms it: `stage1_2d.properties.V_sub = 7.882097967127017` at vs = 0.7786 (prelim) and at vs = 1.075 (modular).
- **Density is not a per-strip constant.**
  - Default: `rho_raw = interp1(config.density_nodes_z, densities_at_nodes, z_orig, 'linear','extrap')` (lines 149-150), i.e. linear interpolation between nodes.
  - Constructability: binned by `strip_edges` (lines 138-147).
- The `config.shell` branch (lines 157-167) is dead, because `config.shell = []` always (`build_config.m:283`).

**CG and GM.**
- `CG_total = [cg_x, 0, cg_z]` in the world frame (lines 170 and 188-193).
- `GM_uncorrected = KM - CG_total(3)` (lines 201-202).
- `props.GM = KM - CG_total(3)*config.k_gm` (lines 204-205).

**Inertia, stiffness and added mass.**
- `Iyy` about the 2-D CG including x terms: `strip_mass*(d(1)^2 + d(3)^2)` (lines 210-275).
- `K55 = mass_total*G*props.GM` uses the k_gm-biased GM (lines 283-287).
- Added mass is evaluated at the raw `CG_total(3)` (lines 302-304).

**Periods:** coupled, the same way as 3-D (lines 353-364).

### 2.3 Added mass: cache and CG transfer (all stages)

- **Building the cache fields in config.**
  - `build_config.m:834-851` takes `added_mass_inf{i}` ("6×6 at origin", line 835), extracts DOFs [1,3,5] and applies `transform_to_cg(A_3x3_origin, B_3x3_origin, z_cg_i)`.
  - `z_cg_i = hydro_table.z_cg(i)` ("CG_z in global frame", line 837). `config.hydro_z_cg` is set at line 826.
- **What z_cg is.** In the C1 config, `hydro_z_cg - hydro_drafts - hull_centroid(3) = 0` for all 11 nodes. So the cache reference CG is the world-frame uniform-density centroid, `centroid_z + vs`.
- **Not verified:** that the WAMIT cache matrices are about the world origin. The cache file is MAT v7.3 and was not opened. The HAMS route comment says `XR = [0,0,0]: HAMS outputs A, B, Fe at the global origin` (`+bem/+hams_mrel/run_at_draft.m:86`).
- **Interpolation** (`+bem/interpolate_at_draft.m`):
  - `draft_clamped = max(min(vertical_shift, max(drafts)), min(drafts))` (line 23).
  - Element-wise linear interpolation (lines 43-59).
  - Then `dz = z_cg_target - z_cg_hams` with `z_cg_hams = interp1(drafts, hydro_z_cg, draft_clamped)` (lines 80-83).
  - If `|dz| > 1e-4`, the matrices are transformed by `transform_to_cg(A_full, B_full, dz)` (lines 94-100).
- **The transform** (`+bem/transform_to_cg.m`): `T = [1 0 -Z_CG; 0 1 0; 0 0 1]` (lines 11-13) and `A_transformed = T' * A_symm * T` (line 41). As a result:
  - A11 and A33 are unchanged.
  - A15 becomes `A15 - dz*A11`.
  - A55 becomes `A55 - 2*dz*A15 + dz^2*A11`.
- **Consequence.** Within the node range [-1, 3.15], which equals `vertical_shift_bounds`, `dz = cg_body - centroid_z` does not depend on vs. So the CG transfer depends only on the mass distribution.
- **`trained` mode only.** `retransform_at_cg` overwrites the nearest node's matrices on a local config copy (`+optim/stage1_trained.m:77-82`; `+bem/retransform_at_cg.m:10-35`). It does not reach Stage 2, because `stage1_trained` returns no config.

### 2.4 Stage-3 evaluators (thin-shell and modular `evaluate_design_point`)

**Geometry grid.**
- `z = linspace(hull_z_min, hull_z_max, n_z)` (`+realise/+thin_shell/build_geometry_grid.m:49`).
- `A_outer` from `interp1(Aw_table_z, Aw_table, z)` (line 50).
- The inner contour is the outer contour offset inward by `offset_dist = t_steel / cos_alpha`, with `cos_alpha` floored at `1/max_slope_factor` (`+realise/+thin_shell/inner_properties_at_z.m:29-39`).
- Modular: the same, per strip, with the wall strip solid (`+realise/+modular_precast/build_geometry_grid.m:56-95`).

**Split at z_fill** (`+realise/+thin_shell/integrate_split.m:71-89`): below z_fill the full outer section is solid; above it there is a jacket annulus with air inside.

**Mass (thin shell)** (`+realise/+thin_shell/evaluate_design_point.m`):
- `M_steel = rho_shell*V_steel + (rho_fill - rho_shell)*V_fill`; `M_total = M_steel + M_air` (lines 47-49).
- `CG_z_body` (line 104).
- `Iyy_about_cg = max(0, Iyy_total_origin - M_total*CG_z_body^2)` (line 122).
- `CG_z_world = CG_z_body + vs` (line 147).

**Hydrostatics: the same tables as Stage 2, at `z_sub_top = min(-vs, hull_z_max)`** (lines 150-157). Then:
- `KM_world = CB_z_world + I_wp_yy/V_sub` (line 168)
- `GM = KM_world - CG_z_world` (line 175)
- `K33 = RHO_WATER*G*Aw` (line 178)
- `K55 = M_total*G*GM` if GM > 0 (lines 180-184)

**Added mass:** `[A11, A33, A55, ~, ~] = interpolate_at_draft(vs, cfg, out.CG_z_world)` (lines 189-191). **The full matrix is discarded, so A15 is not used.**

**Periods: uncoupled.**
- `T_heave = 2*pi*sqrt((M_total + A33)/K33_hydro)` (line 200).
- `T_pitch = 2*pi*sqrt((Iyy_about_cg + A55)/K55_hydro)` (line 205).

**Modular:** identical formulas with a single material density `rho_steel = rho_UHPC` (`+realise/+modular_precast/evaluate_design_point.m:9-127`; periods at lines 119 and 124).

## 3. Stage 1 `sweep`

### 1a. Tier 1 screen (no optimiser)

- **Draft grid.** `vs_grid = config.hydro_drafts(:)'` (`+optim/stage1_sweep.m:17`): the 11 BEM cache nodes. Each node is evaluated with `stage1_screen_draft` (line 42).
- **Inside `stage1_screen_draft`:**
  - `V_sub` comes from probing with `config.initial_densities` (`+optim/stage1_screen_draft.m:10-12`).
  - Default branch: `rho_unif = config.RHO_WATER*V_sub/V_strips` with `V_strips = sum(config.strip_V)` (lines 7 and 36), clamped to the density bounds, giving `rho_bal = repmat(rho_unif,1,N)` (lines 36-39).
  - Constructability branch: the wall is pinned and the rest is spread over the other strips with `rho_free = mass_need/max(V_free,eps)` (lines 23-35).
  - `fval = stage2_objective(x_bal, config)` (line 44), i.e. the 3-D Stage-2 objective.
  - `is_feasible = mass_err < 0.05 && props.GM_L >= config.gm_min` (line 47).
- **Initial densities matter only as a probe.** V_sub depends only on vs (2.1), so in sweep mode `initial_densities` affect only that probe and the no-grid fallback (`stage1_sweep.m:9-15`).

### 1b. Tier 2 refinement (fmincon)

- **Candidate selection.** Nodes are ranked by Tier-1 `|T_heave - T_heave_goal|` and the top `K_refine = config.n_sweep_refine` are kept (`stage1_sweep.m:58-62`). `n_sweep_refine = in.solver.n_sweep_refine = 4` (`WEC_User_Input.m:107`; `build_config.m:1022`).
- **Solver.** fmincon `'sqp'` with: Display none, MaxFunctionEvaluations 150, MaxIterations 25, ConstraintTolerance 1e-4, OptimalityTolerance 1e-4, StepTolerance 1e-6, ScaleProblem true (lines 76-84). One start per candidate; no multi-start.
- **Variables.** x = [vs, rho_1..rho_N].
  - Bounds come from `stage2_bounds` (line 24), but `lb_k(1) = ub_k(1) = vs_k` (lines 92-93), so vs is fixed.
  - Free variables: N (C1 default N = 5), or N-1 with a pinned wall.
  - `x0_k = max(lb_k, min(ub_k, sweep.x{k}))`, the Tier-1 point (line 94).
- **Objective and constraints.** Exactly the Stage-2 `obj_fun`/`con_fun` (lines 34-35); see section 6.
- **Candidate update.** A Tier-2 result replaces Tier 1 only if `fval_k < sweep.fval(k)` (lines 106-112). Its feasibility flag is `ef_k > 0 && mass_err < 0.01 && GM_L >= gm_min` (lines 101-103).
- **Selection.** The minimum `fval` among feasible candidates; if none is feasible, the overall minimum (lines 126-137).
- **Outputs.**
  - `x_opt = sweep.x{best}` (line 139).
  - `props_2d = properties_2d(x_opt, config)`: written to `stage1_2d.properties` but not used for any decision (line 140).
  - `conv_data = struct('final_objective', sweep.fval(best), 'final_exitflag', sweep.exitflag(best), 'sweep', sweep)` (lines 141-143).
  - Tier-2 iteration counts are not captured: only `[x_k, fval_k, ef_k]` are returned (line 97).
- **Physics:** properties_3d throughout (coupled periods; table hydrostatics; strip-table mass).
- **Counts:** same as Stage 2 (6 variables with vs fixed; 1 equality; 9 inequalities without a wall, 8 with one).

## 4. Stage 1 `skip`, `oneshot`, `trained` and the 2-D surrogate

### 1c. `skip`

- `x_opt_2d = x0_2d` (`+optim/run.m:115`).
- `props_2d = properties_3d(x0_2d, config)` (line 116). This is a 3-D property set stored in a field named `stage1_2d`.
- No optimisation.

### 1d. `oneshot`

- Creates a mass PID (`+optim/stage1_oneshot.m:10-12`). **It has no effect:** `solve_2d_surrogate(config, x0, ~)` ignores its third argument (`+optim/solve_2d_surrogate.m:1`).
- Makes one call: `[x_opt, props_2d, exitflag, conv_data] = solve_2d_surrogate(config, x0, mass_pid)` (`stage1_oneshot.m:14-15`).
- If `exitflag <= 0`, x_opt reverts to x0 (lines 17-22).
- `k_vol = k_vol_init` and `k_gm = k_gm_init` (`+optim/run.m:30-31`), both 1.0 (`WEC_User_Input.m:75-76`).

### 1e. `trained`: pre-calibration plus PID loop (not itself an optimiser)

**Pre-calibration** (`+optim/run.m:37-90`):
- `kvol_precal = vsub_3d_0/vsub_2d_0` at x0, clamped to `bounds_kvol` (lines 56-59).
- `kgm_precal = cgz_3d_0/cgz_2d_0`, clamped to `bounds_kgm` (lines 69-72). Only if both `|CG_z| > 0.01` (line 68).

**Loop** (`+optim/stage1_trained.m:41-245`), with `max_iters = config.max_outer_iterations = 50` (`WEC_User_Input.m:106`):
1. Run `solve_2d_surrogate(config, x0, mass_pid)` (line 47).
2. Optional HAMS enrichment, only when `run_HAMS_MREL` is true (`+optim/hams_enrichment_action.m:6-10`). It is false in `WEC_User_Input.m:116`.
3. Evaluate `p3 = properties_3d(x_opt, config)` (line 74). Retransform at the 3-D CG if `hams_dir` is non-empty (lines 77-82). `hams_dir` is always set when `in.bem.hams_dir` exists (`build_config.m:1200-1203`).
4. Volume update:
   - `vol_err_norm = (vsub_3d - vsub_2d)/vsub_3d` (line 122)
   - `k_vol = clamp(k_vol*(1 + damp_v*vol_pid.update(vol_err_norm,1)), bounds_kvol)` (lines 123-135)
5. GM/CG update:
   - `gm_err_norm = cg_z_3d/cg_z_2d - 1` (lines 140-141)
   - `k_gm = clamp(k_gm*(1 + damp_g*gm_pid.update(...)), bounds_kgm)` (lines 142-153)
6. PID details (`+optim/PID_Controller.m`): `P = Kp*e` (line 83); integral with clamp (86-89); `D = Kd*(e - e_prev)` (92-93); output saturation (99); anti-windup (103-105).

**Gains, limits and damping:**
- `vol_gains = [0.6, 0.01, 0.03]` and `gm_gains = [0.3, 0.02, 0.02]` (`WEC_User_Input.m:78-79`).
- Output limits: ±1 for both (lines 81-82).
- `bounds_kvol = [0.3, 3]` and `bounds_kgm = [0.5, 3]` (lines 83-84).
- Damping 0.5/0.7 for volume and 0.4/0.5 for GM, switching after iteration 2 (lines 85-89).

**Convergence:**
- `exitflag > 0 && ((|vol_err%| < 3 && |gm_2d - gm_3d| < 0.05) || (stable_count >= 2 && iteration > 3))` (`stage1_trained.m:159-168` and 200-202).
- `gm_2d = props_2d.GM`, which is k_gm-biased (line 86).
- Early exit after 3 consecutive iterations with both PIDs saturated (lines 186-198).
- Warm start: `x0 = x_opt` for the next iteration (lines 209-211).

**Hand-off:** only `x_opt` reaches Stage 2. The updated `k_vol`, `k_gm` and the retransformed hydro data live in the local config and are not returned (signature at line 1).

### 1f. The 2-D surrogate problem (`+optim/solve_2d_surrogate.m`)

**Solver.** fmincon `'sqp'`, Display off, MaxFunctionEvaluations 10000, ConstraintTolerance 1e-3, StepTolerance 1e-4, OptimalityTolerance 1e-5, MaxIterations 1000, OutputFcn (lines 54-62). Single start from x0 (line 66).

**Variables** x = [vs; rho_1..rho_N]:

| Variable | Units | Lower bound | Upper bound | Source |
|---|---|---|---|---|
| vs | m | `vertical_shift_bounds(1)` | `vertical_shift_bounds(2)` | lines 23-26; see section 6 for where these come from |
| rho_i | kg/m^3 | `ballast_density_bounds(1)`, raised to `per_strip_density_lb(i)` with constructability | `ballast_density_bounds(2)`; wall pinned to `rho_hull` | lines 23-26, 31-40 |

x0: `x0_2d = [config.initial_vertical_shift, config.initial_densities]` (`+optim/run.m:25`), or the previous iterate in `trained`.

**Objective** (lines 132-179):

```
r_gm    = (props.GM - gm_pref)/max(gm_half,1e-6)            % gm_pref = config.gm_target (0.5); props.GM is k_gm-biased
r_heave = (props.periods.heave - T_heave_goal)/max(heave_half,1e-6)   % coupled period
r_pitch = (props.periods.pitch - T_pitch_goal)/max(pitch_half,1e-6)   % coupled period
f = phi(r_gm) + phi(r_heave) + phi(r_pitch)
```

- Half-ranges are `0.5*(range(2) - range(1))` (lines 161-163).
- Guard: `f = 1e4` (a literal, not `config.penalty_guard`) if GM is NaN, `Aw < 1e-6`, a period is NaN or Inf, or an exception occurs (lines 137-147, 177).

**Inequalities c <= 0** (lines 183-248), in this order `[c_density_ratio; c_gm; c_monotonic; c_mass_min]`:
- `c_density_ratio = (max_density/min_density)/config.max_density_ratio - 1.0`. This is one global ratio over the platform strips. `min_density` is the minimum of platform densities above 50, or 100 if none exceed 50 (lines 206-213).
- `c_gm = 1.0 - props.GM/config.gm_min`, on the biased GM (line 222).
- `c_monotonic(j) = (densities(i+1) - densities(i))/rho_max` for adjacent pairs, skipping pairs that touch the wall (lines 226-238).
- `c_mass_min = 1.0 - props.mass_total/config.m_min_constructability`, only if `m_min_constructability > 0` (lines 241-246).

**Equality:** `ceq = props.mass_total/props.mass_buoyant_force - 1.0`, or 1.0 if buoyancy is at most 1e-6 (lines 253-257).

**Physics:** `properties_2d` (section 2.2).

**Counts (N = 5):** 6 variables; 1 equality; 6 inequalities without a wall (1 + 1 + 4), and 6 with a wall (1 + 1 + 3 + 1).

**Exception path:** returns x0 with `exitflag = -99` (lines 117-127).

**Dead code.** The OutputFcn `capture_convergence_data` is a non-nested local function that receives `convergence_data` by value (lines 51-52 and 279). Its updates never return. So:
- `best_feasible.props` stays `[]` and the "use best feasible" block (lines 103-112) never runs.
- `total_fmincon_iterations` is always 0 (line 94).

## 5. Stage 2 Phase A: density pre-conditioning (conditional)

- **Trigger:** `props_x0 = properties_3d(x0_3d)` and `if props_x0.GM_L < config.gm_min` (`+optim/run.m:182-184`).
- **Solver:** fmincon `'sqp'`, MaxFunctionEvaluations 500, MaxIterations 50, ConstraintTolerance 1e-4, OptimalityTolerance 1e-4, StepTolerance 1e-6, ScaleProblem true (lines 196-204). Single start.
- **Variables:** the Stage-2 bounds, but `lb_a(1) = ub_a(1) = x0_3d(1)` (lines 193-194), so vs is frozen. x0 = the Stage-1 point.
- **Objective and constraints:** the same as Stage 2 (`obj_fun_3d`, `con_fun_3d`, lines 175-176 and 206-207).
- **Acceptance:** the result is used only if `ef_a > 0 && props_a.GM_L >= config.gm_min`; then `x0_3d = x_a` (lines 211-219).
- **Counts:** the same as Stage 2.
- **C1:** not triggered in any of the three files. GM_L at the Stage-1 point is 0.2000064 / 0.2000064 / 0.2754 (all ≥ 0.2), taken from the properties_3d record `sweep.props{best}`. `stage2_3d.trajectory.x(:,1)` equals `stage1_2d.x_optimal` in all three.

## 6. Stage 2 main problem

**Solver.** fmincon with `config.stage2_algorithm = 'sqp'` (`WEC_User_Input.m:108`; `build_config.m:1015-1019`; `+optim/run.m:156-160`). Options: Display iter, MaxFunctionEvaluations 50000, MaxIterations 1000, ConstraintTolerance 1e-8, OptimalityTolerance 1e-8, StepTolerance 1e-10, FiniteDifferenceStepSize 1e-6, ScaleProblem true, OutputFcn (lines 163-173). Single start (lines 231-235).

**Variables** x = [vs, rho_1..rho_N], N = `in.geometry.num_ballast_sections = 5` (`WEC_User_Input.m:22`):

| Symbol | Code | Units | Lower bound | Upper bound | Source |
|---|---|---|---|---|---|
| vs | `x(1)` | m | `vertical_shift_bounds(1)` | `vertical_shift_bounds(2)` | `+optim/stage2_bounds.m:7-10`. `in.bounds.vertical_shift_bounds = []` (`WEC_User_Input.m:61`) selects auto bounds `[-hull_z_max+0.1, -hull_z_min-0.1]` (`build_config.m:934-939`). C1: [-1.0, 3.15]. |
| rho_i | `x(1+i)` | kg/m^3 | `ballast_density_bounds(1) = 20`; with constructability `max(20, per_strip_density_lb(i))` | 2500 | `WEC_User_Input.m:59`; `build_config.m:158`; `stage2_bounds.m:7-15` |
| wall rho | `x(1+w)` | kg/m^3 | `rho_hull` | `rho_hull` (pinned) | `stage2_bounds.m:16-20`; `rho_hull = in.materials.modular_precast.rho_hull = 2500` (`WEC_User_Input.m:48`) |

- **Constructability lower bounds.** `per_strip_density_lb` uses the offset-shell rule `rho_min_offset = (V_shell_at_tmin*rho_hull + V_int*rho_fill)/V_strip_i_est` (`build_config.m:525-550`). C1 values: [436.6, 227.4, 174.6, 369.4, 2500].
- **x0.** `x0_3d = x_opt_2d` from Stage 1 (`+optim/run.m:142`), or the Phase-A result.
- **Initial guess built at configuration time:**
  - `initial_vertical_shift = clamp(-hull_centroid(3))` (`build_config.m:1077-1079`).
  - `initial_densities = RHO_WATER*tanh-profile` with shape 3.0 (fallback 1.5), clamped to the bounds and per-strip lower bounds, wall pinned (`build_config.m:1081-1126`).

**Objective** (`+optim/stage2_objective.m`):

```
props = properties_3d(x, config)                                   (l.7)
if isnan(GM_L) or any period NaN/Inf: f = config.penalty_guard (1e4)   (l.10-18)
gm_half = 0.5*(gm_range(2)-gm_range(1)); heave_half=...; pitch_half=...   (l.21-23)
r_gm    = (props.GM_L - config.gm_target)/gm_half                  (l.26)
r_heave = (props.periods.heave - config.T_heave_goal)/heave_half   (l.27)  % coupled
r_pitch = (props.periods.pitch - config.T_pitch_goal)/pitch_half   (l.28)  % coupled
f = range_penalty(r_gm,k_amp)+range_penalty(r_heave,k_amp)+range_penalty(r_pitch,k_amp)   (l.31-34)
```

- **Target values for C1:**
  - `gm_target = 0.5` and `gm_range = [0.2, 0.7]`, so `gm_half = 0.25`.
  - `T_heave_goal = 7.77` and `T_heave_range = [7, 10]`, so `heave_half = 1.5`.
  - `T_pitch_goal = 3.89` and `T_pitch_range = [3, 5]`, so `pitch_half = 1.0`.
  - Sources: `WEC_User_Input.m:64-70`; `build_config.m:876-883`.
- **Flat region of the objective.** When `GM_L <= 0`, K55 is 0 (`properties_3d.m:189-193`). That makes pitch Inf, so `f = 1e4`. This is the source of the flat 1e4 values at Tier-1 nodes with negative GM.

**Equality:** `ceq = props.mass_total / props.mass_buoyant_force - 1.0` (`+optim/stage2_constraints.m:52`).

**Inequalities** (`stage2_constraints.m`), assembled as `c = [c_gm; c_ratio; c_mono; c_mass_min]` (line 51):
- `c_gm = 1.0 - props.GM_L / config.gm_min`, with `gm_min = 0.2` (line 20).
- For each adjacent pair i, excluding pairs where `i == w_idx` or `i+1 == w_idx` (lines 26-31):
  - `c_ratio(j) = densities(i) / (densities(i+1) + 1) - config.max_density_ratio`, with `max_density_ratio = 100` (`WEC_User_Input.m:60`) (line 38). This can bind only when `rho_i/(rho_{i+1}+1) > 100`, e.g. `rho_{i+1} < 24` for `rho_i = 2500`.
  - `c_mono(j) = (densities(i+1) - densities(i)) / rho_max` (line 39): density must not increase upward.
- `c_mass_min = 1.0 - props.mass_total / config.m_min_constructability`, only if `m_min_constructability > 0` (lines 45-49). That value is set only with constructability, from `sum V_strip_i*per_strip_density_lb(i)` (`build_config.m:579-609`); otherwise it is 0 (line 647). C1 modular: 7335.49 kg.
- **On exception:** `c = ones(...)` and `ceq = 1` (lines 53-66).

**Physics:** properties_3d (section 2.1).

**Counts (N = 5):** 6 variables (5 free with a wall); 1 equality; 9 inequalities without a wall (1 + 4 + 4), and 8 with a wall (1 + 3 + 3 + 1).

**Post-solve (diagnostic only).** `check_3d_convergence` (`+optim/check_3d_convergence.m:7-48`) sets `converged` from:
- monotonicity (exact sign)
- `|M - Mb| < 10 kg`
- `GM_L > gm_min` (strict)
- exitflag 1, or exitflag 0/2 with `constrviolation < 1e-6` and `firstorderopt < 1e-2`

Periods (10% of goal) are stored as `periods_acceptable` but are not part of `converged`.

## 7. Stage 3

### 3p. `preliminary`

`results = opt_results; final_props = opt_results.Final3D` (`+realise/+preliminary/run.m:6-7`). There is no optimisation.

### 3t. `thin_shell` (`+realise/+thin_shell/solve.m`; called from `+realise/+thin_shell/run.m:33/36`)

**Solver.** fmincon `'sqp'`, Display iter, StepTolerance 1e-8, OptimalityTolerance 1e-6, ConstraintTolerance 1e-6, MaxIterations 200, MaxFunctionEvaluations 1000 (lines 243-250). Single start (lines 253-254).

**Variables** x = [vs; t_steel; z_fill]:

| Variable | Units | Lower bound | Upper bound | x0 |
|---|---|---|---|---|
| vs | m | `vertical_shift_bounds(1)` = -1.0 (lines 117, 236) | `vertical_shift_bounds(2)` = 3.15 (lines 118, 237) | `vs_opt = x_opt_3d(1)` from Stage 2 (lines 119, 238) |
| t_steel | m | `t_min = config.steel_t_min = in.materials.thin_shell.t_min = 0.025` (`WEC_User_Input.m:44`; `build_config.m:280`; line 49) | `0.95*t_max`, with `t_max = 0.5*min r` over 9 probes in the middle 80% of hull height (lines 76-102, 237). C1: t_max = 0.50249, ub = 0.47737. | `t_init = max(min(0.02, 0.9*t_max), max(t_min, 1e-4))` = 0.025 for C1 (`WEC_User_Input.m:43`; line 110) |
| z_fill | m | `hull_z_min + 1e-3` | `hull_z_max - 1e-3` | analytic seed solving `M(z_fill) = RHO_WATER*V_sub(-vs_opt)` at t_init (lines 147-214) |

x0 is clamped into the bounds (line 241).

**Densities.**
- `rho_shell = config.rho_shell` (line 45) and `rho_fill = config.rho_fill` (line 46). Both are 7500 in C1 (`WEC_User_Input.m:29-31`).
- `rho_air = 1.2` (line 47; `WEC_User_Input.m:30`).

**Objective** (lines 444-463):

```
r_h = (r_.T_heave - T_heave_goal)/max(heave_half,1e-6);  r_p = (r_.T_pitch - T_pitch_goal)/max(pitch_half,1e-6)
f = range_penalty(r_h,k_amp) + range_penalty(r_p,k_amp)        % UNCOUPLED periods; no GM term
```

Guard: `f = 1e4` if more than 5% of grid sections are degenerate or a period is non-finite (lines 449-455).

**Equality:** `ceq = r_.M_total / r_.mass_buoyant_force - 1.0`, where `mass_buoyant_force = RHO_WATER*V_sub(-vs)` (line 479; `evaluate_design_point.m:163`).

**Inequality:** `c = 1.0 - r_.GM / ctx.config.gm_min` (line 481). Guard: `c = 1`, `ceq = 1` (lines 471-477, 483).

**Physics:** section 2.4. Hydrostatics are the same tables as Stage 2; mass, CG and Iyy come from the geometry grid with `n_z_grid = 300` (`WEC_User_Input.m:46`); periods are uncoupled.

**Counts:** 3 variables; 1 equality; 1 inequality.

**Acceptance and outputs.**
- `t_min_active = ((t_star - t_min)/max(t_max - t_min, eps) < 0.01)` (line 260). This is a 1%-of-range proximity flag, not a bound-activity test.
- The final evaluation is packaged into `steel_data` (lines 263-363).
- `steel_data.feasible = realised.feasible = isfinite(GM) && isfinite(T_heave) && isfinite(T_pitch)` (line 361; `evaluate_design_point.m:210`).
- `final_props` is rebuilt if that flag is true (`thin_shell/run.m:45-47`). The rebuild (`+realise/build_realised_properties.m:7-156`):
  - re-interpolates added mass at the realised CG (lines 94-96)
  - recomputes **coupled** periods (lines 145-152)
  - takes GM, mass, CG and Iyy from `steel_data` (lines 46-57)

### 3m. `modular_precast` (`+realise/+modular_precast/solve.m`; called from `solve_and_extract.m:116-117`)

**Pre-check (hard error, not an optimisation).** `final_props.mass_total` from Stage 2 must lie within `[0.95*M_min_uhpc, 1.05*M_max_uhpc]` (`solve_and_extract.m:88-106`).
- `M_min_uhpc` uses the wall solid and all other strips at t_min.
- `M_max_uhpc` uses fully solid UHPC.

**Solver.** fmincon `'sqp'`, Display iter, StepTolerance 1e-8, OptimalityTolerance 1e-6, ConstraintTolerance 1e-6, MaxIterations 300, MaxFunctionEvaluations 2000, OutputFcn logging only (lines 197-205). Single start (lines 208-209).

**Variables** x = [vs; z_fill; t_i for i in nw_idx], with `nw_idx = setdiff(1:N_strips, wall_strip_idx)` (line 24). The wall is the last strip unless the deck name contains `_180` (`solve_and_extract.m:58-62`; `build_config.m:48-52`).

| Variable | Units | Lower bound | Upper bound | x0 |
|---|---|---|---|---|
| vs | m | -1.0 (lines 95, 172) | 3.15 (lines 96, 173) | `vs_opt = x_opt_3d(1)` (lines 70, 174) |
| z_fill | m | `hull_z_min + 1e-3` | `hull_z_max - 1e-3` | analytic seed (lines 107-146) |
| t_i, non-wall strip i | m | `t_min = constructability_t_min = 0.0762` (`WEC_User_Input.m:50`; `build_config.m:325`; `solve_and_extract.m:20`) | `0.95*t_max` (lines 76-88, 173) | `t_init = uhpc_t_init = 0.02*7500/2500 = 0.06` (`WEC_User_Input.m:53-54`; `solve_and_extract.m:24-29`), clamped to `max(min(t_init, 0.9*t_max), max(t_min, 1e-4))` = 0.0762 (line 93), i.e. starts on the lower bound |

**Densities:** `rho_uhpc = opts.rho_steel = rho_UHPC = 2500` (`solve_and_extract.m:18`; `solve.m:56`) and `rho_air = 1.2` (`WEC_User_Input.m:49`).

**Objective:** identical form to thin shell, phi(r_h) + phi(r_p) with **uncoupled** periods and no GM term (lines 376-397). Guard 1e4 (lines 382-389).

**Constraints:**
- `ceq = r_.M_total/r_.mass_buoyant_force - 1.0` (line 414).
- `c = 1.0 - r_.GM/ctx.config.gm_min` (line 416).
- Guards: `c = 1`, `ceq = 1` (lines 405-418).

**Physics:** section 2.4, modular variant. The wall strip is solid; the region below z_fill is solid through `integrate_split` (`modular_precast/evaluate_design_point.m:9-12`).

**Counts:** 2 + (N-1) variables (6 for C1); 1 equality; 1 inequality.

**Post-solve.**
- Strips wholly below `z_fill*` are promoted to solid (`solve.m:222-229`).
- `extract_strip_geometry` inherits `cstr = solve_data` (grep: `extract_strip_geometry.m:8`) and adds per-strip diagnostics (lines 324-385).
- `final_props` is rebuilt with `build_realised_properties` (`modular_precast/run.m:26-27`; `modular_precast/build_realised_properties.m:13`) when `feasible` is true; same finiteness-only flag at `evaluate_design_point.m:129`.

## 8. Dependency maps (as implemented)

**Stage 2 / Phase A / sweep Tier 2 (properties_3d):**

| Quantity | vs | rho_i |
|---|---|---|
| M | no: `strip_V` fixed (l.121-131) | yes |
| CG_z body | no | yes (l.142) |
| CG_z world | yes: `+vs` (l.146) | yes |
| Iyy about CG | no (l.151) | yes |
| Aw, V_sub, I_wp_yy, CB_body | yes: tables at `min(-vs, hull_z_max)` (l.41-45) | no |
| CB world, KM | yes (l.46, 106) | no |
| GM_L | yes, via tables only; the vs translation cancels (l.106, 146, 180) | yes, via CG |
| A11, A33 | yes: interpolation at vs (`interpolate_at_draft.m:43-59`) | no (invariant under T) |
| A55, A15 | yes (interpolation) | yes, via `dz = cg_body - centroid` (`interpolate_at_draft.m:80-100`) |
| K33 | yes (Aw) | no |
| K55 | yes | yes (`M*g*GM`) |
| T_heave | yes | yes (M); equals `2*pi*sqrt((M+A33)/(rho*g*Aw))` in the C1 data |
| T_pitch, coupled | yes | yes; equals `2*pi*sqrt((Iyy+A55-A15^2/(M+A11))/K55)`, verified on the C1 data |

**Stage 1 2-D surrogate (properties_2d):**
- **vs** changes:
  - sub_area, CB, Aw, I_wp, BM and KM (l.18, 28-30, 38-78)
  - V_sub as a step function (l.131-133, 175-178)
  - CG world (l.170)
  - GM via KM, and via `k_gm*CG_world` when k_gm ≠ 1 (l.204-205)
  - A* (l.302-304)
- **vs does not change** M or Iyy (strip body positions are fixed: l.119, 137).
- **rho_i** changes M, CG, Iyy, GM, K55, A55 and A15. Densities are linearly interpolated between nodes (l.149-150).

**Stage 3 thin shell:**
- **vs** changes Aw, V_sub, I_wp_yy, CB, KM (tables, `evaluate_design_point.m:150-157`), CG world (l.147), A11/A33/A55 (l.189-191), GM, T_h and T_p. It does not change M, body CG or Iyy.
- **t_steel and z_fill** change M, CG (body and world), Iyy, GM, A55 (via CG; A11/A33 invariant), K55 and both periods. They do not change Aw, V_sub, CB or KM.
- **A15 is computed but discarded** (`~` at l.189), so it does not enter the Stage-3 objective.

**Stage 3 modular:** as thin shell. t_i acts only on strip i above z_fill (`modular_precast/build_geometry_grid.m:80-94`; `integrate_split.m:71-89`). **A t_i for a strip wholly below z_fill has no effect.** In C1, strip 1 (top -2.6125 < z_fill* -2.1846) is such a strip at the solution.

## 9. Hand-offs

**Stage 1 → Stage 2.** Only `x_opt_2d`, used as `x0_3d` (`+optim/run.m:142`). `props_2d`, `conv_data_2d` and `hist` are only stored in `results.stage1_2d` (`report_assemble_results.m:13-45`). `stage1_2d.iterations = length(hist.mass_errors_history)` (l.16) is 0 in sweep mode.

**Stage 2 → Stage 3.**
- What is passed: `x_opt = x_opt_3d` (`run.m:329`) and `opt_results`. `opt_results.Final3D = final_props_optimiser = properties_3d(x_opt_3d)` plus `cross_section` (`run.m:259-270`; `report_assemble_results.m:66`).
- **What Stage 3 actually reads from Stage 2:**
  - `x_opt(1)` only, as the vs warm start and for the z_fill seed (`thin_shell/solve.m:119, 147-150, 238`; `modular_precast/solve.m:70, 120-123, 174`).
  - Modular only: `final_props.mass_total` for the pre-check (`solve_and_extract.m:95-106`).
  - Reporting only: `final_props.GM_L`, `.periods.heave/.pitch` and `.mass_total`, used as "targets" and "residuals" (`thin_shell/solve.m:124-126, 348-358`; `modular_precast/solve.m:100-102, 250-254, 305-308`).
  - Carried over only: `cross_section`, `densities_at_nodes`, `components`, `A_sub` (`build_realised_properties.m:39-41, 236-245, 306-307`).
- **What Stage 3 does not read: the Stage-2 densities `x_opt(2:end)`.** Thin shell checks only `length(x_opt_3d) >= 2` (`thin_shell/solve.m:37`).

**Stage 3 output.**
- `results.stage2_3d.properties` is overwritten with the realised properties (`thin_shell/run.m:60`; `modular_precast/run.m:37`).
- `results.Final3D` keeps the Stage-2 properties.
- `results.steel_data` holds the thin-shell output (`thin_shell/run.m:59`); `results.constructability` holds the modular output (`modular_precast/run.m:36`).

## 10. C1 stored results (Output/C1_*_results.mat, loaded with scipy)

- Constraint values were recomputed from the stored properties with the code formulas. The recomputed objectives match the stored `fval` to at least 10 significant digits.
- "Active" means within 1e-6 relative of the bound, or |c| < 1e-6.
- `preliminary` and `thin_shell` have identical Stage 1 and Stage 2 results.

| Stage / file | Design vector | Objective | Constraints at solution | exitflag / iterations |
|---|---|---|---|---|
| S1 sweep, prelim and thin | vs = 0.7786 (node 5); rho = [1454.48, 1433.39, 387.98, 387.98, 387.98] | `final_objective` 15.0682 (3-D) | ceq -3.2e-14. c_gm -3.21e-5: within the Tier-2 tolerance, not active at 1e-6. c_mono(3,4) = c_mono(4,5) = 0: active. Ratios ≈ -96 to -99. No bounds active. | `final_exitflag` 1; Tier-2 iterations not stored; `stage1_2d.iterations` 0 |
| S1 sweep, modular | vs = 1.075 (node 6); rho = [1312.405, 1312.405, 703.86, 369.369, 2500] | 28.3742 | ceq 2.9e-15. c_gm -0.377. c_mono(1,2) -2.2e-8: active. rho_4 at its lower bound 369.369: active. rho_5 pinned. c_mass_min -1.659. | 1; not stored |
| S2, prelim and thin | vs = 0.893085; rho = [2500, 1297.195, 310.481, 310.481, 310.481] | 0.792349 (r_gm -0.785, r_h 0.0047, r_p -0.419) | ceq 7.8e-14. c_gm -0.518. rho_1 at upper bound: active. c_mono(3,4) = c_mono(4,5) = 0: active. Ratios -95.8 to -99.0. GM 0.3036; T_h 7.777 (coupled = uncoupled); T_p 3.471 coupled (uncoupled 3.548). | 2; 57 iterations, 422 function evaluations; `converged` 1 |
| S2, modular | vs = 0.982815; rho = [2500, 1430.025, 369.369, 369.369, 2500] | 9.15330 (r_gm -1.2, r_h -1.838, r_p 1.137) | ceq 0. c_gm -2.2e-16: active (GM = 0.2). rho_1 at upper bound and rho_4 at lower bound: active. Wall pinned. c_mono(3,4) = 0: active. c_mass_min -1.779. T_h 5.013; T_p 5.027 coupled. | 1; 32 iterations, 238 function evaluations; `converged` 1 |
| S3 thin shell | vs = 0.893314; t = 0.0272603; z_fill = -2.730626 | `phi_star` 3.38e-12 (uncoupled T_h 7.77000, T_p 3.89000) | ceq ≈ -8.8e-16. c_gm -1.198 (GM 0.4395). No bounds active: t is 9% above t_min = 0.025; the stored `t_min_active = 1` is the 1%-of-range flag. M = 20660.41 kg. Reported (coupled) T_p = 3.8133. | 2; 8 iterations (from `Output/thin_shell/steel_solve.log:17`; not in the .mat) |
| S3 modular | vs = 0.893224; z_fill = -2.184610; t_1..t_4 = 0.0762 | `phi_star` 1.16438 (r_h 0.0019, r_p 1.070; uncoupled T_p 4.960) | ceq -3.3e-16. c_gm -0.154 (GM 0.2308). All four t_i at lower bound: active; t_1 has no effect, strip 1 is below z_fill. M = 20660.68 kg (+1.34% vs Stage 2). Reported coupled T_p 4.9065. | 2; 7 iterations (`iter_history.iter` ends at 7) |

## 11. Discrepancies

### Physics and formulation mismatches between stages (code facts)

1. **Stage 2 and Stage 3 optimise different period definitions.**
   - Stage 2 (and sweep Tier 1/2) uses coupled periods (`properties_3d.m:253-257`; `stage2_objective.m:27-28`).
   - The Stage-3 objectives use uncoupled periods without A15 (`thin_shell/evaluate_design_point.m:189-205`; `modular_precast/evaluate_design_point.m:109-124`).
   - Final reported periods are coupled again (`build_realised_properties.m:145-152`).
   - C1 thin shell: Stage 3 hits pitch 3.8900 exactly (uncoupled), but the reported `final_props` pitch is 3.8133.
2. **Stage 2 has a GM objective term; Stage 3 does not.**
   - Stage 2 includes phi(r_gm) toward `gm_target = 0.5` (`stage2_objective.m:26`).
   - Stage 3 has only the GM floor (`thin_shell/solve.m:481`; `modular_precast/solve.m:416`).
   - Realised GM: 0.4395 (thin shell) and 0.2308 (modular).
3. **Stage-3 "targets/residuals" compare different quantities.** They compare Stage-3 uncoupled periods (`T_heave_realised` / `T_pitch_realised` = uncoupled) with Stage-2 coupled periods (`thin_shell/solve.m:348-356`; `modular_precast/solve.m:250-308`). Example: thin shell `dT_pitch_pct = +12.07%`, which is 3.8900 against 3.4711.
4. **Different mass models.**
   - Stage 2: uniform density per horizontal strip on precomputed strip integrals (`properties_3d.m:120-136`).
   - Stage 3: shell, fill and air on a 300-point grid (`build_geometry_grid.m:49-78`).
   - Stage 3 does not use the Stage-2 densities at all (section 9). Its realised per-strip equivalent densities violate Stage-2 rules:
     - thin shell `strip_rho_eff = [6977, 285.5, 399.6, 2152.8, 2275.4]`: above the 2500 bound and non-monotone
     - modular `[2500, 1638.7, 211.0, 481.2, 2500]`: non-monotone
   - Stage 3 has no monotonicity or density constraint.
   - Hydrostatics are consistent between the stages: the same tables (`properties_3d.m:41-45` against `evaluate_design_point.m:153-156`).
5. **Stage-2 modular feasible set differs from Stage 3's.**
   - Stage 2 stopped at vs = 0.983 with GM at the floor and T_h = 5.01 s (wall pinned, rho_4 at its per-strip lower bound, c_mono(3,4) active).
   - Stage 3 then re-chose vs = 0.893, reaching T_h = 7.77 with GM 0.231.
   - Observation from the data; I did not trace which Stage-2 constraint causes this.
6. **The Stage-1 2-D surrogate (oneshot/trained) differs from the 3-D model.** It has:
   - a 5-strip, step-function V_sub (`properties_2d.m:84, 131-133`)
   - linearly interpolated densities rather than per-strip constants (l.149-150)
   - a k_gm-biased GM in the constraint, the objective and K55 (l.204-205, 283-284), while added mass uses the raw CG (l.302-304)
   - one global max/min density-ratio constraint instead of adjacent pairs (`solve_2d_surrogate.m:206-213` against `stage2_constraints.m:38`)
7. **In sweep mode `stage1_2d.properties` is `properties_2d`, recorded but never used** (`stage1_sweep.m:140`). It contradicts the 3-D numbers at the same x: C1 2-D GM -0.744, mass 16242 against buoyancy 8079; 3-D GM 0.200, mass = buoyancy = 20827. In `skip` mode the same field holds `properties_3d` (`run.m:116`).
8. **Stage-3 acceptance is finiteness only.** `feasible = isfinite(GM) && isfinite(T_heave) && isfinite(T_pitch)` (`thin_shell/evaluate_design_point.m:210`; `modular_precast/evaluate_design_point.m:129`). This gates the property swap (`thin_shell/run.m:45`; `build_realised_properties.m:7`) without checking exitflag, ceq or c.
9. **Guard values differ.** Stage 2 uses `config.penalty_guard` (`stage2_objective.m:11, 16`). The 2-D path and Stage 3 hard-code 1e4 (`solve_2d_surrogate.m:138, 146, 177`; `thin_shell/solve.m:450, 454`; `modular_precast/solve.m:383, 388`). They are equal for C1.

### Defined but no effect (code facts)

- **Mass PID:** `in.pid.mass_gains` and `mass_limits` (`WEC_User_Input.m:77, 80`) build a PID that `solve_2d_surrogate` ignores (`stage1_oneshot.m:10-12`; `stage1_trained.m:23-24`; `solve_2d_surrogate.m:1`).
- **2-D OutputFcn / best-feasible logic** never takes effect (`solve_2d_surrogate.m:51-52, 103-112, 279`).
- **`config.shell` is always `[]`** (`build_config.m:283`), so the shell branches in `properties_2d.m:93-96, 157-167, 259-267` are dead.
- **Unread config fields:** `config.mass_correction_factor`, `gm_correction_factor` and `z_cg_target` (`build_config.m:953-967`) are never read elsewhere in src (grep).
- **Trained-mode corrections never reach Stage 2:** `k_vol`/`k_gm` and the retransformed hydro data are local (section 4).
- **Stage-2 densities have no effect on Stage 3** (section 9).
- **Modular t_i for strips below z_fill have no effect** (section 8).
- **Two-density thin-shell split has no mass effect for C1**, because `rho_fill = rho_shell = 7500` (`WEC_User_Input.m:31`).
- **Unused inputs:** `modular_precast/evaluate_design_point.m` takes `t_offset_strip` and `is_solid_strip` and does not use them (l.139). `B_full` is never used in any period.

### Code comments that disagree with code

- `stage2_objective.m:3` says the residual divides by `max(half_range, eps)`. Code fact: no max (l.26-28).
- `build_realised_properties.m:133-135` says "the periods the objective reads come from the coupled eigenproblem below". Code fact: the Stage-3 objectives read uncoupled periods (`thin_shell/solve.m:452-462`).
- `WEC_User_Input.m:69-70` calls `gm_range` a "reporting/consistency band" and `gm_target` "documentation/objective reference". Code fact: both define the GM objective term in Stages 1 and 2 (`stage2_objective.m:21, 26`; `solve_2d_surrogate.m:161, 166-167`).
- `modular_precast/solve.m:46` cites `solve_and_extract.m:38` for `opts.rho_steel`. Code fact: line 18.

### docs/ against code

- **Cache reference point.** Docs claim: "Cached coefficients are referenced to the BEM/HAMS center of gravity" (`docs/METHODS_ENGINE.md:528`). Code fact: the cache matrices are treated as "6×6 at origin" and transformed to `z_cg` by the code (`build_config.m:835-845`; `run_at_draft.m:86`).
- **`stage1_2d.properties`.** Docs claim: "2-D surrogate properties at `x_optimal`" (`docs/RESULT_SCHEMA.md:69`). Code fact: `skip` stores `properties_3d` (`run.m:116`).
- **`stage1_2d.iterations`.** Docs claim: it is the "Stage-1 execution status" count (`docs/RESULT_SCHEMA.md:70`). Code fact: it is `length(hist.mass_errors_history)` (`report_assemble_results.m:16`), which is 0 for sweep even though Tier-2 fmincon runs.
