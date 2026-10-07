# Result MAT-file schema

Each run writes `Output/<hull>_<realisation_type>_results.mat`. Each MAT file contains exactly two top-level variables:

| Variable | Purpose |
| --- | --- |
| `results` | Complete record of configuration, optimisation, and realisation data. |
| `final_props` | Final design properties for downstream analysis and reporting. |

```matlab
S = load('Output/<hull>_<realisation_type>_results.mat');
results = S.results;
final_props = S.final_props;
```

## Conventions

- SI units apply unless stated otherwise.
- `vertical_shift` is signed [m]. Positive values move the hull upward. `draft = abs(hull_z_min + vertical_shift)`.
- World-frame `z = 0` is the still-water surface, with positive `z` upward.
- Reduced dynamic matrices use `[surge, heave, pitch]` order. `MassMatrix_CG` and `MassMatrix_Origin` use `[surge, sway, heave, roll, pitch, yaw]` order.
- `KM`, `CB(3)`, and `CG_total(3)` are world-frame elevations. `GM_L = KM - CG_total(3)`.
- `Inf` in a period field denotes no restoring stiffness. Empty values denote fields that do not apply to the selected realisation.

## `results`

| Field | Type | Contents | Used for |
| --- | --- | --- | --- |
| `config` | struct | Runtime configuration and precomputed geometry and hydrodynamic data. | Reproducing and interpreting the run. |
| `context` | struct | Site and target metadata. | Reporting metadata. |
| `stage1_2d` | struct | Stage-1 surrogate result and history. | Draft screening and Stage-1 review. |
| `stage2_3d` | struct | Stage-2 solution, convergence data, and final properties. | Optimisation review. |
| `Final3D` | struct | Theoretical Stage-2 design before realisation. | Comparing theoretical and realised designs. |
| `stage3` | struct or `[]` | Stage-3 realised design (thin shell and modular precast); `[]` for `preliminary`. | Realised-design analysis. |
| `optimization_time` | double, s | Stage-1 and Stage-2 elapsed time. | Runtime reporting. |
| `schema_version` | char | Result-schema version. | Schema identification. |
| `realisation_type` | char | `'preliminary'`, `'thin_shell'`, or `'modular_precast'`. | Selecting mode-specific fields. |

### `results.config`

| Fields | Convention | Used for |
| --- | --- | --- |
| `RHO_WATER`, `G`, `water_depth` | Water density [kg/m³], gravity [m/s²], and water depth [m]. `water_depth = -1` denotes deep water. | Hydrostatic and hydrodynamic interpretation. |
| `ms2_file`, `ms2_deck_sha256` | Geometry-deck path and SHA-256 digest. | Geometry provenance. |
| `hull_z_min`, `hull_z_max`, `density_nodes_z`, `strip_edges` | Body-frame geometry elevations [m]. | Draft and strip interpretation. |
| `ballast_density_bounds`, `gm_min`, `gm_range`, `gm_target`, `T_heave_goal`, `T_pitch_goal`, `T_heave_range`, `T_pitch_range` | Design constraints and targets. | Comparing outputs with design settings. |
| `Aw_table_z`, `Aw_table`, `I_wp_xx_table`, `I_wp_yy_table`, `V_sub_table`, `CB_z_table`, `S_wet_table`, `P_table` | Precomputed geometry tables indexed by body-frame elevation. | Hydrostatic interpolation. |
| `strip_V`, `strip_CB_z`, `strip_Ixx`, `strip_Iyy`, `strip_Izz` | Per-strip volume, centroid, and inertia data. | Mass-property reconstruction. |
| `hydro_table`, `hydro_cache` | BEM coefficient cache and metadata. | Added-mass and radiation-damping lookup. |
| `added_mass_diagonal`, `radiation_damping_diagonal`, `added_mass_full`, `radiation_damping_full` | Hydrodynamic coefficients at cached shifts. Reduced DOF order is `[surge, heave, pitch]`. | Dynamic-response calculations. |
| `output` | Output-option snapshot. | Identifying enabled exports. |

### `results.context`

| Fields | Units | Used for |
| --- | --- | --- |
| `T_heave_goal`, `T_pitch_goal`, `T_surge_goal` | s | Recorded period targets. |
| `T_heave_range`, `T_pitch_range`, `T_surge_range` | s, `[lo, hi]` | Recorded target bands. |
| `gm_target` | m | Recorded longitudinal metacentric-height target. |
| `WIS_station`, `data_year` | char | Optional site metadata. |
| `water_depth` | m, positive downward | Recorded BEM water depth. |

### `results.stage1_2d`

| Fields | Type / units | Used for |
| --- | --- | --- |
| `x_optimal` | `[vertical_shift, rho_1, ..., rho_N]` | Stage-1 design vector. |
| `properties` | struct | 2-D surrogate properties at `x_optimal`. |
| `iterations`, `converged`, `pid_saturated` | count and logical values | Stage-1 execution status. |
| `mass_errors`, `gm_errors`, `mass_corrections`, `gm_corrections` | vectors | Stage-1 history. |
| `mass_2d_history`, `mass_3d_history`, `mass_2d_corrected_history`, `gm_2d_history`, `gm_3d_history`, `cg_z_2d_history`, `cg_z_3d_history` | vectors | Surrogate comparison. |
| `R2_mass`, `R2_GM`, `R2_cg`, `MAPE_mass`, `MAPE_GM`, `MAPE_cg`, `ME_mass`, `ME_GM` | double | Stage-1 comparison statistics. |
| `convergence_data`, `convergence_metrics` | struct | Objective, exit, sweep, and convergence records. |

`results.stage1_2d.properties` uses the same principal property names as `results.Final3D`, with `GM` and `GM_uncorrected` for the 2-D surrogate.

### `results.stage2_3d`

| Field | Type / units | Meaning | Used for |
| --- | --- | --- | --- |
| `x_optimal` | `[vertical_shift, rho_1, ..., rho_N]` | Stage-2 design vector. | Optimisation record. |
| `properties` | struct | Final properties for the selected realisation. | Run review. |
| `exitflag`, `fval`, `output` | double, double, struct | SQP exit status, objective value, and solver output. | Optimisation review. |
| `converged`, `quality_metrics` | logical, struct | Mass, GM, period, monotonicity, and solver status. | Acceptance checks. |
| `iteration_errors` | struct | Mass, GM, heave, and pitch errors by iterate. | Convergence plots. |
| `trajectory` | struct | Stage-2 design and property history. | Post-processing. |

`quality_metrics` contains `monotonic`, `mass_balance`, `mass_balance_error_kg`, `GM_satisfied`, `GM_margin`, `periods_acceptable`, `heave_error_pct`, `pitch_error_pct`, `fmincon_optimal`, `fmincon_acceptable`, `exitflag`, `firstorderopt`, `constrviolation`, and `iterations`.

`trajectory` contains `x`, `n`, `vs`, `draft`, `mass_total`, `mass_buoy`, `GM`, `CG_z`, `T_heave`, `T_pitch`, and `fval`, with one value or column per accepted SQP iterate.

### `results.Final3D`

`results.Final3D` stores the theoretical Stage-2 design before realisation. It is the base schema for `final_props`.

| Fields | Units / convention | Used for |
| --- | --- | --- |
| `vertical_shift`, `draft` | m | Hull position and physical draft. |
| `Aw`, `I_wp_xx`, `I_wp_yy`, `V_sub`, `A_sub` | m², m⁴, m³, m² | Hydrostatic geometry. |
| `CB`, `KM`, `GM_L`, `CG_total` | m, world frame | Stability and mass-property reporting. |
| `mass_buoyant_force`, `mass_total`, `mass_discrepancy` | kg | Buoyancy and mass balance. |
| `Ixx`, `Iyy`, `Izz`, `Inertia_Tensor` | kg·m² | Rigid-body inertia. |
| `A11`, `A33`, `A55`, `A_full` | kg, kg, kg·m² | Added mass. |
| `B_full` | N·s/m and N·m·s/rad | Radiation damping. |
| `K_hydro`, `K_pto`, `K_total` | Dynamic stiffness matrices. `K_pto` is zero. | Natural-period calculation. |
| `periods`, `coupled_periods`, `coupled_modes`, `participation_factors`, `surge_per_pitch` | s, mode-shape, and percent data | Coupled-mode analysis. |
| `components`, `densities_at_nodes`, `cross_section` | Density and geometry data | Mass-distribution plots and downstream geometry use. |
| `MassMatrix_CG`, `MassMatrix_Origin` | 6x6 rigid-body matrices | Coupled rigid-body analysis. |

### `results.stage3`

`results.stage3` stores the realised Stage-3 design of both realising types. It always describes
the realised design, accepted or failed; Stage 3 never falls back to the Stage-2 design.

| Field | Type / units | Meaning |
| --- | --- | --- |
| `mode`, `hull_name` | char | Realisation type and deck stem. |
| `status`, `reason` | char | `'accepted'` or `'failed'`; `reason` names the failed checks (empty when accepted). |
| `escalation` | char | Last Stage-3 step run: `'split'`, `'fixed_draft'`, `'spill'` or `'draft_free'`. |
| `vs`, `draft` | m | Vertical shift and draft of the realised design. |
| `stage2` | struct | Stage-2 reference: `vs`, `rho` [kg/m³], `mass` [kg], `Z_CG` (`Final3D.CG_total(3)`, world) [m], `GM` [m], coupled `T_heave`, `T_pitch` [s]. |
| `rho` | struct, kg/m³ | Region densities: `uhpc`, `air` (modular precast); `ballast`, `shell`, `air` (thin shell). |
| `design` | struct | `mode`, `edges` (module edges, body frame) [m], `vs` [m], `t` (shell thickness per module, NaN without void) [m], `z_ballast` (body frame) [m], `solid_modules`. |
| `k_star`, `V_uhpc_target` | index, m³ | Modular precast: ballast module and the UHPC volume of the Stage-2 split, V_i (ρ_i − ρ_air)/(ρ_UHPC − ρ_air). Empty for thin shell. |
| `modules` | struct array | Per module: `z_lo`, `z_hi` [m, body], `t` [m], `h_ballast` (from the module bottom) [m], `V` and `V_<region>` [m³], `mass` [kg], `rho_eff`, `rho_stage2`, `rho_floor` [kg/m³], `CG_world` [m]. |
| `props` | struct | Realised properties (the `final_props` fields). |
| `check` | struct | Comparison with Stage 2: `metrics` (`Z_CG`, `GM`, `T_heave`, `T_pitch`: `value`, `stage2`, `rel_dev`, `limit` = `mass_acceptable_pct`/100, `pass`), `equalities` (`flotation` M/(ρ_w V_sub) − 1 and `GM` GM/GM₂ − 1: `residual`, `tol`, `pass`), `pass` (all metrics and flotation), `failed`, `reason`. |
| `solver` | struct array | Per step run: `step`, `exitflag`, `iterations`, `fval` (Σ((X₃ − X₂)/X₂)² over Z_CG, T_heave, T_pitch), `max_eq_violation`. |
| `fit` | struct array | Fit report of every inner shell surface. |
| `body` | struct | B-rep of the realised body (faces, edges, bodies, shells), body frame, metres. |
| `step_files` | struct array | Written STEP files: `name`, `path`, `bodies`. |

### Field names that changed

Results written before the ballast and material renaming carry the old names. They map as follows.
Results written before `results.stage3` carry the realisation in `results.constructability`
(modular precast) or `results.steel_data` (thin shell), with the names below.

| Old `.mat` field | New field |
| --- | --- |
| `z_fill` (`results.constructability`, `results.steel_data`) | `z_ballast` |
| `rho_fill` in `results.constructability` (air) | `rho_air` |
| `rho_fill` in `results.steel_data` (solid ballast) | `rho_ballast` |
| `rho_steel` in `results.steel_data` (held the ballast density) | removed; use `rho_ballast` |
| `t_steel`, `rho_steel`, `V_steel`, `M_steel`, `z_cg_steel` in `results.constructability` (held UHPC quantities) | `t_uhpc`, `rho_uhpc`, `V_uhpc`, `M_uhpc`, `z_cg_uhpc` |
| `V_fill`, `M_fill`, `z_cg_fill`, `strip_V_fill` in `results.steel_data` | `V_ballast`, `M_ballast`, `z_cg_ballast`, `strip_V_ballast` |
| `rho_fill` in `results.config` (thin shell, solid ballast) | `rho_ballast` |
| `constructability_rho_fill` in `results.config` (air) | `constructability_rho_air` |
| `z_fill`, `V_steel` in `results.constructability.iter_history` | `z_ballast`, `V_uhpc` |
| `in.materials.thin_shell.rho_void`, `in.materials.modular_precast.rho_fill` (inputs) | `in.materials.thin_shell.rho_air`, `in.materials.modular_precast.rho_air` |
| `in.materials.thin_shell.rho_fill` (input) | `in.materials.thin_shell.rho_ballast` |
| `fill_material`, `fill_level` in `results.config.output.style.fill_palette` | `ballast_material`, `ballast_level` |

## `final_props`

`final_props` is the primary exported design record. It uses the base fields listed for `results.Final3D`; realisation fields are added when applicable.

| Field | Type / units | Meaning | Used for |
| --- | --- | --- | --- |
| `fill_method` | char | `'steel_fill'`, `'uhpc_fill'`, or empty for `preliminary`. | Identifying the final design representation. |
| `density_profile_source` | char | `'realised_partition'`: densities of the realised body. | Downstream interpretation. |
| `realised_strip_density` | double, N×1, kg/m³ | Realised effective density (mass/volume) of each module. | Mass-distribution analysis. |
| `realised_strip_edges` | double, (N+1)×1, m | Body-frame module edges. | Geometry reconstruction. |
| `stage3_status` | char | `results.stage3.status`. | Acceptance of the realised design. |
| `stage3_check` | struct | `results.stage3.check`. | Per-metric comparison with Stage 2. |

## Realisation-type availability

| `realisation_type` | `results.stage3` | `final_props` |
| --- | --- | --- |
| `preliminary` | `[]` | Theoretical Stage-2 design. |
| `thin_shell` | Thin-shell realised design | Thin-shell realised design (accepted or failed). |
| `modular_precast` | Modular-precast realised design | Modular-precast realised design (accepted or failed). |
