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
| `constructability` | struct or typed-empty fields | Modular-precast realisation record. | Modular-precast analysis. |
| `steel_data` | struct or `[]` | Thin-shell realisation record. | Thin-shell analysis. |
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

### `results.constructability`

`results.constructability` stores the modular-precast realisation record.

| Fields | Units / convention | Used for |
| --- | --- | --- |
| `t_uhpc`, `t_offset_strip`, `t_UHPC`, `t_min`, `t_min_active` | m. `Inf` thickness denotes a solid strip. | Precast wall and strip geometry. |
| `z_ballast` | m, body frame. Top of the solid ballast region. | Precast ballast level. |
| `wall_strip_idx`, `wall_z_bottom`, `wall_z_top`, `wall_height`, `strip_edges`, `strip_z_lo`, `strip_z_hi` | Indices and body-frame elevations [m]. | Wall and strip layout. |
| `rho_hull`, `rho_uhpc`, `rho_air`, `strip_rho_eff` | kg/m³. `rho_air` is the density of the air in the voids. | Material definition. |
| `V_uhpc`, `V_air`, `V_hull`, `M_uhpc`, `M_air`, `z_cg_uhpc`, `z_cg_air` | m³, kg and m (body frame). | Volume, mass and centroid accounting. |
| `strip_V_total`, `strip_V_UHPC`, `strip_V_void`, `strip_mass_UHPC`, `strip_mass_void`, `strip_mass_total` | m³ and kg. | Per-strip mass accounting. |
| `strip_is_wall`, `strip_is_solid`, `strip_is_feasible`, `is_solid_strip`, `feasibility` | Logical values and summary struct. | Constructability reporting. |
| `strip_Iyy_UHPC`, `strip_Iyy_void`, `strip_z_cg` | kg·m² and m. | Per-strip inertia and centre-of-gravity data. |
| `contours_outer`, `contours_inner` | Cell arrays of `[x,y]` polygon vertices [m]. | Cross-section visualisation. |
| `M_total`, `V_hull`, `V_sub`, `Aw`, `CB_z_world`, `KM_world`, `GM_realised`, `T_heave_realised`, `T_pitch_realised` | Mass, hydrostatic, and period outputs. | Final design reporting. |
| `targets`, `residuals`, `mass_balance_error_pct`, `phi_star`, `feasible`, `exitflag`, `solver`, `elapsed_seconds` | Solver and residual data. | Realisation-solve review. |

### `results.steel_data`

`results.steel_data` stores the thin-shell realisation record.

| Fields | Units / convention | Used for |
| --- | --- | --- |
| `t_steel`, `z_ballast`, `draft`, `vertical_shift`, `draft_optimiser`, `vs_optimiser` | m. `z_ballast` is body-frame elevation. | Thin-shell geometry and Stage-2 comparison. |
| `rho_shell`, `rho_ballast`, `rho_air` | kg/m³. `rho_ballast` is the density of the solid ballast. | Material definition. |
| `V_steel`, `V_shell`, `V_ballast`, `V_air`, `V_hull` | m³. `V_steel` = shell + ballast. | Volume accounting. |
| `M_steel`, `M_shell`, `M_ballast`, `M_air`, `M_total` | kg. `M_steel` = shell + ballast. | Mass accounting. |
| `CG_z_body`, `CG_z_world`, `CB_z_world`, `KM_world`, `GM_realised` | m | Stability reporting. |
| `Ixx_total_origin`, `Iyy_total_origin`, `Izz_total_origin`, `Ixx_about_cg`, `Iyy_about_cg`, `Izz_about_cg` | kg·m² | Inertia reporting. |
| `T_heave_realised`, `T_pitch_realised`, `K33_hydro`, `K55_hydro`, `A11`, `A33`, `A55` | Dynamic properties. | Realised-response reporting. |
| `strip_rho_eff`, `strip_edges`, `strip_V_env`, `strip_V_solid`, `strip_V_void`, `strip_V_ballast`, `strip_V_shell` | Per-strip density, geometry, and volumes. | Strip-level reporting. |
| `targets`, `residuals`, `mass_balance_error_pct`, `phi_star`, `feasible`, `exitflag`, `solver`, `elapsed_seconds` | Solver and residual data. | Realisation-solve review. |

### Field names that changed

Results written before the ballast and material renaming carry the old names. They map as follows.

| Old `.mat` field | New field |
| --- | --- |
| `z_fill` (`results.constructability`, `results.steel_data`) | `z_ballast` |
| `rho_fill` in `results.constructability` (air) | `rho_air` |
| `rho_fill` in `results.steel_data` (solid ballast) | `rho_ballast` |
| `rho_steel` in `results.steel_data` (held the ballast density) | removed; use `rho_ballast` |
| `t_steel`, `rho_steel`, `V_steel`, `M_steel`, `z_cg_steel` in `results.constructability` (held UHPC quantities) | `t_uhpc`, `rho_uhpc`, `V_uhpc`, `M_uhpc`, `z_cg_uhpc` |
| `V_fill`, `M_fill`, `z_cg_fill`, `strip_V_fill` in `results.steel_data`; `V_fill`, `M_fill`, `z_cg_fill` in `final_props` | `V_ballast`, `M_ballast`, `z_cg_ballast`, `strip_V_ballast` |
| `rho_fill` in `results.config` (thin shell, solid ballast) | `rho_ballast` |
| `constructability_rho_fill` in `results.config` (air) | `constructability_rho_air` |
| `z_fill`, `V_steel` in `results.constructability.iter_history` | `z_ballast`, `V_uhpc` |
| `in.materials.thin_shell.rho_void`, `in.materials.modular_precast.rho_fill` (inputs) | `in.materials.thin_shell.rho_air`, `in.materials.modular_precast.rho_air` |
| `in.materials.thin_shell.rho_fill` (input) | `in.materials.thin_shell.rho_ballast` |

## `final_props`

`final_props` is the primary exported design record. It uses the base fields listed for `results.Final3D`; realisation fields are added when applicable.

| Field | Type / units | Meaning | Used for |
| --- | --- | --- | --- |
| `fill_method` | char | `'steel_fill'`, `'uhpc_fill'`, or empty for `preliminary`. | Identifying the final design representation. |
| `density_profile_source` | char | Source of realised strip density data. | Downstream interpretation. |
| `realised_strip_density` | double, N×1, kg/m³ | As-built effective density by strip. | Mass-distribution analysis. |
| `realised_strip_edges` | double, (N+1)×1, m | Body-frame strip boundaries. | Geometry reconstruction. |
| `realised_strips` | struct | Per-strip realised geometry, volumes, masses, and flags. | Detailed fabrication and post-processing use. |

### `final_props.realised_strips`

| Field | Units / convention | Used for |
| --- | --- | --- |
| `fill_method` | char | Identifying the realisation method. |
| `t_offset`, `is_solid`, `is_wall` | m and logical flags | Strip geometry and classification. |
| `z_lo`, `z_hi` | m, body frame | Strip elevations. |
| `V_uhpc`, `V_void`, `mass_uhpc`, `mass_void` | m³ and kg | Material quantities by strip. |
| `uhpc_volume_fraction` | dimensionless | Per-strip UHPC fraction. |
| `contours_outer`, `contours_inner` | Cell arrays of `[x,y]` polygons [m] | Cross-section visualisation. |

## Realisation-type availability

| `realisation_type` | `results.steel_data` | `results.constructability` | `final_props` |
| --- | --- | --- | --- |
| `preliminary` | Empty | Typed-empty fields | Theoretical Stage-2 design. |
| `thin_shell` | Thin-shell record | Typed-empty fields | Thin-shell realised design. |
| `modular_precast` | Empty | Modular-precast record | Modular-precast realised design. |
