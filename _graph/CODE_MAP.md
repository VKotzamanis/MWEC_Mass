# C1_graph code map

Target archive: `MASS_SUITE_C1_COMPLETE_R2026a.zip`
Target SHA-256: `5ee90b4c86c1c19a9c9348e5c65af4f89b923de4a21aa56368797eac903bd079`
Inventory: **197 files**; 158 MATLAB, 6 Markdown, 33 artifacts.

This is a graph-only publication map. Paths are relative to the extracted target archive; no source or data payload is included.

## Directory and file inventory

| path | kind | bytes | functions/classes |
|---|---|---:|---|
| `CHANGELOG.md` | Markdown | 1589 | — |
| `Input/C1.ms2` | artifact | 2547 | — |
| `Input/C1_wamit_cache.mat` | artifact | 227748 | — |
| `Input/WAMIT/WAMIT_Cache_Builder.m` | MATLAB | 9980 | WAMIT_Cache_Builder:1, discover_drafts:60, init_hydro_cache_struct:150, parse_wamit_1_local:77, parse_wamit_3_local:124 |
| `Input/WAMIT/WAMIT_Input_Writer.m` | MATLAB | 7531 | WAMIT_Input_Writer:1, discover_vertical_shifts:68, write_6x6_block:139, write_value_block:131, write_wamit_cfg:115, write_wamit_frc:98, write_wamit_pot:82 |
| `Input/WAMIT/rename_cache_fields.m` | MATLAB | 2079 | rename_cache_fields:1 |
| `Output/C1_modular_precast_results.mat` | artifact | 7797093 | — |
| `Output/C1_preliminary_results.mat` | artifact | 249496 | — |
| `Output/C1_thin_shell_results.mat` | artifact | 262796 | — |
| `Output/modular_precast/WEC_Constructability_Strips.fig` | artifact | 123697 | — |
| `Output/modular_precast/WEC_Constructability_Strips.pdf` | artifact | 91689 | — |
| `Output/modular_precast/WEC_Constructability_Strips.png` | artifact | 226691 | — |
| `Output/modular_precast/WEC_Constructability_XZ.fig` | artifact | 310324 | — |
| `Output/modular_precast/WEC_Constructability_XZ.pdf` | artifact | 158392 | — |
| `Output/modular_precast/WEC_Constructability_XZ.png` | artifact | 234936 | — |
| `Output/modular_precast/WEC_Final_3D_CrossSection_UHPC.fig` | artifact | 51415 | — |
| `Output/modular_precast/WEC_Final_3D_CrossSection_UHPC.pdf` | artifact | 32579 | — |
| `Output/modular_precast/WEC_Final_3D_CrossSection_UHPC.png` | artifact | 190196 | — |
| `Output/modular_precast/diagnostic_panels.log` | artifact | 10618 | — |
| `Output/modular_precast/final_results.log` | artifact | 757 | — |
| `Output/modular_precast/stage2_summary.log` | artifact | 280 | — |
| `Output/preliminary/WEC_Final_3D_CrossSection.fig` | artifact | 44948 | — |
| `Output/preliminary/WEC_Final_3D_CrossSection.pdf` | artifact | 27838 | — |
| `Output/preliminary/WEC_Final_3D_CrossSection.png` | artifact | 165717 | — |
| `Output/preliminary/diagnostic_panels.log` | artifact | 8662 | — |
| `Output/preliminary/final_results.log` | artifact | 757 | — |
| `Output/preliminary/stage2_summary.log` | artifact | 280 | — |
| `Output/thin_shell/Steel_Solve.fig` | artifact | 112930 | — |
| `Output/thin_shell/Steel_Solve.pdf` | artifact | 67235 | — |
| `Output/thin_shell/Steel_Solve.png` | artifact | 202725 | — |
| `Output/thin_shell/WEC_Final_3D_CrossSection_Steel.fig` | artifact | 45015 | — |
| `Output/thin_shell/WEC_Final_3D_CrossSection_Steel.pdf` | artifact | 28139 | — |
| `Output/thin_shell/WEC_Final_3D_CrossSection_Steel.png` | artifact | 165612 | — |
| `Output/thin_shell/diagnostic_panels.log` | artifact | 8662 | — |
| `Output/thin_shell/final_results.log` | artifact | 757 | — |
| `Output/thin_shell/stage2_summary.log` | artifact | 280 | — |
| `Output/thin_shell/steel_solve.log` | artifact | 1379 | — |
| `README.md` | Markdown | 12440 | — |
| `WEC_Output_Options.m` | MATLAB | 8892 | WEC_Output_Options:1 |
| `WEC_User_Input.m` | MATLAB | 10500 | WEC_User_Input:1 |
| `docs/HAMS_MREL_ROUTE.md` | Markdown | 6945 | — |
| `docs/METHODS_ENGINE.md` | Markdown | 50640 | — |
| `docs/RESULT_SCHEMA.md` | Markdown | 58816 | — |
| `docs/RUNTIME_GUIDE.md` | Markdown | 46261 | — |
| `src/+mwecmass/+bem/+hams_mrel/HamsWriter.m` | MATLAB | 9720 | HamsWriter:1, control_file:15, hydrostatic_file:108, matrix_6x6:5, pnl_file:150 |
| `src/+mwecmass/+bem/+hams_mrel/band_averaged_damping.m` | MATLAB | 627 | band_averaged_damping:1 |
| `src/+mwecmass/+bem/+hams_mrel/default_hams_params.m` | MATLAB | 3643 | default_hams_params:1 |
| `src/+mwecmass/+bem/+hams_mrel/generate_draft_nodes.m` | MATLAB | 7047 | generate_draft_nodes:1 |
| `src/+mwecmass/+bem/+hams_mrel/hams_constants.m` | MATLAB | 479 | hams_constants:1 |
| `src/+mwecmass/+bem/+hams_mrel/parse_wamit_1_file.m` | MATLAB | 5088 | parse_wamit_1_file:1 |
| `src/+mwecmass/+bem/+hams_mrel/parse_wamit_3_file.m` | MATLAB | 4602 | parse_wamit_3_file:1 |
| `src/+mwecmass/+bem/+hams_mrel/run.m` | MATLAB | 1884 | run:1 |
| `src/+mwecmass/+bem/+hams_mrel/run_at_draft.m` | MATLAB | 11285 | run_at_draft:1, setup_hams_directory:217, write_panelizer_hull_pnl:246 |
| `src/+mwecmass/+bem/+hams_mrel/run_solver.m` | MATLAB | 2968 | run_solver:1 |
| `src/+mwecmass/+bem/+wamit/restore_mesh_sizing_from_cache.m` | MATLAB | 3048 | restore_mesh_sizing_from_cache:1 |
| `src/+mwecmass/+bem/empty_hydro_cache.m` | MATLAB | 1652 | empty_hydro_cache:1 |
| `src/+mwecmass/+bem/get_or_run_hydro.m` | MATLAB | 2943 | get_or_run_hydro:1 |
| `src/+mwecmass/+bem/interpolate_at_draft.m` | MATLAB | 4901 | interpolate_at_draft:1 |
| `src/+mwecmass/+bem/interpolate_matrix.m` | MATLAB | 1101 | interpolate_matrix:1 |
| `src/+mwecmass/+bem/load_hydro_cache.m` | MATLAB | 4582 | load_hydro_cache:1 |
| `src/+mwecmass/+bem/rebuild_config_hydro.m` | MATLAB | 4488 | rebuild_config_hydro:1 |
| `src/+mwecmass/+bem/retransform_at_cg.m` | MATLAB | 1607 | retransform_at_cg:1 |
| `src/+mwecmass/+bem/transform_to_cg.m` | MATLAB | 2332 | transform_to_cg:1 |
| `src/+mwecmass/+driver/build_config.m` | MATLAB | 61015 | build_config:1 |
| `src/+mwecmass/+driver/build_hydrostatic_tables.m` | MATLAB | 11154 | build_hydrostatic_tables:1 |
| `src/+mwecmass/+driver/build_strip_geometry_tables.m` | MATLAB | 9061 | build_strip_geometry_tables:1 |
| `src/+mwecmass/+driver/load_inputs.m` | MATLAB | 339 | load_inputs:1 |
| `src/+mwecmass/+driver/parse_hull_deck.m` | MATLAB | 3278 | compute_surface_z_range:50, parse_hull_deck:1 |
| `src/+mwecmass/+driver/run.m` | MATLAB | 6388 | run:1 |
| `src/+mwecmass/+geometry/MS2Parser.m` | MATLAB | 76173 | MS2Parser:1, arc_evaluate:1894, bspline_basis_all:1830, bspline_curve_eval:1808, bspline_curve_eval_with_deriv:1868, classify_visible_surfaces:1938, clear_cache:2146, dfs:2101, eval_any_point:803, eval_arc:996, eval_arc_deriv:1055, eval_bcurve:914, eval_bcurve_deriv:1301, eval_bloft_surf:1722, eval_bloft_surf_at_u:1741, eval_bloft_surf_derivs:1570, eval_bsub_curve:978, eval_bsub_curve_deriv:1374, eval_bsub_snake:1181, eval_bsub_snake_deriv:1441, eval_conic:931, eval_conic_deriv:1319, eval_copy_curve:951, eval_copy_curve_deriv:1343, eval_curve:865, eval_curve_or_snake:848, eval_curve_with_deriv:1250, eval_dev_surf:1762, eval_dev_surf_derivs:1627, eval_edge_snake:1148, eval_edge_snake_deriv:1414, eval_line:968, eval_line_deriv:1363, eval_mirr_surf:1773, eval_mirr_surf_derivs:1643, eval_point:759, eval_polycurve2:1008, eval_polycurve2_deriv:1067, eval_proj_curve:1031, eval_proj_curve_deriv:1094, eval_rev_surf:1672, eval_rev_surf_derivs:1505, eval_ruled_surf:1662, eval_ruled_surf_derivs:1492, eval_snake:1119, eval_snake_with_deriv:1391, eval_surface:1201, eval_surface_grid:1235, eval_surface_with_derivs:1458, geometry_summary:2118, get_boundary_entities:2020, get_required_entities:2095, make_clamped_knots:1793, parse:24, parse_abs_bead:543, parse_abs_ring:566, parse_arc:444, parse_bcurve:305, parse_bloft_surf:668, parse_bsub_curve:415, parse_bsub_snake:589, parse_conic:334, parse_copy_curve:365, parse_dev_surf:698, parse_edge_snake:520, parse_entity_line:185, parse_frame_point:247, parse_line_entity:393, parse_mirr_point:281, parse_mirr_surf:720, parse_polycurve2:471, parse_proj_curve:496, parse_rev_surf:641, parse_ruled_surf:619, resolve_mirror_chain:2002, trace_to_edge_snake_of:2053 |
| `src/+mwecmass/+geometry/bspline_param_at_z.m` | MATLAB | 2177 | bspline_param_at_z:1, eval_scalar_bspline:50 |
| `src/+mwecmass/+geometry/chain_open_edges.m` | MATLAB | 4126 | chain_open_edges:1 |
| `src/+mwecmass/+geometry/clip_polygon_at_z.m` | MATLAB | 1689 | clip_polygon_at_z:1 |
| `src/+mwecmass/+geometry/compute_rmin_at_z.m` | MATLAB | 1183 | compute_rmin_at_z:1 |
| `src/+mwecmass/+geometry/detect_open_edges.m` | MATLAB | 2973 | detect_open_edges:1 |
| `src/+mwecmass/+geometry/extract_isocurve_at_z.m` | MATLAB | 1651 | extract_isocurve_at_z:1 |
| `src/+mwecmass/+geometry/extract_midplane_profile.m` | MATLAB | 4722 | extract_midplane_profile:1 |
| `src/+mwecmass/+geometry/find_strip_v_limits.m` | MATLAB | 3153 | find_strip_v_limits:1 |
| `src/+mwecmass/+geometry/find_submerged_v_limits.m` | MATLAB | 3019 | find_submerged_v_limits:1 |
| `src/+mwecmass/+geometry/find_waterline_intersections.m` | MATLAB | 1154 | find_waterline_intersections:1 |
| `src/+mwecmass/+geometry/isocurve_bloftsurf.m` | MATLAB | 2906 | isocurve_bloftsurf:1 |
| `src/+mwecmass/+geometry/isocurve_devsurf.m` | MATLAB | 902 | isocurve_devsurf:1 |
| `src/+mwecmass/+geometry/isocurve_revsurf.m` | MATLAB | 1352 | isocurve_revsurf:1 |
| `src/+mwecmass/+geometry/precompute_boundary_cache.m` | MATLAB | 2640 | precompute_boundary_cache:1 |
| `src/+mwecmass/+geometry/profile_param_at_z.m` | MATLAB | 1202 | profile_param_at_z:1 |
| `src/+mwecmass/+hydrostatics/build_mass_matrix.m` | MATLAB | 1253 | build_mass_matrix:1 |
| `src/+mwecmass/+hydrostatics/cap_contribution.m` | MATLAB | 2220 | cap_contribution:1 |
| `src/+mwecmass/+hydrostatics/compute_hull.m` | MATLAB | 8876 | compute_hull:1 |
| `src/+mwecmass/+hydrostatics/compute_perpendicular_shell_volume.m` | MATLAB | 4296 | compute_perpendicular_shell_volume:1 |
| `src/+mwecmass/+hydrostatics/compute_strip.m` | MATLAB | 9554 | compute_strip:1 |
| `src/+mwecmass/+hydrostatics/compute_submerged.m` | MATLAB | 13593 | compute_submerged:1 |
| `src/+mwecmass/+hydrostatics/compute_wetted_surface_area.m` | MATLAB | 2711 | compute_wetted_surface_area:1 |
| `src/+mwecmass/+hydrostatics/coupled_periods_by_share.m` | MATLAB | 4026 | coupled_periods_by_share:1 |
| `src/+mwecmass/+hydrostatics/hydrostatic_inputs.m` | MATLAB | 4172 | hydrostatic_inputs:1 |
| `src/+mwecmass/+hydrostatics/polygon_properties.m` | MATLAB | 2094 | polygon_properties:1 |
| `src/+mwecmass/+hydrostatics/properties_2d.m` | MATLAB | 18591 | compute_submerged_properties:423, properties_2d:1 |
| `src/+mwecmass/+hydrostatics/properties_3d.m` | MATLAB | 14693 | properties_3d:1 |
| `src/+mwecmass/+hydrostatics/surface_integral.m` | MATLAB | 2806 | — |
| `src/+mwecmass/+hydrostatics/surface_integral_strip.m` | MATLAB | 2166 | surface_integral_strip:1 |
| `src/+mwecmass/+hydrostatics/surface_integral_submerged.m` | MATLAB | 3731 | — |
| `src/+mwecmass/+hydrostatics/waterplane_properties.m` | MATLAB | 1465 | waterplane_properties:1 |
| `src/+mwecmass/+internal/clip_z.m` | MATLAB | 1133 | clip_z:1 |
| `src/+mwecmass/+internal/gauss_legendre.m` | MATLAB | 747 | gauss_legendre:1 |
| `src/+mwecmass/+internal/integrate_piecewise_cubic.m` | MATLAB | 687 | integrate_piecewise_cubic:1 |
| `src/+mwecmass/+internal/offset_polygon.m` | MATLAB | 2069 | offset_polygon:1 |
| `src/+mwecmass/+internal/option_or_config.m` | MATLAB | 399 | option_or_config:1 |
| `src/+mwecmass/+internal/ternary.m` | MATLAB | 199 | ternary:1 |
| `src/+mwecmass/+mesh/close_open_edges.m` | MATLAB | 2490 | — |
| `src/+mwecmass/+mesh/compute_panel_normals.m` | MATLAB | 1109 | compute_panel_normals:1 |
| `src/+mwecmass/+mesh/detect_source_quadrant.m` | MATLAB | 2326 | detect_source_quadrant:1 |
| `src/+mwecmass/+mesh/generate.m` | MATLAB | 32601 | enforce_min_spacing:762, generate:1, grid_to_quads:784, merge_vertices:807 |
| `src/+mwecmass/+mesh/hull_waterline_polygon.m` | MATLAB | 5476 | hull_waterline_polygon:1 |
| `src/+mwecmass/+mesh/panel_grid_counts.m` | MATLAB | 4431 | panel_grid_counts:1 |
| `src/+mwecmass/+mesh/split_panel_at_z.m` | MATLAB | 1975 | split_panel_at_z:1 |
| `src/+mwecmass/+mesh/trim_at_waterline.m` | MATLAB | 1960 | — |
| `src/+mwecmass/+mesh/waterplane_mesh_structured.m` | MATLAB | 5863 | waterplane_mesh_structured:1, wp_arc_len:147, wp_interp:157 |
| `src/+mwecmass/+mesh/waterplane_mesh_unstructured.m` | MATLAB | 8131 | quad_scaled_jacobian:172, resample_polygon:185, waterplane_mesh_unstructured:1 |
| `src/+mwecmass/+optim/PID_Controller.m` | MATLAB | 5124 | PID_Controller:1, PID_Controller:31, reset:124, update:64 |
| `src/+mwecmass/+optim/check_3d_convergence.m` | MATLAB | 2945 | check_3d_convergence:1 |
| `src/+mwecmass/+optim/hams_enrichment_action.m` | MATLAB | 490 | hams_enrichment_action:1 |
| `src/+mwecmass/+optim/range_penalty.m` | MATLAB | 554 | range_penalty:1 |
| `src/+mwecmass/+optim/report_assemble_results.m` | MATLAB | 3525 | report_assemble_results:1 |
| `src/+mwecmass/+optim/run.m` | MATLAB | 19877 | init_history_struct:395, run:1, save_stage2_iteration_local:348, unpack_stage2_errors:433 |
| `src/+mwecmass/+optim/solve_2d_surrogate.m` | MATLAB | 18935 | calculate_convergence_rate:348, calculate_feasibility_score:396, capture_convergence_data:279, constraint_func:183, get_field:455, objective_func:132, range_penalty_2d:465, solve_2d_surrogate:1 |
| `src/+mwecmass/+optim/stage1_oneshot.m` | MATLAB | 3482 | — |
| `src/+mwecmass/+optim/stage1_screen_draft.m` | MATLAB | 2358 | stage1_screen_draft:1 |
| `src/+mwecmass/+optim/stage1_sweep.m` | MATLAB | 6462 | — |
| `src/+mwecmass/+optim/stage1_trained.m` | MATLAB | 14023 | compute_MAPE:289, compute_R2:279, sat_tag:298 |
| `src/+mwecmass/+optim/stage2_bounds.m` | MATLAB | 930 | stage2_bounds:1 |
| `src/+mwecmass/+optim/stage2_constraints.m` | MATLAB | 2852 | stage2_constraints:1 |
| `src/+mwecmass/+optim/stage2_objective.m` | MATLAB | 1580 | stage2_objective:1 |
| `src/+mwecmass/+output/+figures/apply_axes_style.m` | MATLAB | 2182 | apply_axes_style:1 |
| `src/+mwecmass/+output/+figures/apply_layout_style.m` | MATLAB | 501 | apply_layout_style:1 |
| `src/+mwecmass/+output/+figures/build_silhouette_profile.m` | MATLAB | 1610 | build_silhouette_profile:1 |
| `src/+mwecmass/+output/+figures/clip_profile_to_z_range.m` | MATLAB | 521 | clip_profile_to_z_range:1 |
| `src/+mwecmass/+output/+figures/clip_to_halfplane.m` | MATLAB | 1218 | clip_to_halfplane:1 |
| `src/+mwecmass/+output/+figures/compute_inner_profile.m` | MATLAB | 1702 | compute_inner_profile:1 |
| `src/+mwecmass/+output/+figures/compute_strip_equivalent_density.m` | MATLAB | 744 | compute_strip_equivalent_density:1 |
| `src/+mwecmass/+output/+figures/draw_composite_strips.m` | MATLAB | 2548 | draw_composite_strips:1 |
| `src/+mwecmass/+output/+figures/draw_density_strips.m` | MATLAB | 1922 | draw_density_strips:1 |
| `src/+mwecmass/+output/+figures/draw_hatch_strips.m` | MATLAB | 1212 | draw_hatch_local:6, draw_hatch_strips:1 |
| `src/+mwecmass/+output/+figures/export_figure.m` | MATLAB | 1525 | export_figure:1 |
| `src/+mwecmass/+output/+figures/figure_colormap.m` | MATLAB | 754 | figure_colormap:1 |
| `src/+mwecmass/+output/+figures/new_figure.m` | MATLAB | 621 | new_figure:1 |
| `src/+mwecmass/+output/+figures/output_options.m` | MATLAB | 1495 | output_options:1, root_defaults:29 |
| `src/+mwecmass/+output/+figures/plot_complete_convergence.m` | MATLAB | 7013 | plot_complete_convergence:1, plot_stage2_panel:121 |
| `src/+mwecmass/+output/+figures/plot_draft_landscape.m` | MATLAB | 7943 | local_finalize_axes:180, local_legend:168, local_limits:152, plot_draft_landscape:1 |
| `src/+mwecmass/+output/+figures/plot_equivalent_density_2d.m` | MATLAB | 9684 | escape_latex_local:194, plot_equivalent_density_2d:1 |
| `src/+mwecmass/+output/+figures/plot_equivalent_density_3d.m` | MATLAB | 6983 | plot_equivalent_density_3d:1 |
| `src/+mwecmass/+output/+figures/plot_hydrodynamics.m` | MATLAB | 16649 | plot_hydrodynamics:1, shade_T_band:367 |
| `src/+mwecmass/+output/+figures/plot_inner_spline.m` | MATLAB | 1036 | plot_inner_spline:1 |
| `src/+mwecmass/+output/+figures/plot_mesh_diagnostic.m` | MATLAB | 7406 | plot_mesh_diagnostic:1 |
| `src/+mwecmass/+output/+figures/plot_modular_precast.m` | MATLAB | 28617 | build_profile_from_outer_contours:572, compute_inner_offset_local:512, escape_latex_local:503, plot_modular_precast:1 |
| `src/+mwecmass/+output/+figures/plot_optimised_cross_section.m` | MATLAB | 11917 | ordered_plot_profile:183, plot_optimised_cross_section:1 |
| `src/+mwecmass/+output/+figures/plot_panel_normals.m` | MATLAB | 7075 | plot_panel_normals:1 |
| `src/+mwecmass/+output/+figures/plot_steel_solve.m` | MATLAB | 14920 | clip_polygon_above_z:246, clip_polygon_below_z:241, draw_hatch:256, plot_steel_solve:1, save_figure:297 |
| `src/+mwecmass/+output/+figures/presentation_style.m` | MATLAB | 548 | presentation_style:1 |
| `src/+mwecmass/+output/+figures/style_colorbar.m` | MATLAB | 687 | style_colorbar:1 |
| `src/+mwecmass/+output/+figures/style_legend.m` | MATLAB | 362 | style_legend:1 |
| `src/+mwecmass/+output/+figures/style_line.m` | MATLAB | 1284 | style_line:1 |
| `src/+mwecmass/+output/+figures/style_text.m` | MATLAB | 1308 | style_text:1 |
| `src/+mwecmass/+output/Report.m` | MATLAB | 28361 | Report:1, all_panels:419, constructability:207, density_profile:155, exception:354, final_results:382, header:5, matrix:368, stability:113, stage1:27, stage2:75, verdict:278 |
| `src/+mwecmass/+output/check_export_schema.m` | MATLAB | 8169 | check_export_schema:1, local_assert_no_object_leaf:157, local_is_typed_empty:181 |
| `src/+mwecmass/+output/cividis_map.m` | MATLAB | 548 | cividis_map:1 |
| `src/+mwecmass/+output/close_log.m` | MATLAB | 157 | close_log:1 |
| `src/+mwecmass/+output/dispatch.m` | MATLAB | 10322 | dispatch:1, offer:177, report_fids:202, report_stage2_summary:231, save_flag:194, validate_save_flags:211 |
| `src/+mwecmass/+output/emit.m` | MATLAB | 194 | emit:1 |
| `src/+mwecmass/+output/export_results.m` | MATLAB | 5419 | export_results:1, local_assert_no_object_leaf:111, local_config_as_data:59, local_file_sha256:93, mass_typed_empty:137 |
| `src/+mwecmass/+output/export_schema.m` | MATLAB | 16799 | export_schema:1, local_add3:179, local_add:174 |
| `src/+mwecmass/+output/format_axis_publication.m` | MATLAB | 1277 | format_axis_publication:1 |
| `src/+mwecmass/+output/load_results.m` | MATLAB | 805 | load_results:1 |
| `src/+mwecmass/+output/open_log.m` | MATLAB | 888 | open_log:1 |
| `src/+mwecmass/+output/output_dir.m` | MATLAB | 618 | output_dir:1 |
| `src/+mwecmass/+output/print_property_comparison.m` | MATLAB | 4410 | fmt_or_dash:130, print_property_comparison:1, print_row:115 |
| `src/+mwecmass/+output/save_hydro_cache.m` | MATLAB | 156 | save_hydro_cache:1 |
| `src/+mwecmass/+realise/+modular_precast/build_geometry_grid.m` | MATLAB | 4448 | build_geometry_grid:1 |
| `src/+mwecmass/+realise/+modular_precast/build_realised_properties.m` | MATLAB | 974 | build_realised_properties:1 |
| `src/+mwecmass/+realise/+modular_precast/evaluate_design_point.m` | MATLAB | 5233 | evaluate_design_point:1 |
| `src/+mwecmass/+realise/+modular_precast/extract_strip_geometry.m` | MATLAB | 20514 | extract_strip_geometry:1 |
| `src/+mwecmass/+realise/+modular_precast/run.m` | MATLAB | 2304 | run:1 |
| `src/+mwecmass/+realise/+modular_precast/solve.m` | MATLAB | 22469 | con_fn:399, obj_fn:376, out_fn:422, solve:1, unpack_dvs:367 |
| `src/+mwecmass/+realise/+modular_precast/solve_and_extract.m` | MATLAB | 6680 | solve_and_extract:1 |
| `src/+mwecmass/+realise/+preliminary/run.m` | MATLAB | 470 | run:1 |
| `src/+mwecmass/+realise/+thin_shell/build_geometry_grid.m` | MATLAB | 4273 | build_geometry_grid:1 |
| `src/+mwecmass/+realise/+thin_shell/evaluate_design_point.m` | MATLAB | 9069 | evaluate_design_point:1 |
| `src/+mwecmass/+realise/+thin_shell/hull_slope_cos_at_z.m` | MATLAB | 1252 | hull_slope_cos_at_z:1 |
| `src/+mwecmass/+realise/+thin_shell/inner_properties_at_z.m` | MATLAB | 2857 | inner_properties_at_z:1 |
| `src/+mwecmass/+realise/+thin_shell/integrate_split.m` | MATLAB | 4614 | integrate_split:1, tz:92 |
| `src/+mwecmass/+realise/+thin_shell/run.m` | MATLAB | 3698 | close_fid_if_open:63, run:1 |
| `src/+mwecmass/+realise/+thin_shell/solve.m` | MATLAB | 26527 | con_fn:465, invert_piecewise_linear_cumulative:490, obj_fn:444, solve:1 |
| `src/+mwecmass/+realise/+thin_shell/strip_partition_volumes.m` | MATLAB | 5629 | strip_partition_volumes:1 |
| `src/+mwecmass/+realise/build_realised_properties.m` | MATLAB | 14626 | build_realised_properties:1 |
| `src/+mwecmass/+realise/empty_realised_properties.m` | MATLAB | 2108 | empty_realised_properties:1 |
| `validation/diagnostics/hull_at_draft.m` | MATLAB | 9098 | hull_at_draft:1 |
| `validation/diagnostics/reconstruct_live_config.m` | MATLAB | 3030 | local_sha256:37, reconstruct_live_config:1 |
| `validation/diagnostics/stage_animations.m` | MATLAB | 61751 | append_frame:1206, band:1173, clip_poly:1134, crossings:1142, draw_density_section:882, draw_uhpc_from_fraction:608, draw_uhpc_section:1003, get_or:1238, gif1_stage1:46, gif2_stage2:152, gif2b_stage2_uhpc:250, gif2c_stage2_pair:379, gif3_stage3:793, inner_offset:1154, ks_valid:595, new_fig:1200, pad:1233, rho_rgb:1129, stage_animations:1, tick_vector:602, trace_panel:1178 |
| `validation/diagnostics/uhpc_mass_balance.m` | MATLAB | 12254 | uhpc_mass_balance:1 |

## C1 runtime modes

| mode | outputs |
|---|---:|
| `preliminary` | 7 |
| `thin_shell` | 11 |
| `modular_precast` | 13 |

The two C1 inputs are linked to all three runtime modes; each output is linked to its producing mode.

## Edge provenance

- `source_scan`: MATLAB file/function/class declarations.
- `textual`: direct package-local invocation scan; this is evidence, not an AST.
- `textual_handle`: `@name` function-handle reference evidence; never a call assertion.
- `doc_reference`: relative path or public filename references in Markdown.
- `package_mode`: deterministic C1 input/output/runtime-mode relationships.

Graph totals: 540 nodes, 936 edges.
