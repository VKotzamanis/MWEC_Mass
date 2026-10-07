# AGENTS.md — MWEC_Mass working notes for coding agents

Read this file before you touch the code. It records how the code works today, the corrections
the project owner made during review, the defects found, and the behaviour the code must have
when the current work is finished. The implementation plan lives in
`docs/plans/2026-10-06-exact-geometry-stage3-step.md`.

Status: written 2026-10-06, updated 2026-10-07, against commit `8da92d5` (v1.0, `Example_C1`). Work happens on branch
`claude/lucid-cray-7o9442`. Every number quoted below as "measured" was computed in Python from
the v1.0 reference results `Output/C1_modular_precast_results.mat` and
`Output/C1_thin_shell_results.mat`; MATLAB was not available in the cloud container. Floor and
mass figures added on 2026-10-07 are Python estimates on the exact C1 sections.

---

## 1. Ground rules (set by the project owner — follow them without exception)

Terms, as the owner uses them (§3 item 29): the **wall** is the slender upper part of the hull (in
UHPC the solid wall module set by `wall_height`, C1 module 5, 1.8 m; in thin shell the neck the
shell passes through). The **shell** is the layer of thickness t (UHPC or steel) around air. Never
call a shell a wall.

1. **Geometry comes from the parametric definition, always.** Every geometric operation must use
   the `.ms2` parametric surfaces and splines (through `mwecmass.geometry.MS2Parser`) or splines
   fitted to points computed from them. Never assume a shape: no circles, no equivalent radius
   (`sqrt(A/pi)`), no inscribed or centroid-to-boundary radius as a section proxy, no
   axisymmetric or homothetic formula, no `dA/dz ÷ P` slope estimate. Delete any function,
   branch, comment or document sentence that does, in the same change that replaces it.
2. **A `RevSurf` is a spline revolved about an axis, not a circle.** In C1, `surface1` is the
   profile `curve6` (a 0.1 m arc plus half of the quadratic B-spline `curve1`) revolved 90°
   about the vertical line x = 0, y = 1. Its horizontal sections are arcs about that axis whose
   radius r(z) the spline sets. Evaluating it exactly is correct; treating the hull section as a
   circle is not.
3. **Shell thickness is a normal thickness.** The inner surface of the shell is the outer surface
   offset by t along the surface normal n = S_u × S_v / |S_u × S_v|. Build it from offset nodes
   and fit a new spline through them (owner's preferred method). Trim the fold where the offset overlaps itself
   (convex regions with principal curvature radius < t; on C1 only from t = 0.100 m on, at the top
   arc and the neck ends: the smallest convex radius is 0.100 m, so no fold at 25.4 or 76.2 mm).
4. **t_min is a floor, not a fixed value.** Per-module shell thickness t_i ≥ t_min stays a design
   variable. Below the ballast level the section is solid and t_min does not apply there.
5. **No invented tolerances.** Stage-3 acceptance uses `in.pid.mass_acceptable_pct` (10 %).
   The only other gate is the spline-fit band ε = 0.01 · t_min (one tenth of that acceptance
   band, so fitting can never use up more than a tenth of it; see §5 item 9). Any further numeric
   gate needs the owner's approval first.
6. **No side quests.** Work only on correctness of the code. Do not add before/after studies,
   paper comparisons or unrelated refactors.
7. **Thin shell is not UHPC.** The thin-shell (steel) realisation has its own logic. Rebuild its
   geometry on the exact kernel, but do not copy the modular-precast process into it.
8. **Testing runs in the cloud container** with GNU Octave (and Python `gmsh` only to verify
   written STEP files). Commit and push after every task; the container is ephemeral.
9. Do not push the `STEP_Producer` zip the owner uploaded. It stays outside the repository.
10. **Names say what the material is.** The solid region at the bottom is the *ballast*:
    `z_fill` becomes `z_ballast` everywhere (code, `.mat` fields, docs). Air is `rho_air` in both modes;
    the thin-shell solid ballast density is `rho_ballast`; UHPC quantities carry `uhpc` names. No
    variable may carry the name of one material and the value of another (rename table in the plan,
    task T0b; old `.mat` field names map to the new ones in `docs/RESULT_SCHEMA.md`).
11. **Count every volume once.** Each height z in a module belongs to exactly one description:
    below the ballast level the full outer section is solid (no separate shell term); above it,
    shell + air. A 2D shell surface never overlaps a solid body. The wall module (UHPC) and the
    ballast zone are full solid sections; no t_min shell is added to them. Tests assert volume
    closure to machine precision: per module V_UHPC + V_air = V_module (thin shell: V_ballast +
    V_shell + V_air = V_hull), and the module volumes sum to the hull volume.
12. **Judge by metrics.** The owner: "you and any agent judges from metrics." Every judgement (the
    orchestrator's, an implementer's, a grader's, and every report to the owner) rests on numbers
    computed from the code's outputs: exact checks, measured values, the gates of rule 5. Never
    judge by how a figure or a text looks. A figure is judged by its section data: every drawn
    boundary has a face of the solid behind it, the drawn areas equal the body's section areas, and
    the drawn ballast level equals `z_ballast`. The owner reviews real XZ (y = 0) and YZ (x = 0)
    sections of the realised C1 bodies, with these metrics printed beside them.

---

## 2. What the code does today

### 2.1 Pipeline

`WEC_User_Input.m` → `mwecmass.driver.run` → `build_config` → `mwecmass.optim.run`
(Stage 1, Stage 2) → `mwecmass.realise.<type>.run` (Stage 3) → `mwecmass.output.dispatch`
(`.mat`, logs, figures). `realisation_type` is `preliminary`, `thin_shell` or `modular_precast`.

### 2.2 Geometry and hydrostatic tables

- `MS2Parser` reads the MultiSurf deck and evaluates curves and surfaces exactly
  (`eval_surface(name,u,v)`; supported: `RuledSurf`, `RevSurf`, `BLoftSurf`, `DevSurf`, `MirrSurf`).
- C1 (`Input/C1.ms2`) is two patches mirrored in x and y: `surface1` (RevSurf, the rounded ends)
  and `surface2` (RuledSurf from an edge of `surface1` to its projection on y = 0, the flat
  sides). Every horizontal section is a stadium: arcs of radius x_w(z) centred at (0, ±1 m)
  joined by straight sides (measured deviation 7e-16 m over all stored contours).
- `extract_isocurve_at_z` slices the surfaces at a height and returns contour points. The
  driver builds tables `A_w(z)`, `I_wp(z)`, `P(z)`, `V_sub(z)`, `CB_z(z)` from these slices
  (`build_hydrostatic_tables.m`). Stages 1–2 interpolate these tables linearly.

### 2.3 Stage 1 and Stage 2

- Design vector `x = [vertical_shift, rho_1 … rho_N]`: one equivalent (bulk) density per module
  ("strip"). C1 has 5 modules; in modular precast, module 5 is the solid UHPC wall module, pinned to
  2500 kg/m³ (thin shell has no wall module).
- Stage-2 objective: range penalties on GM (to `gm_target`), heave and pitch periods.
  Constraints: flotation equality, GM ≥ `gm_min`, adjacent density ratio, bottom-heavy
  monotonicity (`c_mono`, `optim/stage2_constraints.m`), and for modular precast a minimum
  constructable mass (`c_mass_min`, same file). The Stage-1 2-D surrogate solver
  (`optim/solve_2d_surrogate.m`) holds a Stage-1 copy of each: `c_monotonic` and `c_mass_min`.
- For modular precast, each module's density has a lower bound ("floor") meant to equal the
  density of that module built with a t_min UHPC shell and an air void
  (`build_config.m` ≈ 355–556 → `config.per_strip_density_lb`). **This floor is computed on a
  circle** (see issue I2).
- C1 Stage-2 result (`results.Final3D`): vs = 0.9828 m, ρ = [2500, 1430, 369.4, 369.4, 2500],
  M = 20 387.7 kg, z_CG = −1.040 m (world), GM = 0.200 m, T_pitch = 5.03 s.

### 2.4 Stage 3, modular precast (UHPC) — current behaviour

- `solve_and_extract` → `solve.m` runs a new fmincon optimisation with variables
  `[vs, z_ballast, t_2 … t_N_nonwall]`, objective = heave and pitch period penalties (evaluated
  with uncoupled periods, I20), equality = flotation, inequality = GM ≥ `gm_min`. **It reads only the Stage-2 vertical shift** (and the
  mass for a pre-check); the Stage-2 densities are never used.
- Pre-check (`solve_and_extract.m:64–106`): the Stage-2 mass must lie in [0.95 · M_min, 1.05 · M_max]
  of a UHPC hull built with t_min shells and the wall module solid (M_min) and of a fully solid hull
  (M_max); otherwise the function
  stops with an error (`MassTooLight`, `MassTooHeavy`).
- Inner (void) geometry: each horizontal slice is offset in its own plane by t·L, with
  L = 1/cos α estimated from `dA/dz ÷ P` and capped at 5 (`inner_properties_at_z`,
  `hull_slope_cos_at_z`). One global `z_ballast`: everything below is solid UHPC.
- `extract_strip_geometry` repeats the slicing per module to store `contours_outer/inner` and
  per-module volumes. For C1 it used the `R_eq = sqrt(A/π)` circle fallback, because `P_table`
  (92 entries) and `Aw_table_z` (96) differ in length.
- `build_realised_properties` overwrites `final_props` with the realised design when the solve
  returns finite numbers, otherwise keeps Stage 2.
- C1 result: vs = 0.8932 m, M = 20 660.7 kg, z_CG = −1.182 m, GM = 0.231 m, t = 76.2 mm in every
  hollow module (all at the floor), z_ballast = −1.291 m (world, inside module 2), exitflag 2.

### 2.5 Stage 3, thin shell (steel) — current behaviour and how it differs

- Stage 2 for thin shell: 5 modules with half-height end modules (C1 edges, body frame:
  −3.25, −2.706, −1.619, −0.531, 0.556, 1.10 m); **no wall strip; no density floor**. Densities
  are bounded only by `in.bounds.ballast_density_bounds` = [20, 2500] kg/m³ (the floors in
  `build_config.m` run only when `enable_constructability`, i.e. modular precast).
- `in.materials.thin_shell.t_min` = 25 mm and `t_init` = 20 mm. Both are used **only** in Stage 3
  (lower bound and start of the thickness variable). No 15 mm value exists in the code.
- `thin_shell/solve.m`: fmincon over `[vs, t, z_ballast]` — the draft, **one** uniform shell
  thickness for the whole hull, and one ballast level free to rise through any module.
  Objective: heave and pitch period penalties against the configured goals. Equality:
  flotation. Inequality: GM ≥ `gm_min`. Reads only the Stage-2 vertical shift.
- Materials: shell `rho_shell` (7500), ballast `rho_ballast` (defaults to `rho_shell`), air `rho_air`
  (1.2). Below `z_ballast` the **full** outer section is ballast; above it, a steel annulus encloses air.
- Same in-plane `t/cos α` inner offset as UHPC, plus a forced `cos α = 0.1` (offset 10 t) just
  below the top of the hull. `t_max = 0.5 × compute_rmin_at_z` came out as **0.5025 m** in the
  C1 run: the function measures centroid-to-*vertex* distances, and the neck's flat sides carry
  no vertices, so it returned ≈1.005 m instead of the neck half-width 0.10 m. `compute_rmin_at_z`
  is probed at 9 heights from 10 % to 90 % of the hull height (`thin_shell/solve.m:76–81`); the
  thickness bound is 0.95 · t_max (lines 123, 237), the start is clamped to 0.9 · t_max (line 110),
  and the flag `t_min_active` is set when t lies within 1 % of the t_min to t_max span (line 260).
- If the solve is infeasible, `final_props` keeps the Stage-2 properties.
- C1 thin-shell example: Stage-2 densities [2500, 1297, 310.5, 310.5, 310.5]; Stage 3 built
  t = 27.3 mm, `z_ballast` = −2.731 m (body, inside module 1), realised densities
  [6977, 286, 400, 2153, 2275]. A 15 mm steel shell alone already gives modules 4–5 at least
  ≈1196 / 1336 kg/m³ (≈1978 / 2184 at 25 mm; Python estimate on the exact C1 geometry), so the
  Stage-2 values 310.5 cannot be built.
- Differences from UHPC: one uniform thickness instead of one per module, no wall strip, ballast
  free to pass module edges, steel ballast up to 7500 kg/m³, plate modelled as a surface for
  structural use.

### 2.6 Outputs

`Output/<hull>_<type>_results.mat` holds `results` and `final_props` (see `docs/RESULT_SCHEMA.md`).
Figures: `WEC_Constructability_XZ`, `WEC_Constructability_Strips` (precast), `Steel_Solve`
(thin shell), `WEC_Final_3D_CrossSection_*`. No STEP output exists in the repository today.

---

## 3. Corrections the project owner made during review

1. A `RevSurf` is a revolved spline, not a circle (§1 rule 2). Plot of C1 `surface1` alone was
   produced to confirm this.
2. The current in-plane slope offset fails at the near-horizontal shoulder; the normal must come
   from the parametric surface.
3. Per-module thickness stays a variable with t ≥ t_min; ballast removes the shell only where it is.
4. Stage 3 must start from the **whole** Stage-2 solution (mass, draft, densities, GM, periods),
   not a fresh optimisation:
   1. Split each module's Stage-2 density into V_UHPC and V_air:
      V_UHPC,i = V_i (ρ_i − ρ_air)/(ρ_UHPC − ρ_air).
   2. Build the real shells and ballast that hold V_UHPC,i in each module.
   3. Compute mass, z_CG, Iyy, GM and periods **on the actual 3D geometry** (the same solid that
      is written to STEP) and compare them with Stage 2.
   4. If z_CG, GM and periods are all within `mass_acceptable_pct` of Stage 2, accept
      (flotation itself: see OD10).
   5. Otherwise optimise Stage 3 from that initial split. It may change draft (for mass balance)
      and the mass distribution to come closer to the Stage-2 solution.
5. GM is an **equality**: GM_realised = GM_Stage2 (not GM ≥ `gm_min`), so the optimiser cannot
   default to ballast everywhere.
6. Ballast first, inside one module only. The ballast may enter the next module only if flotation
   and the GM equality cannot be met otherwise. When it does, the owner asked for the more
   robust of (i) rebuilding the module edges or (ii) measuring the ballast level from that module's own
   bottom; the agent chose (ii) (reasons under OD6).
7. If Stage 3 cannot meet its requirements, return the closest realisable design, flag it as
   failed, report which metric failed and why, **plot it, store it in `final_props` and the
   `.mat`, and export its STEP files** — never fall back to Stage 2.
8. Use the PID's `mass_acceptable_pct` as the Stage-3 acceptance tolerance; invent no others.
9. Correct floors will change Stage 2 results; that is expected. Use the correct ρ_air.
10. Delete every circle/shape assumption across the code, including later UHPC and density code.
11. Rebuild the thin-shell mode on the exact kernel without copying the UHPC logic.
12. Stage 3 writes STEP files: UHPC — one per module plus one fused solid of all modules; steel —
    the ballast solid, the shell as a 2D surface, and (if possible) one file with both.
13. Rename `z_fill` → `z_ballast` and remove the density naming trap, in code and docs (done in T0b;
    the old names remain only where v1.0 data are read or recorded: `tools/`, `tests/baseline/*.json`,
    `Output/`, `docs/plans/`, and the rename statements in this file).
14. Thin shell is a separate pipeline:
    1. Stage-2 densities start from (are floored by) a thin-shell minimum shell thickness that the
       user sets in the input file (value fixed later at 25.4 mm, item 19).
    2. Stage 2 runs again for this mode, without a solid wall module.
    3. Stage 3 finds **one uniform thickness** for all modules — built with care in the slender
       neck — and `z_ballast`.
    4. `z_ballast` may rise past module 1 into other modules without penalty; it is not a last
       resort as in UHPC.
15. Stage 3 must evaluate the **coupled** periods, as Stage 2 and the report do. Its objective
    is item 27.
16. The centre-of-gravity check uses the Stage-2 value as the reference:
    |Z_CG,3 − Z_CG,2| / |Z_CG,2| ≤ `mass_acceptable_pct`/100, with Z_CG = `CG_total(3)`, the
    CG elevation in the world frame the code already reports (z = 0 at the still-water line).
17. Thin shell also uses the closest-fail rule: never fall back to Stage 2.
18. Spline fitting is adaptive: the algorithm judges explicit metrics and refines locally where
    curvature, slender sections or folds need more nodes (§5 item 9).
19. Thin shell: t_min = t_init = 25.4 mm (one inch), one value set by the user; it also sets the
    Stage-2 floors.
20. The draft is a Stage-3 variable in both modes but the **last resort**: first find the best
    design at the Stage-2 draft; change the draft only if equilibrium (mass = displaced mass)
    cannot be reached there. UHPC order: ballast within k* → ballast spill → draft → closest fail.
    Thin shell order: (t, z_ballast) at fixed draft → draft → closest fail.
21. Coupled periods are fine to use (owner). Evidence for keeping them: with no surge stiffness
    the body surges freely as it pitches, and the surge–pitch added mass A15 lowers the effective
    pitch inertia by A15²/(M + A11). Heave is uncoupled (A13 = A35 = 0, fore–aft symmetry).
22. CG means the centre of gravity, Z_CG = `CG_total(3)`; compare it against Z_CG of Stage 2
    (item 16). Take the owner's wording literally; ask instead of reinterpreting.
23. The steel shell STEP surface is the exterior parametric surface, from `z_ballast` up only:
    below it the ballast solid already contains the plate (no double counting).
24. Thin-shell Stage-2 upper density bound = 7500 kg/m³ (solid steel); keep ρ_air correct.
25. GM = GM_Stage2 in both modes, for consistency.
26. Delete the Stage-2 monotonic density constraint. Owner: "I see what you're saying. The issue I
    had in the past is that the optimizer would not 'think' to add material first at the lower
    strips and it would provide me un-optimized results. Delete c_mono." Stability stays enforced
    by GM ≥ `gm_min` in Stage 2 and GM = GM_Stage2 in Stage 3. To answer the owner's concern,
    Stage 2 also runs from a bottom-filled start at the Stage-1 draft (agent's proposal, owner
    informed): every module at its floor, then modules filled to the mode's solid density from the
    keel up until mass = ρ_w · V_sub; the UHPC wall module stays pinned. Stage 2 keeps the start
    with the lower objective and logs both starts. Evidence: I22.
27. Stage-3 objective, both modes (owner's consistency rule, OD11). Owner: "For now make it come
    closer to Stage 2 converged solution." Stage 3 minimises Σ ((X3 − X2)/X2)² over
    X ∈ {Z_CG = `CG_total(3)` (world frame), coupled T_heave, coupled T_pitch}, with flotation and
    GM = GM_Stage2 as equalities. After every solve the `mass_acceptable_pct` check on Z_CG, GM,
    T_heave and T_pitch decides accepted or failed. This replaces the heave and pitch range
    penalties for Stage 3 (items 15 and 16 keep the coupled periods and the Z_CG reference).
    Stage-1 and Stage-2 objectives do not change. At the Stage-2 draft the equalities fix Z_CG and
    T_heave, so the objective acts through T_pitch only. Closest fail: if the equalities hold but
    the check fails, the optimum is the closest design (it minimises the deviation); if the
    equalities cannot be met, the iterate with the smallest equality violation is the closest
    design.
28. UHPC shell above the ballast (agent's decision under rule 4, owner informed; owner on the
    wording: "What do you mean? The wall is supposed to be solid in UHPC. What wall are you
    referring to even? ... make sure we dont double count on the wall section the solid wal +
    t_min."). In the ballast module k*, the UHPC shell around the air above the ballast has
    thickness t_k* ≥ t_min and is a Stage-3 variable, like the shells of the hollow modules above
    k*. The initial split starts it at t_min. The wall module and the ballast zone are full solid
    sections with no t_min shell added (rule 11).
29. Terms, as the owner uses them: the **wall** is the slender upper part of the hull (UHPC: the
    solid wall module, C1 module 5; thin shell: the neck); the **shell** is the layer of thickness
    t around air. Never call a shell a wall.
30. Thin-shell `t_max`. Owner: "to make the upper limit of t be the thickness of the wall/2
    (because offset from both edges). I agree." `t_max` is half the thickness of the slender wall
    (the neck), where the two offset shells meet, from the exact geometry: 0.10 m for C1.
31. Draft. Owner (earlier): "a change in draft, for both cases, should be the last case option
    (i.e., if we have created the best possible solution but the mass balance isn't satisfied
    within tolerance." The draft is released only when mass balance cannot be met at the Stage-2
    draft (UHPC: after the spill step). A failed `mass_acceptable_pct` check at the Stage-2 draft
    does not release it; that design ends as the closest fail.
32. Mass balance (resolves OD10). Owner: "check if mass balance is satisfies and if the Z_CG, GM
    and periods are the same... within the acceptable bound". Flotation is an equality held to
    the solver's constraint tolerance in every reported state. `mass_acceptable_pct` applies only
    to Z_CG, GM, T_heave and T_pitch.
33. Delete the Stage-2 minimum-mass constraint `c_mass_min` (agent's decision, owner informed):
    `m_min_constructability` = Σ V_i · `per_strip_density_lb`(i) (`build_config.m:579–609`) is
    implied by the density bounds, and the constraint was inactive at the C1 optimum (−1.779).
    The m_min computation has one other use, the feasibility check at `build_config.m:620–645`,
    so it stays for that check.
34. Delete the UHPC Stage-3 pre-check (agent's decision, owner informed). It uses invented factors
    (0.95, 1.05), aborts instead of applying the closest-fail rule (item 7), and is redundant once
    the Stage-2 floors from the kernel exist. Evidence: I23.
35. Run time (owner, 2026-10-07). In Octave, `build_config` alone takes 34.5 min for C1, mostly
    name lookups of parser entities. Owner: "Then make it SAVE that and reload it" and, on
    resolving each curve and surface once instead of on every point, "yes definetely do that."
    Plan tasks T0c (parser) and T0d (geometry cache); both leave every number identical.
36. General hulls (owner, 2026-10-07). On the proposal to build the geometry kernel for hulls like
    C1 only and stop every other hull with a named error, the owner: "No it needs to be
    generalized." The kernel takes any entity type `MS2Parser` evaluates and any parametrisation
    for hulls whose horizontal section is one closed loop at every height (one body, no holes;
    at a horizontal shelf, the loops just below and above bound it); a section of several loops
    raises `mwecmass:solid:SectionNotClosed`, and a hull the parser would misread (it references an
    entity of a type the parser skips, or has a mirror plane other than x = 0 or y = 0) stops with a
    named error; deck lines the hull does not reference are ignored. Each outer patch takes the
    exact path (its entity converts exactly to NURBS and z depends on one parameter (F1 orders it as
    u), monotone: all of C1) or the general path: faces fitted through exact points of the parametric definition with z as
    one parameter, split where z turns back and at creases, so every face is an untrimmed patch and
    every horizontal cut is a parameter line. Inner faces whose offset keeps no structure of the
    outer patch are fitted the same way. Metrics M1–M3 are judged between the faces as written;
    fitted outer faces stay within ε/4 of the exact surface (§5 item 9). Interface contract §8.

---

## 4. Issues found (evidence for C1 unless stated)

| # | Issue | Where | Evidence |
|---|---|---|---|
| I1 | Shell built by in-plane offset t/cos α with cap L ≤ 5; fails near horizontal surfaces and next to them | `thin_shell/inner_properties_at_z.m`, `hull_slope_cos_at_z.m`, `modular_precast/extract_strip_geometry.m:138–218` | Module 4: 41 of 100 sections below 76.2 mm, minimum 18 mm; up to 155 mm above the shoulder |
| I2 | Stage-2 density floor computed on a circle of radius `compute_rmin_at_z` | `hydrostatics/compute_perpendicular_shell_volume.m`, `build_config.m` 405–556 | Floors 437/227/175/369 kg/m³ reproduced from the formula; buildable minimum with a true 76.2 mm shell ≈ 608/284/211/610 |
| I3 | Stage-2 optimum is not buildable, so Stage 3 cannot reproduce it | consequence of I2 | Module 4 placed at 369.4 kg/m³, below its buildable minimum |
| I4 | Stage 3 ignores the Stage-2 densities and re-optimises from scratch | `modular_precast/solve.m:70` | Realised design differs: +273 kg, draft +9 cm, z_CG −142 mm |
| I5 | Equivalent-radius `sqrt(A/π)` slope fallback used for C1's stored contours | `extract_strip_geometry.m:146–149` | `P_table` 92 vs `Aw_table_z` 96 entries |
| I6 | Solver (table areas + contour offset) and stored contours (separate slicing) are two geometries | `build_geometry_grid.m` vs `extract_strip_geometry.m` | Stored contours: +0.85 % mass, z_CG +7 mm vs solver |
| I7 | `t_max = 0.5 × compute_rmin_at_z` (a radius) bounds the thickness | both `solve.m` files | — |
| I8 | Figures do not show the solved geometry: voids squeezed to a volume ratio, `z_ballast` ignored, shell lines from 2D offsets or t_min | `plot_modular_precast.m`, `compute_inner_profile.m`, `plot_inner_spline.m`, `plot_steel_solve.m`, `stage_animations.m` | Read as a 3D body, the XZ figure weighs 23 687 kg with z_CG 71 mm higher |
| I9 | Zero-thickness points: the offset leaves seam vertices unmoved | `internal/offset_polygon.m:47–50` | 932 inner vertices lie on the outer hull |
| I10 | `strip_scale_factor`, `strip_r_min`, `feasibility.s_max` assume shape scaling | `extract_strip_geometry.m`, `export_schema.m`, `Report.m` | — |
| I11 | Name trap (fixed by T0b): one density name meant air in one mode and solid ballast in the other | `WEC_User_Input.m`, `build_config.m` | Values were correct; the name invited a wrong density |
| I12 | `BM_L` taken about x = 0 without the parallel-axis term | `hydrostatics/properties_3d.m:95–100` | Exact for C1 (symmetric); wrong for a hull whose waterplane centroid is off x = 0 |
| I13 | Silhouette for figures = widest x mirrored, labelled "midplane"; midplane points sampled with a 5 cm band | `build_silhouette_profile.m`, `extract_midplane_profile.m` | — |
| I14 | Dead `config.shell` branches (always empty) | `compute_strip_equivalent_density.m`, figure and `properties_2d` branches | `build_config.m:283` is the only assignment |
| I15 | No STEP output in the code | — | Owner requirement |
| I16 | Thin-shell Stage 2 has no density floor, so it can request unbuildable densities | `build_config.m:358` (floors only for precast) | C1: modules 4–5 at 310.5 kg/m³; a 15 mm shell alone gives ≥ ≈1196 / 1336 |
| I17 | `compute_rmin_at_z` measures distances to vertices only; flat faces have none | `geometry/compute_rmin_at_z.m:36–37` | C1 thin shell: `t_max` = 0.5025 m although the neck closes at t = 0.10 m |
| I18 | Stage-2 density upper bound is 2500 kg/m³ in every mode, although steel ballast reaches 7500 | `WEC_User_Input.m` `ballast_density_bounds`, `stage2_bounds.m` | C1 thin shell: Stage 2 capped module 1 at 2500; Stage 3 built 6977 |
| I19 | More material-name traps (fixed by T0b): the UHPC path used steel names for UHPC density, volume and thickness; thin shell stored the ballast density under a steel name | `modular_precast/solve_and_extract.m:18`, `solve.m`, `thin_shell/solve.m` packaging | — |
| I20 | Stage 3 optimises **uncoupled** periods, while Stage 2 and the reported results use coupled periods | `modular_precast/evaluate_design_point.m:119–126`, `thin_shell/evaluate_design_point.m:200–207` vs `properties_3d.m:254–257`, `build_realised_properties.m` | Pitch only (heave identical): precast 4.960 s optimised vs 4.907 s reported; thin shell 3.890 s optimised (the target) vs 3.813 s reported. A15²/(M+A11) = 2.1 % and 3.9 % of the pitch inertia |
| I21 | C1 floats with ≈95 % of its volume submerged (V_sub 20.16 of 21.16 m³). Mass balance therefore cannot carry a percentage tolerance | hydrostatics of C1 | From the C1 tables: +1 % mass raises the waterline 160 mm; +5 % submerges the hull completely; −10 % lowers the waterline 256 mm |
| I22 | The monotonic density constraint `c_mono` makes Stage 2 infeasible for thin shell with the true floors, and over-constrains UHPC | `optim/stage2_constraints.m:35,39,51,64`; Stage-1 copy `c_monotonic` in `optim/solve_2d_surrogate.m:224–238,248,261–264,269` | Python estimate on the exact C1 sections. Thin shell (25.4 mm steel shell, air inside): module floors ≈ 722/268/451/2009/2217 kg/m³, so `c_mono` forces every module to ≥ 2217 kg/m³: minimum mass 46 934 kg against 21 702 kg displaced with the hull fully submerged; no feasible point. UHPC (76.2 mm): floors ≈ 608/284/211/610 (wall module 2500); `c_mono` lifts modules 1–3 to ≥ 610 kg/m³, ≈ 3.1 t more in module 3 than its floor requires |
| I23 | The UHPC Stage-3 pre-check uses invented factors and aborts instead of returning the closest design | `modular_precast/solve_and_extract.m:64–106` (factors 0.95 at line 95, 1.05 at line 101; errors `MassTooLight`, `MassTooHeavy`) | Rule 5 allows no further numeric gate; §3 item 7 requires the closest-fail rule; the check is redundant once the Stage-2 floors from the kernel exist |

Known approximations **not** in scope (report, do not change without approval): linear
interpolation of the hydrostatic tables in Stages 1–2; angular sorting of contour points about
their mean (assumes star-shaped sections); CG_x = CG_y = 0 (declared symmetric-body model).

---

## 5. Expected outcome

1. **One geometry kernel.** Outer surfaces come from `MS2Parser`; inner surfaces of the shells
   from the normal offset of those surfaces (offset nodes → fold trimming → fitted spline surface). Slices,
   volumes, centroids and inertias are computed from these surfaces. Stage-2 floors, Stage-3
   realisations, stored contours, figures and STEP files all use this kernel.
2. **Stage-2 floors** equal the true density of each module with a t_min normal shell and air.
3. **Stage 3, modular precast** follows §3 item 4: split → build → check against Stage 2 →
   optimise if needed (variables: ballast level in the ballast module k*, t_k* and the t_i of the
   hollow modules above k*, and the draft as the last resort; equalities: flotation and
   GM = GM_Stage2; objective: item 11) → spill into the next module only if needed →
   closest-fail report if the equalities cannot be met or the check fails.
4. **Thin shell** keeps its own pipeline: Stage-2 floors from the user-set minimum shell
   thickness (25.4 mm) on the exact geometry, no wall strip, upper bound 7500 kg/m³; Stage 3
   solves one uniform thickness and `z_ballast` (free to pass module edges), with the draft as
   the last resort, flotation and GM = GM_Stage2 as equalities, objective item 11, the
   exact-normal inner surface built correctly in the slender neck, and an exact-geometry `t_max`:
   half the thickness of the neck, where the two offset shells meet (0.10 m for C1). The probe
   range, the 0.95 · t_max and 0.9 · t_max limits and the 1 % `t_min_active` flag of today's
   `thin_shell/solve.m` are deleted.
5. `final_props` and the `.mat` always describe the realised (or closest-fail) 3D design, with a
   status flag and a per-metric report.
6. **Figures** draw true sections of the realised solid and honour the ballast level.
7. **STEP files** written by Stage 3 (METRE units, named solids):
   - UHPC: `<hull>_UHPC_module_<i>.step` for each module, and `<hull>_UHPC_all.step` — all
     modules as one connected solid (built directly as one B-rep, no boolean needed).
   - Steel: `<hull>_STEEL_ballast.step` (solid: full outer section below `z_ballast`),
     `<hull>_STEEL_shell.step` (2D surface: the exterior parametric surface from `z_ballast` to the
     deck), and a combined file with both bodies sharing the junction curve at `z_ballast`.
   - Verified on import with `gmsh`/OpenCASCADE: closed solids, correct solid count, volume equal
     to the kernel's volume.
8. Docs (`METHODS_ENGINE`, `RUNTIME_GUIDE`, `RESULT_SCHEMA`) describe the new methods only.
9. **Adaptive spline fitting with judged metrics** (error-bounded fitting with local knot
   insertion and knot removal, after Piegl & Tiller, *The NURBS Book*, 2nd ed., 1997, ch. 5 and 9):
   1. Place initial nodes densely where the outer surface curves sharply or the void is narrow.
   2. Offset the nodes by t + ε/2 along the exact normal, trim the fold, and **split the face
      along any crease** the trimming leaves, and along every C0 seam of the outer surface,
      instead of forcing one smooth spline across it — a smooth spline across a crease overshoots.
      A concave crease gets a face of its own: the crease curve offset by t + ε/2 along the fan of
      normals of its two sides. Where concave creases meet at a vertex whose cone of normals spans
      a solid angle (e.g. the foot of a re-entrant edge of an L-, T- or cross-shaped column on a
      wider pontoon, where three concave creases meet), the gap between their fan faces is closed by
      a face of its own: the vertex offset by t + ε/2 along every normal of its cone (part of a
      sphere about the vertex), fitted with z as one parameter, its boundaries shared bitwise with
      the adjacent fan faces.
   3. Fit cubic B-splines through the nodes; check them on dense points *between* the nodes.
   4. Metrics the algorithm judges at every check point, between the faces as written (outer as
      written, inner as fitted):
      - **M1 (hard):** local normal thickness t_local ≥ t_min.
      - **M2 (band):** t ≤ t_local ≤ t + ε, with ε = 0.01 · t_min. Fitting to t + ε/2 centres
        the error band, so the floor holds while the error stays below ε.
      - **M3 (hard):** the inner surface stays inside the outer surface, never self-intersects,
        and every horizontal slice of a void is one simple closed loop; opposite shell layers of a
        slender section never meet.
   5. Insert knots only in the spans that fail; refit and recheck. Stop with an error report,
      never silently, if the iteration cap is reached.
   6. Remove knots that are not needed while M1–M3 still hold, so the STEP stays lean.
   7. Mass properties are computed on the fitted surfaces themselves (the STEP geometry), so any
      fit deviation is already in the reported mass, CG and inertia.
   Outer faces are exported as exact NURBS conversions of the `.ms2` entities where the entity
   type allows and z depends on one parameter (F1 orders it as u), monotone (all of C1: B-spline curves,
   arcs, revolution, ruled surface); otherwise they are fitted with the same metrics through exact
   points of the parametric definition, with z as one parameter of every face, split where z turns
   back and at creases, so every face stays an untrimmed patch (§3 item 36). Where a section lies
   wholly on one patch with no seam or crease point (a smooth dome, a crowned deck, a revolution
   about a non-vertical axis), the face is cut along one z-monotone curve on the exact surface (the
   steepest-ascent line of z), used bitwise as both of its v-boundaries, and a band that ends at a
   single highest or lowest point inside a patch ends in a pole row there. Bands end at one set of
   heights for the whole hull, and a point where a seam, crease or cut ends inside another face's
   row is joined by a cut to a vertex of that face's opposite row, so every face boundary has one
   neighbour and is shared as the same curve (contract S1; bitwise between fitted faces; no T-junctions); a flat region (one connected constant-z area
   at one height, merged across seams and mirror planes) is one plane face bounded only by rows of
   lateral faces. The mirror of a fitted face is its source's face with the control points
   flipped, exactly 0 in the flipped coordinate on a boundary in the mirror plane, and its source
   is cut at the vertices of both. A fitted outer face stays within ε/4 of the exact
   surface: to first order the thickness between the written faces is t + ε/2 minus the outer and
   the inner fitting error, so M1 and M2 hold for every sign when the two errors sum to at most
   ε/2, and the outer fit, made once per hull before any t, leaves each inner fit the same half.
   Inner faces whose offset keeps no structure of the outer patch are fitted the same way
   (interface contract §8, general hulls).
10. **Stage-2 constraints and starts.** Constraints: flotation equality, GM ≥ `gm_min`, adjacent
    density ratio. There is no monotonic-density constraint (`c_mono`) and no minimum-mass
    constraint (`c_mass_min`). Stage 2 also runs from the bottom-filled start of §3 item 26,
    keeps the start with the lower objective and logs both starts. Stability is enforced by
    GM ≥ `gm_min` in Stage 2 and GM = GM_Stage2 in Stage 3.
11. **Stage-3 objective, acceptance and closest fail, both modes.** Minimise
    Σ ((X3 − X2)/X2)² over X ∈ {Z_CG = `CG_total(3)` (world frame), coupled T_heave, coupled
    T_pitch}; equalities: flotation (to the solver's constraint tolerance) and GM = GM_Stage2.
    After every solve, the `mass_acceptable_pct` check on Z_CG, GM, T_heave and T_pitch decides
    accepted or failed. Escalation: UHPC — ballast within k* → ballast spill → draft; thin shell
    — (t, `z_ballast`) at the Stage-2 draft → draft. The draft is released only when mass balance
    cannot be met at the Stage-2 draft; a failed check does not release it. Closest fail: if the
    equalities hold and the check fails, the optimum itself; if the equalities cannot be met, the
    iterate with the smallest equality violation. Status flag and per-metric report as item 5.

---

### What the two equalities imply (derived, read before T5–T7)

At the Stage-2 draft, flotation fixes the mass (M = ρ_w V_sub) and GM = GM_Stage2 fixes Z_CG
(the draft fixes KM, and GM = KM − Z_CG). Heave then cannot change either: T_heave depends only
on M, A33 and A_w, all set by the draft. Only Iyy, and through it the pitch period, remains free.

The Stage-3 objective (§5 item 11) therefore acts through T_pitch only at the Stage-2 draft.

- **UHPC (C1):** variables at fixed draft are the ballast level in k*, t_k* and t_3, t_4 — four
  unknowns, two equalities, two degrees of freedom left for the pitch objective.
- **Thin shell:** variables at fixed draft are t and `z_ballast` — two unknowns, two equalities,
  **no freedom left**. Stage 3 becomes a solve, and the objective acts only if the draft is
  released (the last resort).

---

## 6. Testing and environment

- Cloud container (Ubuntu 24.04), branch `claude/lucid-cray-7o9442`. Commit and push after
  every task.
- `tools/install_toolchain.sh` installs GNU Octave 8.4.0 (`apt-get install -y octave
  gnuplot-nox octave-optim`) and Python `gmsh` 4.15.2 (`pip install gmsh`). It is idempotent;
  rerun it in every fresh container.
- Run Octave with `tests/octave_shims/` on the path **before** `src/`. Verified on 2026-10-06:
  with the `startsWith`/`endsWith` shims, Octave parses C1 (`MS2Parser.parse`, a static
  method), slices it, and reproduces the stored `.mat` contours to 1.3e-15 m.
- Missing or different in Octave 8.4 (need test-only shims, task T0): `fmincon` and
  `optimoptions` (Octave's optim `fmincon` rejects infeasible starts, so wrap core `sqp`),
  `contains`, `discretize` (Stage-1 `properties_2d`), `datetime` (`optim/run.m`), `save -v7.3`.
  Figure-only gaps: `tiledlayout`, `nexttile`, `exportgraphics`, `polyshape`/`polybuffer`,
  `xline`/`yline`, `sgtitle`, `delaunayTriangulation`, `histogram`. Keep these out of tested
  data functions; the owner confirms figures in MATLAB.
- Python `gmsh` (bundles OpenCASCADE) is a **test-only** tool for checking STEP files.
- Octave fidelity (T0 baseline, 2026-10-07): thin shell reproduces the MATLAB v1.0 results;
  modular precast does not (Stage 2 stops with exitflag −2 at vs = 0.409 against MATLAB 0.983).
  Owner: "You can't expect to match my matlab results when we're running different geometry
  and initilization." Octave results are not compared with MATLAB v1.0 numbers; Stage-2 and
  Stage-3 tests check the formulation directly (constraints and equalities at the returned
  point, solver exit flags), and the owner's MATLAB run is the final check. A full C1 run takes
  24–94 min in Octave, so full-pipeline tests run only with `MWEC_REGRESSION=1`, and no
  whole-pipeline run happens before Waves B–E are implemented.

---

## 7. Directory list

```text
MWEC_Mass/
├── AGENTS.md                  this file
├── README.md, CHANGELOG.md, CITATION.cff, LICENSE
├── WEC_User_Input.m           author inputs and entry point
├── WEC_Output_Options.m       output selection and figure style
├── Input/
│   ├── C1.ms2                 C1 hull (MultiSurf deck)
│   ├── C1_wamit_cache.mat     hydrodynamic coefficients
│   └── WAMIT/                 WAMIT preprocessing helpers
├── Output/                    reference results of the v1.0 C1 run (per realisation type)
├── docs/
│   ├── METHODS_ENGINE.md, RUNTIME_GUIDE.md, RESULT_SCHEMA.md, HAMS_MREL_ROUTE.md
│   └── plans/                 implementation plans (this work)
├── _graph/                    generated code map (Graphify)
├── validation/diagnostics/    optional post-processing scripts
└── src/+mwecmass/
    ├── +driver/               run, build_config, hydrostatic and strip tables, deck parsing
    ├── +geometry/             MS2Parser, iso-z slicing, boundary cache, profile helpers
    ├── +hydrostatics/         section and strip integrals, properties_2d/3d, mass matrix
    ├── +bem/                  hydrodynamic cache, interpolation, CG transfer, HAMS-MREL route
    ├── +mesh/                 BEM panel mesh and normals
    ├── +optim/                Stage 1 modes, Stage 2 objective/constraints/bounds
    ├── +realise/
    │   ├── +preliminary/      pass-through of Stage 2
    │   ├── +thin_shell/       steel shell + ballast solve
    │   ├── +modular_precast/  UHPC modules solve and strip extraction
    │   └── build_realised_properties.m, empty_realised_properties.m
    ├── +output/               export, schema, reports, logs, +figures/
    └── +internal/             small numerical helpers
```

New folders created by the plan: `src/+mwecmass/+solid/` (exact geometry kernel),
`src/+mwecmass/+output/+step/` (STEP writer), `tests/`, `tools/`.

---

## 8. Open decisions (ask the owner; do not decide alone)

- **OD1** Resolved: Stage 3 minimises the squared relative deviation from Stage 2 of Z_CG, coupled
  T_heave and coupled T_pitch (§3 item 27); the periods are coupled. Stage-1 and Stage-2
  objectives do not change.
- **OD2** Resolved: |Z_CG,3 − Z_CG,2| / |Z_CG,2| with Z_CG = `CG_total(3)` (world frame).
- **OD3** Resolved: the steel shell is the exterior surface as written (exact NURBS where the
  patch takes the exact path, fitted faces otherwise; interface contract §8) **from `z_ballast` up to the deck only**. Below `z_ballast` the mass model
  counts the full section as ballast (plate included), so a shell surface there would count the
  plate twice. The ballast solid is the full outer section below `z_ballast`. The two bodies
  meet along one closed curve: the section of the exterior surface at `z_ballast`, which is also
  the outer edge of the ballast's top face. The combined file holds both bodies with that curve
  identical in each, so a mesher can merge them; a single manifold solid cannot contain a 2D
  sheet.
- **OD4** Resolved: thin shell uses the closest-fail rule.
- **OD5** Resolved: adaptive fitting with metrics M1–M3 (§5 item 9).
- **OD6** Resolved (agent's decision, owner informed): option (ii), the ballast level measured
  from the next module's bottom with module edges unchanged. For the same ballast level, (i) and
  (ii) give the identical body (solid from the keel to `z_ballast`, shells above); they differ
  only in where the precast joint sits. (ii) keeps every module equal to its Stage-2 strip, so
  each keeps its Stage-2 density target and its own STEP file; it keeps one continuous variable
  whose upper bound simply widens from the top of k* to the top of k*+1 when spill is allowed
  (smooth for SQP); and it never changes module heights, which are fabrication inputs (including
  the 1.8 m wall module). (i) would move a module edge during the solve: a change of topology
  (non-smooth for SQP), reassigned Stage-2 targets, and silently altered module heights.
- **OD7** Resolved: t_min = t_init = 25.4 mm for thin shell; it also sets the Stage-2 floors.
- **OD8** Resolved: the draft stays a variable in both modes, used last (§3 items 20 and 31). It
  is released only when mass balance cannot be met at the Stage-2 draft (UHPC: after the spill
  step); a failed `mass_acceptable_pct` check does not release it.
- **OD9** Resolved: the Stage-2 upper density bound is the solid density of the mode's ballast —
  7500 kg/m³ (solid steel) for thin shell, 2500 kg/m³ (solid UHPC) for modular precast — taken
  from the material inputs, not from the shared `ballast_density_bounds`. ρ_air is 1.2 kg/m³ in
  both modes (`config.rho_air` from `thin_shell.rho_air`; `config.constructability_rho_air` from
  `modular_precast.rho_air`); T4a asserts it.
- **OD10** Resolved by the owner's wording (§3 item 32): mass = displaced mass is an equality
  held to the solver's constraint tolerance (1e-6) in every reported state; `mass_acceptable_pct`
  applies only to the Stage-2 comparison of Z_CG, GM, T_heave and T_pitch (see I21).
- **OD11** Resolved: GM = GM_Stage2 (equality) in both modes.
- **OD12** Resolved: delete `c_mono`; Stage 2 also runs from a bottom-filled start and keeps the
  start with the lower objective (§3 item 26, §5 item 10; agent's proposal for the safeguard, owner informed).
- **OD13** Resolved (agent's decision, owner informed): in the ballast module k*, the shell above
  the ballast has thickness t_k* ≥ t_min, a Stage-3 variable that starts at t_min (§3 item 28).
- **OD14** Resolved: thin-shell `t_max` = half the neck thickness (§3 item 30).
- **OD15** Resolved (agent's decisions, owner informed): delete `c_mass_min` and the UHPC
  Stage-3 pre-check (§3 items 33 and 34).
