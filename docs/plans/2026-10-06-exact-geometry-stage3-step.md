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
- **Interfaces:** `2026-10-07-interfaces.md` is the contract every task builds and tests against
  (structures S1–S8, functions F1–F14, invariants I1–I9, stand-in kit SK, file ownership map).
- **Lanes.** Wave A (T0, T1, T9) is done; from there tasks run in lanes and meet at joins
  (contract §6). Tasks in different lanes run in parallel; inside a lane they run in order.
  - J0: merge T0, T1, T9 and the contract.
  - Lane K: SK stand-in kit, then T2a (exact path), then T3; T2b (general path) runs in parallel
    with T3, both starting from T2a's branch with disjoint files; join J1 merges T2a, T2b and T3.
    SK2 (stand-ins follow the general-hulls amendment) runs right after that amendment merges,
    before J1. **Owner checkpoint after T3.**
  - Lane N: T0b, merged as soon as accepted (every later lane starts from the new names); T0c in
    parallel (it edits only `MS2Parser.m`).
  - Lane P: after T0b, SK and T0c: T0d, then T4a, then T4b (all edit `build_config.m`); T4b
    starts after SK2 (it calls F1 with three inputs) and merges after J1. Lanes U, S and O need not
    wait for SK2 (contract §6: its changes are additive for them).
  - Lane U: after T0b and SK: T5, then T6. Lane S: after T0b and SK: T7. Lane O: after T0b and
    SK: T8 and T10.
  - J2: after J1 and lane P, the group T5, T6, T7, T8, T10 merges as one chain in that order; the
    first whole-pipeline runs follow. **Owner checkpoints after T6 and after T10**, on those runs.
  - Lane G, after J2: T11 cleanup, then T12 docs, then T13 final review.
  - Before a producer is merged, its consumers test against the stand-ins of SK (contract §3).
    A shared file is edited in the order of the contract's file ownership map (§4), each edit
    after rebasing on the earlier ones.
- **Baseline:** the Octave regression baseline (T0, `tests/baseline/`) may change only through
  tasks that intend to change results (T4a, T4b, T5, T6, T7). Renames (T0b), parser lookups (T0c)
  and the geometry cache (T0d) leave every number identical. Full-pipeline tests run only with
  `MWEC_REGRESSION=1` (a C1 run takes 24–94 min in Octave). Owner (2026-10-07): no
  whole-pipeline run before T5–T7 are implemented; components are tested on their own. Tasks are
  therefore graded on component tests, and acceptance items that need a whole-pipeline run are
  checked after J2: the baseline on the integration commit after T0d must reproduce exactly
  (T0b, T0c, T0d), then C1 runs of both modes on the J2 head serve T4a, T4b, T5, T6, T7 and T10,
  and the baseline is regenerated once with the changed quantities printed per task.
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

Contract: implements the test harness every item uses (contract §0, Tests); consumes nothing.

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

Contract: implements the names of contract §0 and I5 in existing code; first in the contract §4 order of every
shared file.

Lane N, merged before lanes P, U, S and O start, so every Stage-2 and Stage-3 task uses the new
names. Pure renaming: no change to any
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

Contract: no contract item; `MS2Parser` API unchanged (consumed by F1, F2).

Lane N, parallel to T0b; owner request (2026-10-07). Files: `src/+mwecmass/+geometry/MS2Parser.m`,
new tests and fixtures.

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

Contract: no contract item; its cached block later holds `config.hull_solid` (S1b) and F8 (added by T4b).

Lane P, after T0b and T0c; owner request (2026-10-07). Files: `src/+mwecmass/+driver/build_config.m`, one new helper
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

Contract: done; `outer_rows` and `surface_normals` are consumed by F2 and by the I9 test.

Files: `src/+mwecmass/+solid/outer_rows.m`, `surface_normals.m`, tests.

- Sample each `.ms2` patch so that rows are horizontal sections at requested heights (module
  edges, fill level, waterline, quadrature nodes) and columns follow the patch parameter.
- Normals from S_u × S_v using `MS2Parser.eval_surface` (analytic derivative where the parser
  provides it, central differences otherwise), oriented outward; handle mirror patches, seams
  and the undefined normal at the keel point and the neck top.

Acceptance: sections equal `extract_isocurve_at_z` and the Python reference; normals are unit,
outward, continuous across patch seams; exact symmetry under the deck's mirror planes.

### SK2 — Stand-ins follow the general-hulls amendment (Sonnet high)

Contract: §3 SK2 and §6. Edits `tests/standins/*` and the SK tests that mirror the changed items
(`tests/standin_kit/test_sk_outer_nurbs.m`, `test_sk_offset_slice.m`); nothing in `src/`. Runs right
after the general-hulls amendment (owner, 2026-10-07: "No it needs to be generalized.") merges,
before J1.

- F1 stand-in: `outer_nurbs(model)` and `outer_nurbs(model, cache, opts)`, `cache` and `opts`
  ignored. F6b stand-in: `body_section(body, z)` and `body_section(body, z, side)`; `'above'` (the
  default) is its present half-open rule, `'below'` takes the faces and the module below a module
  edge or `z_ballast` (contract F6b). S1 entries carry `visible`, `fit` (`[]`) and `swap_uv`
  (false): in `geo.outer`, `visible` = k (one entry per visible surface on both fixtures); in the
  F2 stand-in's S2 set, each piece's
  `visible` is that of the outer entry it is offset from (the cylinder: the three crease pieces of a
  quarter take that quarter's index, so the 12 pieces take 1 to 4). S1b and the F2 stand-in's S2 set
  carry `flat` (empty).
- F2 stand-in: the box returns `sti_inner_box(geo, t, opts.t_min)` (the set the real F2 must equal,
  T2a) instead of `FitNotConverged`; that branch and the assertion on it are retired.

Acceptance: every SK test passes with the field-set assertions extended by the new fields and the
box assertion replaced by a comparison with `sti_inner_box`; F6b at a module edge of a stand-in
body gives the module above by default and the module below with `'below'`, the same outer loop;
`outer_nurbs(model, [], struct())` equals `outer_nurbs(model)` for both fixtures; the SK acceptance
items of contract §3 still hold.

### T2a — Kernel B, exact path: normal offset, fold trimming, spline surfaces (Sonnet high)

Contract: implements S1, S1b, S2, S2r, F1, F2, F2b, F3, F3b, F4 and invariants I3, I4, I9 on the exact
path; consumes T1 (`surface_normals`, `outer_rows`) and the SK fixtures. Starts after SK in lane K; the
critical path to T3. Scope per contract F1, F2 and §8 (general hulls; owner, 2026-10-07: "No it needs
to be generalized."): offset at d = t + eps_fit/2; the exact path of F1 (deck check with
`UnsupportedEntity` only for an unparsed entity the hull references, `UnsupportedMirror`,
`HullNotClosed`; u and v exchanged where z depends on v alone, `swap_uv`; split at constant-z
intervals inside the z range; `FitInputMissing`; the options `max_passes`, `force_general`);
structured inner pieces (`rev_z`, `ruled_parallel`); creases split along u and v (convex: trimmed;
concave: a face of the crease curve offset along its fan of normals); at a vertex whose cone of
normals spans a solid angle, the sphere face that closes the gap between the fan faces (contract
F2; found by T2a's F2, fitted by T2b's `fit_z_faces`); seams shared bitwise;
`c0_u`/`c0_v` knot rows in S1 and S2; F2 `opts.knots_from` (refit at a new t on fixed knots); the
cylinder fixture is its join test, the box is not (T2a still runs the real F2 on it, below).

Interface to T2b: F1 and F2 hand every patch or piece of the general path to
`[faces, flat, rep] = mwecmass.solid.fit_z_faces(model, cache, faces, general, opts)` (T2b's file):
`faces` is the S1 array of the exact path (F1) or of the inner pieces (F2; it also holds one entry
per fan face that keeps no structure and per vertex sphere face, with `surf` empty and a field
`offset_of` naming what it offsets: `crease`, the crease curve and the outer entries on its two
sides, or `vertex`, the point and the outer entries around it, from which `fit_z_faces` computes the
exact offset points), `general` the logical mask of the entries that take the general path, `opts` `t_min`, `max_passes`, `kind` (`outer` |
`inner`) and, for `inner`, `d` and the outer `geo`. It returns the final S1 array (exact entries split
at the band ends, general entries replaced by fitted faces, further entries moved to the general path
by the Cuts and flat-region, tie-break and seam rules, in F1's order), the S1b/S2 `flat` regions
and the fit report of its faces.
F1 and F2 call it only when an entry takes the general path; T2a tests only decks with none (on T2a's
branch alone such a deck stops with the undefined-function error until J1).

Files: `src/+mwecmass/+solid/outer_nurbs.m`, `offset_surface.m`, `trim_fold.m`,
`fit_bspline_surface.m`, `eval_bspline_surface.m`, `eval_bspline_curve.m`,
`split_bspline_surface.m`, `slice_bspline_surface.m`, `void_closing_distance.m`, tests,
`tests/solid/fixtures/stepped_spar.ms2` and `box_swapped.ms2`.

- Adaptive, error-bounded fitting as specified in AGENTS.md §5 item 9: initial nodes dense where
  curvature is high or the void is narrow; offset by t + ε/2 along the exact normal; trim the
  fold; split faces along creases; fit cubic B-splines; check M1–M3 on points between the nodes;
  insert knots only in failing spans; remove unneeded knots at the end; stop with an error report
  at the iteration cap.
- Slice the fitted surface at any height into an ordered closed contour.

Acceptance: M1–M3 pass on a dense check grid not used for fitting (report min / max t_local, the
number of refinement passes and knots per face); fold trimming and crease splitting tested where
folds occur (C1 at t ≥ 0.100 m, and a stand-in fixture with a convex radius below t); C1 at
t = 0.0762 m compared with the Python reference half-widths per height and with the independent
erosion result for module 4 (void 3.548 m³ with the v1.0 module edges and `z_ballast`); a
slender-section case (thin-shell neck, t = 0.025 m) passes M3. C1 keeps its exact-path oracles (all 8
patches exact, contract §3; `outer_nurbs(model)` and `outer_nurbs(model, cache, opts)` give the same
`geo`). The box: the real F2 set equals `sti_inner_box` (convex C0 v-seams, trimmed along parameter
lines, structure kept), compared and printed.

Exact-path additions (contract F1). `stepped_spar.ms2` (format of SK's `cylinder.ms2`): FramePoints
K (0, 0, −3), P1 (1.5, 0, −3), P2 (1.5, 0, −1), P3 (0.75, 0, −1), P4 (0.75, 0, 1), T (0, 0, 1); a
PolyCurve2 of the Lines K–P1, P1–P2, P2–P3, P3–P4, P4–T revolved 0° to 90° about the Line K–T;
`Symmetry: x y`. Verified 2026-10-07: `MS2Parser` parses it (4 patches) and T1 `outer_rows` (grid 60)
returns one closed loop at eight heights off the step in [−2.95, 0.95] m (shoelace areas 2.25π and 0.5625π m² less the inscribed-polygon error, 1.2e-4 relative).
Asserted: each quarter is `rev_z` and exact, split by F1 at the two PolyCurve2 joints
that bound the step (existing C0 knots: no control point changes) into three consecutive S1 entries
with `z_range` [−3 −1], [−1 −1] and [−1 1], the cut rows shared bitwise; F3b cuts the outer pieces
at z = −2 and z = 0, and on the unsplit converted patch at z = −1 raises `ZNotMonotonic`; F4 at
z = −1 gives the circle of radius 1.5 from the faces below and of radius 0.75 from the faces above;
F2 at t = t_min trims the convex rim P2, gives the concave corner P3 a face of its own, splits the
inner shelf as F1 does, and passes M1–M3; inner volume printed next to the closed form of the
offset profile. `box_swapped.ms2`: SK's `box.ms2` with each side face a RuledSurf between its bottom
and top edge Lines (new Lines B1–B2, T1–T2, …) instead of its vertical Lines. Verified 2026-10-07:
the parser's z of `box_side1` depends on v alone, and T1 `outer_rows` gives the area 3 m² of
`box.ms2` at z = −2, −1, 0. Asserted: z of the
side faces depends on v alone, F1 exchanges u and v (`swap_uv` true; control points and weights
transposed bitwise), they take the exact path, the parser's edge e is S1 boundary 5 − e, the seams
pair up bitwise (I2), and F4 sections at z = −1 equal those of `box.ms2` bitwise. Deck check, on
copies written by the test: `box.ms2` plus a `Variable` line and a line of an unknown type that
nothing references gives the same `geo`; `FramePoint T1` given an unknown type (referenced by `V1`)
raises `UnsupportedEntity`; `box_top` given an unknown type raises `HullNotClosed`; `cylinder.ms2`
with `Symmetry: z` raises `UnsupportedMirror`.

### T2b — Kernel B, general path: faces fitted with z as a parameter (Sonnet high)

Contract: implements §8 (general hulls) through `fit_z_faces.m` (interface in T2a): band ends for
the whole hull, faces fitted through exact points (T1 sections and normals) with z as one parameter,
split where z turns back and at creases, closed-loop cuts and vertex cuts (no T-junctions), merged
flat regions, mirrors (flipped faces, vertex union), seams with exact neighbours, outer faces within
ε/4 of the exact surface, M1–M3 between the written faces; the same path for inner pieces whose
structure is not kept; F1 `opts.force_general`. Starts from T2a's branch, in parallel with T3; its
files are disjoint from T2a's and T3's; J1 merges T2a, T2b and T3.

Files: `src/+mwecmass/+solid/fit_z_faces.m` (and helpers named `fit_z_*.m`), tests,
`tests/solid/fixtures/tilted_revolution.ms2`, `tilted_revolution_two_loops.ms2`,
`lying_revolution.ms2`, `split_side_cylinder.ms2` and `cross_column.ms2`.

- General path (contract §8): outer patches that do not take the exact path, and inner pieces whose
  structure is not kept, become untrimmed faces fitted through T1 sections and normals, with z as
  one parameter; outer fitted faces within ε/4 of the exact surface (contract §8 derivation); M1–M3
  judged between the faces as written. Where a face's part of a section is a closed loop with no
  seam or crease point, the face is cut along one z-monotone curve on the exact surface
  (steepest-ascent line of z), used bitwise as both v-boundaries (self-seam), one curve through
  consecutive closed-loop bands; a band that ends at a single highest or lowest point inside a patch
  ends in a pole row there. Band ends are one set of heights for the whole hull (exact patches split
  there by F3b, or at their existing knots where the band end is the height of a constant-z interval),
  and every vertex that would lie inside a face row (the end or z-extreme of a seam or crease, the end
  of a closed-loop cut, a vertex of a mirror face's rows flipped back) is joined by a cut to a vertex
  of the face's opposite row, so every face boundary has one neighbour and is shared bitwise (contract
  §8 Cuts, Mirrors). A flat region is one connected constant-z area at one height, merged across seams
  and mirror planes, bounded only by lateral-face rows, and written by F5 as one plane face. The mirror
  of a fitted face is its source's face with the control points flipped, exactly 0 in the flipped
  coordinate on a boundary in the mirror plane (contract §8). The fan faces of concave creases that
  keep no structure and the sphere faces of vertices whose cone of normals spans a solid angle (F2)
  are fitted through their exact offset points in the same way.

Acceptance:

Non-C1 hull (general path). `tests/solid/fixtures/tilted_revolution.ms2`: the profile of BCurves
through, in the axis frame (r, h) [m], T (0, 1), D (0.5, 1), N (0.5, 0.8), S (1.0, 0.4), Q (1.0, −1),
K (0, −1), one RevSurf of 360° per segment about the Line K–T, every point rotated by 15° about the
y axis (x' = x cos 15° + z sin 15°, z' = −x sin 15° + z cos 15°), written in the format of
`capped_cylinder.ms2`. z depends on both surface parameters, the rims at D, S and Q are convex C0
creases, N is a concave one, and the creases are tilted circles, so bands end at their z-extremes.
Verified 2026-10-07: T1 `outer_rows` (grid 60) returns one closed loop at each of 38 heights through
its z range [−1.2247, 1.0952] m. Oracles: the solid is a rigid rotation of a solid of revolution, so
the distance of any point from the exact outer surface, and t_local, are 2D distances to the profile
in the meridian plane of the axis frame (lines; the inner profile is the 2D offset by d, trimmed at the
convex corners, with an arc of radius d about N). Asserted: F1 marks every patch fitted; each fitted
outer face has `fit.dev_max` ≤ ε/4 and the oracle distance at the check points is ≤ ε/4; M1–M3 pass
for a shell at t = t_min (the offset folds at every convex corner) as the kernel reports them between
the written faces (S2r); at the same check points the oracle t_local (to the exact profile) lies in
[t − e, t + ε + e], e the largest oracle distance of the written outer faces from the exact surface
(≤ ε/4, asserted above), and its difference from the kernel's t_local is printed per face; every
face is `z_of_u` with monotone z and F3b cuts it at a module-edge height; seams pass the I2 bitwise
test and the seam fields pair up (contract I2). No T-junction (contract §8 Cuts): at each of the
highest points of the N, S and Q rims (φ = 180°, where each rim's crease branches end and no seam
runs) and at the disk centres T and K (where the self-seams of the top and bottom disks end), every
face whose closure contains the point has it as a corner, never inside a row, and the boundaries
that meet there pair one to one and bitwise. The closed outer body: a test helper assembles the
S1 faces into a T9 brep (one edge per paired boundary, its curve the shared row), which passes
`validate_brep`, is written with `write_step` and passes `tests/step_check.py` (one closed solid,
METRE); its volume is printed next to the closed-form volume of the solid of revolution,
π∫r² dh = 1.68333·π = 5.2884 m³. Printed: deviations, t_local range, passes and knots per face. The same deck with
N (0.5, 0.6) and S (1.0, 0.5) (`tilted_revolution_two_loops.ms2`) raises `SectionNotClosed`
(verified: two loops at z = 0.738 m).

Closed loops and poles inside a patch. `tests/solid/fixtures/lying_revolution.ms2`: FramePoints
A (0, −1, 0), B (0.6, −1, 0), C (0.6, 1, 0), E (0, 1, 0); `BCurve prof` of degree 3 on { A B C E };
`Line AX` from A to E; `RevSurf hull` of `prof` about AX, 0° to 360° (one patch, format of
`capped_cylinder.ms2`). A smooth body of revolution about the horizontal y axis: radius
1.8 s(1 − s) at profile parameter s, poles on the axis at y = ±1 (z = 0); its self-seam (φ = 0) and
both poles lie at z = 0; its top and bottom, z = ±0.45 (s = 1/2, v = 3/4 and 1/4), are single
points inside the patch; smallest principal radius 0.253 m, so no fold at t_min. Verified
2026-10-07: T1 `outer_rows` (grid 60) returns one closed loop on the one patch with no seam point at
each of 37 heights in [−0.449, 0.449] m. Asserted: two bands, [−0.45, 0] and [0, 0.45]; each face's
part of a section is the whole loop, so each face has one self-seam (v0 and v1 the same curve
bitwise, z-monotone, its edge once in each direction in the face's loop, `validate_brep` passes),
the two self-seams are one chain cut (contract §8 Cuts: they meet at z = 0 in one point, the
bottom face's top row and the top face's bottom row starting there, bitwise equal) and each face
has a pole row at z = ±0.45 (all control points bitwise equal); every face `z_of_u` with monotone z,
cut by F3b at z = ±0.2; oracle as above (2D distances to the profile in the meridian plane, closest
point by Newton on the exact cubic; inner profile the 2D offset by d), `fit.dev_max` ≤ ε/4, M1–M3 at
t = t_min with the oracle check of the previous paragraph. Printed as above.

Horizontal-tangent rows, against the exact path. C1 with `opts.force_general` true (contract F1):
every patch fitted; band ends at C1's horizontal rows, the shoulder z = −1 (verified 2026-10-07:
dz/du of `surface1` vanishes there while dr/du does not) and z_max = 1.1 (ridge row of `surface2`,
top pole of `surface1`), and at the keel z_min = −3.25; the oracle is the exact path's `geo` (exact
NURBS of the same entities, I9): distance of every fitted face from it ≤ ε/4 at the check points,
F4 sections, hull volume and CG of both printed side by side; mirrors bitwise in the symmetry
planes (I2); M1–M3 at t = 0.0762 m between the written faces.

Flat regions (contract §8). `split_side_cylinder.ms2`: the points and Lines of SK's `cylinder.ms2`;
RevSurfs `side_a` (side Line, 0° to 90°), `side_b` (90° to 180°), `bottom` and `top` (the disk Lines,
0° to 180°) about the Line K–T; `Symmetry: y`. Every patch converts exactly, but the end of the seam
`side_a`–`side_b` (90°) lies inside the rim rows of `bottom` and `top`, so these and their mirrors go
general and become flat regions (contract F1 order of the rules, step 1a); the side faces have no
vertex inside a row and share their vertical seams end to end, so they stay exact whatever the
patch names (the tie-break of step 2 never applies). Verified 2026-10-07: `MS2Parser` parses it (8
patches with the
mirrors) and T1 `outer_rows` (grid 60) returns one closed loop (shoelace area 2.25π m² less 1.2e-4 relative) at nine heights in
[−2.95, 0.95] m. Asserted: `geo.outer` holds the four exact side faces; `geo.flat`
holds two regions, z = −3 with `normal_z` −1 and z = 1 with +1, each covering a disk and its mirror
(`visible`), merged across the mirror plane y = 0; every top and bottom row of a side face names its
region (`[0 j]`) and no other boundary does; the test helper's brep (one plane face per region, its
loop the chain of those rows) passes `validate_brep`, is written with `write_step` and passes
`tests/step_check.py` (one closed solid, METRE); its volume is printed next to 9π m³.

Concave vertices (contract F2). `cross_column.ms2` (`Symmetry: x y`; RuledSurfs between Lines, the
quarter x, y ≥ 0): a square pontoon of half-width 1.5 m from z = −3 to −1 (bottom, walls x = 1.5 and
y = 1.5, its top as the rectangles [1, 1.5] × [0, 0.3], [0.3, 1.5] × [0.3, 1] and [0, 1.5] × [1, 1.5])
carrying a cross-shaped column of arm half-width 0.3 m and arm half-length 1.0 m from z = −1 to 1
(walls x = 1 and y = 0.3 of one arm, x = 0.3 and y = 1 of the other, each a RuledSurf between two
vertical Lines; its deck as the rectangles [0, 1] × [0, 0.3] and [0, 0.3] × [0.3, 1]). Verified
2026-10-07: `MS2Parser` parses it (48 patches with the mirrors) and T1 `outer_rows` (grid 60) returns
one closed loop at ten heights off z = −1 in [−2.95, 0.95] m, area 9 m² below z = −1 and 2.04 m²
above. At the foot of each of the four re-entrant column edges, e.g. (0.3, 0.3, −1), three concave
creases meet (the re-entrant edge and the two foot lines), so the cone of normals there is the octant
n_x, n_y, n_z ≤ 0. Asserted: the walls and the bottom stay exact, whatever the patch names (the
pontoon-top rectangles have corners inside the wall x = 1.5's top row and inside each other's rows,
e.g. (1, 0.3, −1) inside a row of [0.3, 1.5] × [0.3, 1], so they merge into one flat region (contract
F1 order of the rules, step 1a) before the vertices are counted, and their corners then put no vertex
on the walls' rows (step 1b); likewise the deck); the pontoon top and the deck are
flat regions (`geo.flat`: z = −1, `normal_z` +1, whose loops are the pontoon walls' top rows and,
as a hole loop, the column walls' bottom rows; z = 1, +1); F2 at t = t_min gives a fan face at every
concave crease and one sphere face at each of the four vertices, its points at distance d from the
vertex, its boundaries the end arcs of the three adjacent fan faces (bitwise, I2) and a pole row at
its lowest point (0.3, 0.3, −1 − d); no such face at the convex column corners, e.g. the foot
(1.0, 0.3, −1), where the fans overlap and are trimmed; M1–M3 pass between the written faces; the
oracle t_local (distance from a check point to the outer polyhedron, exact: to its planar faces) lies
in [t, t + ε], since the outer faces are exact (e = 0), and its difference from the kernel's
t_local is printed per face. The test helper's breps of the outer faces and of the inner faces each
pass `validate_brep`, are written with `write_step` and pass `tests/step_check.py` (one closed
solid, METRE); the outer volume is printed next to 22.08 m³, the inner volume printed.

### T3 — Kernel C: bodies and exact properties (Sonnet high)

Contract: implements S3–S7, F5, F6, F6b, F7 (lateral faces split at planes and `c0_u`/`c0_v` rows,
S4) and invariants I1, I2, I6–I8 (including collapsed-row loops and
a `write_step` + `step_check` import of C1 bodies); consumes F1–F4 (T2a; the general path of T2b joins at J1), T9 `write_step`/`validate_brep`.
Methods per contract §3: divergence-theorem surface integrals (F6) and hydrostatics on the cut outer
patches (F7, including `full`/`none` submersion and `S_wet`) replace the section integration below.

Files: `src/+mwecmass/+solid/build_body.m`, `body_properties.m`, `hydrostatics_at_draft.m`,
tests, `tests/solid/fixtures/stepped_box.ms2`.

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

Planes at the height of a flat part (contract F5; exact path, so T3 tests them before J1).
`stepped_spar.ms2` (T2a), a module edge at the shelf height z = −1: F4 gives the circle of radius 1.5
from the faces below and of radius 0.75 from the faces above; the cap there is the disk of radius
0.75 (the intersection of the two loops' areas), and the shelf piece (`z_range` [−1 −1], the ring
between the circles) is written once, as an outer face of the module below; with both modules solid
and with both hollow at t = t_min (the cap then the annulus between the radius 0.75 and the inner
loop), each module and the fused hull pass `validate_brep` and `step_check.py`, I1 holds, and the
solid module volumes are printed next to 4.5π and 1.125π m³ (rational faces). F6b at z = −1 (the
module edge on the shelf) raises no error: the default side `'above'` gives the module above and the
outer loop of radius 0.75, `'below'` the module below and the loop of radius 1.5 (contract F6b),
each equal to the F4 loop of that side.
`tests/solid/fixtures/stepped_box.ms2` (`Symmetry: y`; RuledSurfs between Lines; the half y ≥ 0): a
box x ∈ [−1, 1], y ∈ [−0.75, 0.75] from z = −2.5 whose deck lies at z = 0.5 for x < 0 and at z = 1.5
for x > 0, with a riser at x = 0; eleven patches: the bottom in two pieces split at x = 0 (so that
the common corner (0, 0.75, −2.5) of the two lower pieces of the wall y = 0.75 is a corner of every
patch through it, not a point inside a bottom row), the wall x = −1, the wall y = 0.75 in three
pieces (below z = 0.5 for x < 0 and for x > 0, above it for x > 0), the wall x = 1 below and above
z = 0.5, the riser, the low and the high deck; every patch exact, every seam one curve shared end to
end, no corner of any patch inside another patch's boundary (contract §8 Cuts). Verified
2026-10-07: `MS2Parser` parses it (22 patches with the mirrors; no patch corner lies strictly inside
an edge of another patch, checked on the parser's corner points) and T1 `outer_rows` (grid 60)
returns one closed loop at eight heights off z = 0.5 in [−2.45, 1.45] m, area 3 m² below z = 0.5
and 1.5 m² above. With a module edge at the step z = 0.5, where the loops from below and from
above share the
rows along x = 1 and along y = ±0.75 for x > 0: no error; the cap is the area of the loop from above,
bounded by the shared rows (each one edge, taken once) and the riser's bottom rows; the low deck is
written once, as an outer face of the module below; with both modules solid, `validate_brep`,
`step_check.py` and I1 pass, and the module volumes equal 9 and 1.5 m³ (polynomial faces, asserted
exact). (Its inner sets have a fan face at the concave riser foot, which may take the general path,
so the hollow and thin-shell cases run at J1.)

**Owner checkpoint after T3.**

### J1 — Join of lane K (orchestrator)

Merges T2a, T2b and T3 in that order, deletes the stand-ins of F1–F7 and F6b, and reruns every
consumer test on the real code (the join test, contract §3). Adds `tests/solid/test_join_general.m`
(J1's file), which runs the general-path fixtures of T2b through the real F1, F2, F5, F6 and F7:
`split_side_cylinder.ms2` with module edges at z = −2 and z = 0: the two flat regions of `geo.flat`
are F5 plane faces whose loops are the chains of the rows that name them; each module and the fused
hull pass `validate_brep` and `step_check.py` (closed solids, solid count); I1 holds; the hull volume
is printed next to 9π m³. `cross_column.ms2` with a module edge at the pontoon top z = −1 (the loops
from below and from above, square and cross, strictly nested): the cap there is the cross's area and
the pontoon-top flat region is written once; hollow modules at t = t_min carry the inner sphere faces
at the four concave vertices; `validate_brep`, `step_check.py` and I1 pass; the hull volume is
printed next to 22.08 m³. `stepped_box.ms2` (T3) with the module edge at the step z = 0.5 and both
modules hollow at t = t_min: the cap is the area of the loop from above less the void, the shared
rows are one edge each, and `validate_brep`, `step_check.py` and I1 pass; thin shell with
`z_ballast` = 0.5: the ballast body is closed by the low deck and a `ballast_top` over the area of
the loop from above, the shell sheet starts at that loop, `V_ballast` = 9 m³ (asserted: polynomial
faces), and I1 holds.

### T4a — Stage-2 changes that do not need the kernel (Sonnet high)

Contract: no kernel item; edits shared files in the contract §4 order (after T0d). Stage-2 runs are checked after J2;
the per-task baseline update below is superseded by contract §6 (baseline regenerated once, after J2).

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

Contract: implements F8 and `config.hull_solid` (S1b, inside the T0d cache); consumes F1, F2, F5, F6
(SK stand-ins until J1; join test on the cylinder fixture). Merges after J1. The per-task baseline update below is superseded by contract
§6 (baseline regenerated once, after J2).

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

Contract: implements S8 and its schema for both modes, F9, F10 and F14 for modular precast (including the
F11 and F13 calls); consumes F2, F2b, F5–F7 (SK until J1; join test on the cylinder fixture),
`config.hull_solid`. Paths of F9, F10 per contract §3; `tools/baseline_run.m` reads `results.stage3` (§4).

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

Contract: extends F14 (precast escalation, closest fail) and fills S8 `escalation`, `solver`, `check`;
consumes S8, F9, F10 (T5), F2 (with `opts.knots_from`), F2b, F5–F7 (join test on the cylinder
fixture); caching per contract §7; t_max,i = d_close − eps_fit/2 (F2b,
contract §8).

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

Contract: implements F14 for thin shell (S8 with `ballast`, `shell`, `air`, including the F11 and F13
calls); consumes F2 (with `opts.knots_from`), F2b, F5–F7 (SK until J1; join test on the cylinder
fixture), F9 and F10 (SK until J2), `config.hull_solid`; caching per
contract §7; t_max = d_close − eps_fit/2 (F2b, contract §8).

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

Contract: implements F12, F13; consumes S8 (SK `sti_realised`), F6b.

Files: `src/+mwecmass/+output/+figures/` (new `realised_section_data.m` returning the y = 0 slice
and plan sections of the realised body; `plot_modular_precast.m`, `plot_steel_solve.m`,
`plot_optimised_cross_section.m`), `validation/diagnostics/stage_animations.m`; delete D8, D10,
D13.

- Figures draw true sections of the realised (or closest-fail) solid, show the ballast level,
  and label the status. Data functions run in Octave; drawing stays MATLAB code.

Acceptance: data functions tested in Octave against T3 sections; the owner confirms the figures
in MATLAB.

### T9 — STEP writer (Sonnet high) — wave A (separate files, done)

Contract: done; its B-rep struct is the `brep` of S4, consumed by F11 and the T3 tests.

Files: `src/+mwecmass/+output/+step/write_step.m` and helpers, tests.

- Outer faces: exact NURBS conversions of the `.ms2` entities where the type allows (C1: B-spline
  curves, arcs, revolution, ruled surface); inner faces: the adaptive fits from T2a and T2b.
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

Contract: implements F11 and `out.save.stage3.step` (also in `dispatch.m` `validate_save_flags`, after T8,
§4); consumes S4, S8 (SK `sti_realised`), T9 `write_step`.

Files: realisation `run.m` files, `WEC_Output_Options.m` (new `out.save.stage3.step_*`
switches), `src/+mwecmass/+output/+step/` builders.

- UHPC: one STEP per module; one STEP with all modules as one connected solid (built directly as
  one B-rep with the cavity as a void shell).
- Steel: ballast solid = full outer section below `z_ballast`; shell = the exterior parametric
  surface (exact NURBS where the patch takes the exact path, fitted faces otherwise, contract §8)
  from `z_ballast` to the deck only (no double counting of the plate);
  combined file with both bodies sharing the identical junction curve at `z_ballast`. For C1 the
  split at `z_ballast` is exact (z depends only on the profile parameter of the RevSurf and the
  RuledSurf, so the cut is an iso-parameter line found by knot insertion).
- Files go to `Output/<type>/step/`; paths stored in the `.mat`.

Acceptance: C1 files pass `tests/step_check.py`; imported volumes equal the kernel volumes
(reported). **Owner checkpoint after T10.**

### T11 — Remaining cleanup (Sonnet medium)

Contract: edits shared files last in the contract §4 order (before T12); checks I5 over the whole code.

`offset_polygon.m` if orphaned (D11); remaining D9 fields in `export_schema.m` and `Report.m`;
check that no shape-assumption pattern remains (`grep` for `sqrt(.*/pi)`, `pi\s*\*\s*r`,
`compute_rmin`, `R_eq`, `max_slope`, `cos_alpha`, `L_k`, `homothetic`).

### T12 — Documentation (Sonnet medium)

Contract: documents S1–S8 and F1–F14 as implemented; last in the contract §4 order of every doc file.

Rewrite the methods, runtime and schema sections listed under D12; document the kernel, the new
Stage 3, the status/report fields and the STEP outputs; update `CHANGELOG.md` and `AGENTS.md`.

### T13 — Final review (Opus high)

Contract: checks every invariant I1–I9 and the contract §4 order on the merged branch.

Whole-branch review against AGENTS.md §1 and §5; full Octave test run; a short list of what the
owner must confirm in MATLAB (real `fmincon` runs, figures).
