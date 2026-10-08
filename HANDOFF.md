# HANDOFF — MWEC_Mass exact-geometry refactor

Pause point of 2026-10-08 (UTC), written for the owner, who continues locally in MATLAB, and for any
agent that picks the work up. Every branch named here is on GitHub (`VKotzamanis/MWEC_Mass`).
Section 9 holds the branch heads at the moment of the pause.

## 1. Goal and binding documents

Make the MATLAB code correct by the owner's rules: every geometric operation from the `.ms2`
parametric definition (no shape assumptions), Stage 3 realised from the whole Stage-2 solution on
one exact geometry kernel (UHPC modular precast and steel thin shell as separate pipelines), STEP
output of the realised bodies. In one line: Stage 3 = the Stage-2 design built from real material,
checked against Stage 2.

Binding documents, in this order:
- `AGENTS.md`: owner rules (§1, rules 1–12), the v1.0 code (§2), owner corrections and decisions
  (§3 items 1–36), issues I1–I23 (§4), expected outcome (§5), decision log (§8).
- `docs/plans/2026-10-07-interfaces.md`: interface contract (structs S1–S8, functions F1–F14,
  stand-ins, shared-file order §4, invariants I1–I9, lanes §6, general hulls §8). It wins over the
  plan where they differ. Section 6 below lists decisions not yet written into it.
- `docs/plans/2026-10-06-exact-geometry-stage3-step.md`: task sections T0…T13 and the deletion
  register D1–D17.

## 2. What is where

### 2.1 Main branch `claude/lucid-cray-7o9442` (merged, tests green)

| Merged | Content |
|---|---|
| T0 | Octave test harness `tests/run_tests.m`, Octave shims, baselines, Python reference evaluator `tests/reference/c1_reference.py` |
| T1 | exact outer rows and normals of the `.ms2` surfaces |
| T9 | AP214 STEP writer `+output/+step/write_step.m`, `validate_brep`, `tests/step_check.py` (gmsh import check) |
| T0c | `MS2Parser` resolves each entity once (4.3× faster, bit-identical outputs) |
| T0b | renames: `z_ballast`, `rho_air`, `rho_ballast`, `uhpc` names (rule 10) |
| SK, SK2 | stand-in kit `tests/standins`, `tests/standin_kit` (closed-form cylinder and box bodies) |
| T0d | `build_config` saves its geometry products (`Output/cache/`, gitignored; key = deck, inputs, source) and reloads them (C1: 0.2–0.44 s instead of 6–7 min in Octave) |
| T4a | Stage 2: `c_mono` and `c_mass_min` deleted in both copies; second, bottom-filled start; upper density bound per mode from the material (7500 thin shell, 2500 UHPC); thin-shell t_min = t_init = 25.4 mm |
| — | `tests/run_tests.m` now runs under MATLAB too (shims only in Octave; an interactive session stays open) |

Last full suite on main (Octave, after T4a): 30 passed, 0 failed, 2 skipped (the two pipeline
regressions, which run only with `MWEC_REGRESSION=1`).

### 2.2 Task branches (not merged)

| Task | Branch | Status at the pause | What it holds |
|---|---|---|---|
| T5 | `task/T5` `31f470f` | **accepted** (9) | UHPC Stage 3 part a: split V_UHPC,i = V_i(ρ_i − ρ_air)/(ρ_UHPC − ρ_air), shells at t_min in k*, ballast for flotation, F9 properties on the exact body, F10 check, S8 record, `final_props` = realised design; pre-check D16 deleted |
| T7 | `task/T7` `4c310df` | **accepted** (8 + quick check) | thin-shell Stage 3 rebuilt: one t, `z_ballast` free across module edges, flotation and GM = GM_Stage2 equalities, item-27 objective, draft last, t_max = d_close − ε/2, closest fail |
| T8 | `task/T8` `3b0a0be` | **accepted** (9) | figures from the realised solid: `realised_section_data` (F12) with exact y = 0 crossings, void outlines with holes, joint steps, ballast level; `plot_modular_precast`, `plot_steel_solve`, `plot_optimised_cross_section` |
| T10 | `task/T10` `07856a2` | **accepted** (9) | `export_stage3` (F11): UHPC per-module and fused STEP, steel ballast/shell/combined STEP |
| T2a | `task/T2a` (head §9) | round 2 at 5, round-3 fixes in progress | exact geometry kernel: F1 `outer_nurbs`, F2 `offset_surface` + `trim_fold` + `fit_bspline_surface` (normal offset, fold trimming, adaptive fit, metrics M1–M3), F2b `void_closing_distance`, F3/F3b/F4 B-spline evaluate/split/slice |
| T3 | `task/T3` (head §9) | implementing (base T2a; merges T2a's fixes) | F5 `build_body`, F6 `body_properties`, F6b `body_section`, F7 `hydrostatics_at_draft`; C1 tests; shared C1 kernel cache `tests/c1_kernel_cache.m` (in progress) |
| T6 | `task/T6` (head §9) | round 1 at 5, fixes in progress (base T5) | UHPC Stage 3 part b: optimisation from the split, spill, draft last, closest fail, report |
| T4b | `task/T4b` (head §9) | round 1 at 6, fixes in progress (base T4a) | Stage-2 floors from the kernel (F8 `density_floors`), floors in Stage 1 and 2 for both modes; deletes `compute_perpendicular_shell_volume` (D1) |
| T2b | — | not started | general path (`fit_z_faces`) for hull patches that do not convert exactly (not needed for C1) |
| T11–T13 | — | not started | cleanup, docs (`METHODS_ENGINE`, `RUNTIME_GUIDE`, `RESULT_SCHEMA`, code map), final review |

Every checker verdict is saved in `docs/plans/orchestration/pause/` (`<task>_verdict_r<n>.json`:
summary, required fixes with file, line, issue, fix and how to verify), next to the implementers'
reports and `launch_inputs.json` (every task's brief).

### 2.3 What already works on C1 (measured on the real kernel, in branch joins)

- Normal-offset shells, checked with the independent Python evaluator on 16 000 random points per
  thickness: t_local 76.208–76.917 mm at t = 76.2 mm and 25.419–25.654 mm at t = 25.4 mm, inside the
  fitting band t ≤ t_local ≤ t + ε (ε = 0.01 t_min; the shell is built at t + ε/2).
- Neck closing distance d_close = 0.100 m.
- Stage-2 floors from the kernel (T4b), kg/m³, module 1 at the keel:

  | Mode | 1 | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|---|
  | UHPC, t_min 76.2 mm | 610.4 | 285.7 | 212.1 | 612.8 | 2500 (wall module) |
  | Thin shell, t_min 25.4 mm | 725.1 | 269.2 | 454.6 | 2018.2 | 2228.7 |

  v1.0's circle floors were 437 / 227 / 175 / 369 for the UHPC hollow modules. V_solid + V_air = V
  per module exactly. Hull volume 21.1784 m³ on the exact surfaces (hydrostatic table: 21.1608 m³).

## 3. Working locally in MATLAB

1. `git clone https://github.com/VKotzamanis/MWEC_Mass.git`, `git checkout claude/lucid-cray-7o9442`,
   `git fetch origin 'refs/heads/task/*:refs/remotes/origin/task/*'`.
2. Run the suite from the repository root: `run('tests/run_tests.m')` at the prompt, or
   `matlab -batch "run('tests/run_tests.m')"`. The environment variable `TESTS_FILTER` restricts it,
   e.g. `setenv('TESTS_FILTER','tests/solid')`. Under MATLAB the Octave shims stay off the path.
   The tests were written and run under Octave 8.4; this is their first MATLAB run, so a failure can
   be a MATLAB/Octave difference in the test itself.
3. Production code targets MATLAB R2022b. Figures were never rendered in the container (Octave
   lacks `tiledlayout` and others); their data (F12) is tested, the drawing is yours to confirm.
4. `build_config` cache: `Output/cache/` (or `in.files.geometry_cache_dir`). It rebuilds by itself
   when the deck, a keyed input or a keyed source file changes; delete it to force a rebuild.
5. Python is needed only for `tests/step_check.py` (gmsh 4.15) and the reference evaluator.
6. The pipeline regressions (`MWEC_REGRESSION=1`) compare with Octave baselines that are stale by
   design (regenerated once at J2); do not use them as a gate before then.

### Agents (Claude Code, local)

`CLAUDE.md` loads `AGENTS.md` into every session. `.claude/agents/` defines three subagents:
`mwec-implementer` (Opus, high: implements or fixes one task on its branch), `mwec-checker` (Opus,
xhigh, read-only: strict review; marks findings blocking or minor) and `mwec-quick-check` (Sonnet,
medium, read-only: confirms minor fixes). Loop per task: implementer → checker → implementer fixes
→ checker again for blocking findings, quick check for minor-only findings → merge. One MATLAB
test run at a time (or `-singleCompThread` per run, never more runs than cores).

## 4. Merge plan (what remains, in order)

1. **Finish T2a** (§5) and **T3**, then **J1**: merge `task/T2a`, then `task/T3` (T3 already holds
   T2a's commits). For C1, T2b is not needed: C1 takes the exact path everywhere. At J1:
   - retire or adapt the 6 tests in `tests/standin_kit` that fail only because the real F1–F4 now
     shadow the stand-ins (they pass at the base);
   - **owner checkpoint**: real C1 UHPC XZ (y = 0) and YZ (x = 0) sections from the real kernel with
     the rule-12 metrics beside them (every drawn boundary has a face of the solid behind it, drawn
     areas equal the body's section areas, the drawn ballast level equals `z_ballast`).
2. **T4b** after J1 (it needs the real kernel on main).
3. **J2**, in this order: T5, T6, T7, T8, T10. Each task recorded "deferred" edits to shared files it
   could not edit on its base; apply them as that task's last commit when it merges (exact files,
   lines and text in its last report). Known ones:
   - T5: `optim/report_assemble_results.m` (drop `constructability`, `steel_data`; add
     `results.stage3 = []`), `optim/run.m` (drop the matching variables), `build_config.m` lines
     325–334 (`uhpc_t_init`, `uhpc_max_slope_factor`, `uhpc_n_z_grid`), `WEC_User_Input.m` 53–56.
   - T7: d1–d3, d6 in `WEC_User_Input.m` (43–46, 53–56: one thin-shell t_min input, t_init removed)
     and `build_config.m` (275, 277–278, 306–309), plus its test lines 33–35 and 187–188.
   - T8: `output/dispatch.m` line 24 becomes
     `mwecmass.output.figures.plot_optimised_cross_section(final_props, config, x_opt, results.stage3);`
     (`config.shell = []` in `build_config.m` goes to T11; `_graph/CODE_MAP.md` to T12).
   - T10: the `validate_save_flags` edit in `dispatch.m`.
   - T6: the `tests/octave_shims/fmincon.m` change (empty multipliers when Octave's QP fails) is
     test-only and accepted.
   Line numbers refer to each task's base; check them against the merged file before applying.
   After J2: the first whole-pipeline runs of both modes on C1 (Stage 1 → Stage 3, STEP, figures),
   the Octave baseline regenerated once, owner checkpoints after T6 and T10.
4. **G**: T11 (cleanup: the `config.shell` assignment; dead files `empty_realised_properties.m`,
   `thin_shell/integrate_split.m`, `build_realised_properties.m` per contract §4), T12 (docs,
   `METHODS_ENGINE.md` lines 205–208, code map, `tests/README.md` folders table; delete this file and
   `docs/plans/orchestration/`), T13 (final review).

## 5. Open work per task

**T2a** (round 2: 5; the C1 results are correct, general hulls not yet):
- F1 must give every patch end row the height MS2Parser gives exactly (C1 keel −3.25, not
  −3.2500000000000004), so bitwise module-edge comparisons never cut a sliver.
- F2b misses a void that closes at the hull's own end face without a double normal (test deck: an
  inverted frustum); it must find opposite offset layers that meet there.
- Minor: a comment calls the inner profile's largest axis distance "the void's half-width".
- The round-1 fixes (general decks reach `fit_z_faces`, decimal arcs, fold trimming on a fillet
  deck, derived rounding bounds, M3 at every check height) are in and were verified in round 2.

**T3**: F5–F7 committed; remaining: C1 tests through the shared kernel cache, volume closure per
module and for the hull with the ballast inside module 2, at a module edge and spilled; C1 CG
x = y = 0; `step_check` import of the C1 STEP.

**T6** (round 1: 5). The optimiser must:
1. establish both equalities (flotation, GM = GM_Stage2) before optimising (e.g. Newton steps on the
   two equalities, or a restart from the best feasible evaluated point);
2. keep phase-2 progress whenever it holds the equalities (it was discarded);
3. release the draft only when no point at the Stage-2 draft can hold flotation (item 31), never
   because one rebuild failed;
4. call F7 once while the draft is fixed (contract §7 item 3);
5. have tests that assert the outcome: case C must reach a design at least as good as the
   checker's feasible one (objective 0.0084, T_pitch deviation 9.2 %); case D must float (it ended
   submerged, flotation residual 0.29, although the lightest design floats at vs = 0.937 m).

**T4b** (round 1: 6; the floors are correct):
- thin-shell Stage 2 with the new floors crashes under Octave at the bottom-filled start (Octave
  `sqp` returns empty multipliers when its QP subproblem is unbounded; MATLAB's `fmincon` is not
  affected): Stage 2 must treat such a start as failed, as Stage 1 already does;
- solid modules must return ρ_solid exactly (assign it; (ρ·V)/V is not always ρ in floating point);
- rule 10: in thin-shell mode do not store the shell density in a field named `ballast`;
- minor: a docstring claim; Stage-1 `sweep` mode still screens without the floors.

**T2b**: not started (`fit_z_faces`, general path, mirrors); needed only for hulls that are not
exact-path like C1.

## 6. Decisions not yet written into the contract (apply them)

- F10 `pass` = the four metrics within `mass_acceptable_pct` and flotation; the GM equality residual
  is reported beside the check without entering it (AGENTS §3 item 27, §5 item 11).
- F11 `out_dir` is the type folder; `export_stage3` writes into its `step` subfolder.
- Thin-shell floors live in `config.per_strip_density_lb`, like the UHPC ones.
- t_max,i per UHPC module = d_close over the module's whole range − ε/2 (fixed, conservative);
  inner sets are keyed by t and z-range.
- Kernel errors (`VoidClosed`, `JointNotNested`, `FitNotConverged`) are failed evaluations in the
  Stage-3 optimisation, never an abort.
- Closest fail: points that hold flotation rank first, by the item-27 objective; if none holds it,
  the smallest equality violation wins.
- Stage-2 start choice: a start that holds the constraints beats one that does not; then the lower
  objective; if none holds them, the smallest violation (owner agreed).
- A module whose t_min shell holds more UHPC than its V_UHPC split is reported in the Stage-3 log.
- Kernel: patch end rows take the parser's exact heights; on the exact path F2 offsets the F1 NURBS
  as written; F1 owns the classification rules 1a/1b and `fit_z_faces` makes the mirrors; rounding is
  bitwise wherever the construction makes values identical, otherwise a bound derived in one comment
  line from the coordinate magnitudes and the operation count.
- Fold trimming is tested on decks that fold; C1 cannot fold before its neck closes (smallest convex
  radius 0.100 m = d_close), so C1 tests that d ≥ 0.1 raises `VoidClosed`.
- From the stand-in kit: a void end counts as a flat part for F6b `side`; the bad-side error is
  `mwecmass:solid:BadSide`; a row cut at a plane takes that height bitwise.

## 7. Octave-only behaviour (ignore in MATLAB)

- Octave's `sqp` is weaker than MATLAB's `fmincon` SQP: v1.0's precast Stage 2 stops at exitflag −2
  (vs 0.409 against MATLAB's 0.983), and `sqp` can return empty multipliers (the shim handles it).
- Octave builds C1's geometry in 6–7 min.
- One Octave process per core was enforced in the container with a wrapper (owner rule,
  `docs/plans/orchestration/octave_one_core.sh`); not needed locally.

## 8. How the work was run (only if agents continue it)

Scripts: `docs/plans/orchestration/lane.js` (implementer → Opus xhigh checker → fixer; minor findings
go to a quick check instead of a full round; targeted tests while working and the whole suite once
before acceptance; one Octave per core; judge by metrics, rule 12) and `merge.js` (merge `--no-ff`,
suite, push, checker). Lessons: commit and push at least every 30 min (the container restarted about
every 4 h); report measured numbers to the owner, never pictures of test shapes; verify a claim on
the code and the exact deck before writing it into a spec.

## 9. Branch heads at the pause (2026-10-08 ~01:40 UTC; all pushed, nothing running)

| Branch | Head | Last commits since its last verdict |
|---|---|---|
| `claude/lucid-cray-7o9442` (main) | this commit | T0, T1, T9, T0c, T0b, SK, SK2, T0d, T4a merged; MATLAB-safe `run_tests.m` |
| `task/T5` | `31f470f` | accepted |
| `task/T7` | `4c310df` | accepted |
| `task/T8` | `3b0a0be` | accepted |
| `task/T10` | `07856a2` | accepted |
| `task/T2a` | `270d0a6` | `57e87e0` keel row at −3.25 bitwise (fix 1 of round 2, done); `270d0a6` unverified work in progress on the F2b end-face closing (fix 2) — review before keeping |
| `task/T3` | `3b9fb6b` | `tests/c1_kernel_cache.m` added (folder `/home/user/geomcache`; set the environment variable `MWEC_GEOMCACHE` to a local folder) |
| `task/T6` | `d56c817` | outcome assertions for cases C, S, D, G and a known-optimum test (`81bff4c`, written but not yet run), keep the step's start when the chosen point ranks worse (`2b67fef`), solve returns the evaluation behind `solver(end)` (`d56c817`); not yet graded |
| `task/T4b` | `1023873` | round-1 fixes: shim survives a core-sqp QP failure, solid-module floors set by definition, thin-shell floors take the real ballast density, Stage-1 screen respects the floors (`a359046`), shim test (`1023873`); not yet graded |

Commits titled "work in progress ... (unverified)" hold an agent's edits saved when it was stopped
(container restart or this pause); review them before building on them.
