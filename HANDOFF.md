# HANDOFF — MWEC_Mass exact-geometry refactor

Updated 2026-10-07 13:45 UTC on `claude/lucid-cray-7o9442`. Process file for agents: the
orchestrator updates it at every milestone and at least every 100k tokens; delete it in T12.

## 1. Goal

Make the MATLAB code correct per the owner's rules: every geometric operation from the `.ms2`
parametric definition (no shape assumptions), Stage 3 realised from the whole Stage-2 solution on
one exact geometry kernel (UHPC modular precast and steel thin shell as separate pipelines), and
STEP output of the realised bodies. The binding documents, in this order:
- `AGENTS.md` — owner rules (§1), current code (§2), owner corrections and decisions (§3 items
  1–36), issues (§4), expected outcome (§5), decisions log (§8).
- `docs/plans/2026-10-07-interfaces.md` — interface contract: structs and functions (§1–2),
  ownership and stand-ins (§3), shared-file order (§4), invariants (§5), lanes and joins (§6).
  It wins over the plan where they differ.
- `docs/plans/2026-10-06-exact-geometry-stage3-step.md` — task sections T0…T13, deletion register.

## 2. Progress tracker

Main `claude/lucid-cray-7o9442` = `e7e9114` (pushed; holds the general-kernel amendment). Every task branch `task/<id>` is on GitHub.

| Task | What | Lane | Status | Grade (rounds) | Head |
|---|---|---|---|---|---|
| T0 | test harness, Octave shims, baselines | A | merged | 9 (3) | `c2234c7` |
| T1 | exact outer rows and normals | A | merged | 9 (3) | `85443fb` |
| T9 | AP214 STEP writer, `tests/step_check.py` | A | merged | 9 (2) | `53067fd` |
| spec | owner decisions, plan, contract | — | merged `7fcd02d` | 9 (5) | `e60d6d7` |
| T0c | parser resolves entities once (4.3× faster, bit-identical) | N | merged `a935491` | 9 (2) | `80ef938` |
| SK | stand-in kit (box, cylinder, closed forms) | K | merged `8eb4552` | 9 (3) | `a03bc4b` |
| T0b | renames (`z_ballast`, `rho_air`, `rho_ballast`, `uhpc`) | N | merged `b4523aa` | 9 (4) | `c64a4d0` |
| spec2 | contract amendment: general kernel (owner) + errata | K | merged `e7e9114` | 9 (round 8; errata 9; merge 10) | `5ce48fc` |
| T5 | UHPC Stage 3 a: split, build, check, store | U | fixing | r1 –, r2 8 | task/T5 |
| T6 | UHPC Stage 3 b: optimisation, spill, closest fail | U | queued after T5 | — | — |
| T7 | thin-shell rebuild | S | fixing | r2 7 | task/T7 |
| T8 | figures from the realised solid | O | extra fix rounds 5–6 (void outlines at module edges) | r1–r4 8 | task/T8 |
| T10 | Stage-3 STEP exports | O | accepted, waits for J2 | 9 (3) | `07856a2` |
| T0d | geometry cache of `build_config` | P | fixing | r1 8 | task/T0d |
| T4a | Stage-2 changes (delete `c_mono`, `c_mass_min`; bottom-filled start; bounds) | P | queued after T0d | — | — |
| T4b | Stage-2 floors from the kernel (merges after J1) | P | queued after T4a | — | — |
| SK2 | stand-ins updated to spec2 | SK2 | implementing | — | task/SK2 |
| T2a | exact-path offset, fold trim, adaptive fit | K1 | implementing (Opus high) | — | task/T2a |
| T2b | general path (refits, flat regions, mirrors) | K | after T2a, ∥ T3 | — | — |
| T3 | bodies and exact properties | K | after T2a | — | — |
| J1 | merge T2a, T2b, T3; owner checkpoint | — | pending | — | — |
| J2 | merge T5, T6, T7, T8, T10 (in order); first whole-pipeline runs; owner checkpoints after T6, T10 | — | pending | — | — |
| G | T11 cleanup → T12 docs → T13 final review | — | pending | — | — |

## 3. How the work is run

- The orchestrator delegates; implementers and graders are subagents run by the Workflow tool.
  Scripts: `docs/plans/orchestration/lane.js` and `merge.js` (copies of the scripts in use).
- **Lane workflow** (`lane.js`, args `{lane, mode: 'chain'|'parallel', tasks: [{id, wt, branch,
  base, model, effort, brief}]}`): per task, implementer → Opus 5.5 xhigh grader → fixer, at most 4
  grading rounds. The grader accepts only with no rule violation, every acceptance item reproduced
  by itself, all tests passing, its own exact test written, and the minimum of five scores ≥ 9; it
  reports `spec_issues` instead of enforcing a wrong spec. Implementers commit after every
  deliverable and push their task branch after every commit. Workflows never merge.
- **Merge workflow** (`merge.js`, args `{branches: [...]}`): merge `--no-ff`, run the suite, push
  only if green; an Opus xhigh grader checks diff = union of the branches, reruns the suite, checks
  the remote. One merge at a time.
- Models: Opus 5.5 high for kernel and Stage-3 tasks (SK, T2a, T2b, T3, T5, T6, T7, contract
  edits); Sonnet 5.5 high for T0b, T0c, T0d, T4a, T4b, T8, T10. Graders: Opus 5.5 xhigh.
- Worktrees `/home/user/wt/<task>` on `task/<task>`. A chain's next task branches from the
  previous accepted branch. Shared files (contract §4): a task edits one only when the earlier tasks
  in its row are merged into its base; otherwise it records the exact edit as `deferred`, applied as
  its last commit at J2.
- Testing: GNU Octave 8.4 (`bash tools/install_toolchain.sh` in a fresh container), suite
  `octave --no-gui --quiet tests/run_tests.m` (24 pass, ~210 s). Owner: no whole-pipeline run before
  J2; `MWEC_REGRESSION=1` enables the pipeline regression (24–94 min per mode in Octave).
- A workflow runs at most 2 agents at once (4 CPUs − 2); run several workflows for parallelism.

## 4. Decisions not yet in the merged contract

- spec2 (general kernel, owner: "No it needs to be generalized"): exact path where an entity
  converts exactly and z depends on one parameter (F1 orders it as u), monotone; otherwise faces
  refitted through exact points with z as a parameter, split at z-turning rows and C0 seams,
  untrimmed; the only number is ε = 0.01·t_min (fit share ε/4 derived); one closed loop per section.
  Orchestrator decisions O1–O7 (OD3 wording; F1 swaps u/v; T2 → T2a + T2b; SK2; shared-edge
  equality = same degree, bitwise control points and weights, knots equal after the affine map to
  [0, 1]; column splits at patch-corner heights; split_wall_box scope). Parser facts: unknown entity
  types are skipped silently and an unknown mirror plane is read as y = 0 with a warning.
- Contract errata to apply right after spec2 merges (SK grader): F3 exempt from `NotAnalytic`; F9
  stand-in identifies the fixture by `geo.analytic` or `hull_name`, F10 by `props.analytic`; a file-map
  row for `tests/standin_kit/*` (J1, J2 owners); thin shell: `z_ballast` equal to the inner `z_lo`
  cuts the inner set (no zero-thickness layer); S4: precast joint faces take precedence over inner
  ends on a module edge; `tests/README.md` folders table → T12.
- Carry-overs: `build_config` time after T0c (T0d prints it); `METHODS_ENGINE.md` lines 205–208 and
  `_graph/CODE_MAP.md` regeneration → T12; T0b deferrals for T5 (`rho_UHPC`/`t_UHPC` capitals,
  `config.rho_steel` in precast `solve.m`). The old names stay only in `tools/`,
  `tests/baseline/*.json`, `Output/`, `docs/plans/` and AGENTS' rename statements.

## 5. What worked

- The grader loop finds real defects every round (e.g. a STEP checker that dropped sheet bodies,
  NaN directions in STEP, a parser fixture with absolute paths, unbuildable stepped decks).
- Contract plus stand-ins let lanes U, S, O, P run in parallel against closed-form bodies.
- Checking claims on the exact deck: the C1 curvature figure (smallest convex radius 0.100 m; no fold
  at 25.4 or 76.2 mm) corrected a spec error.
- Pushing main only after a green suite, or when the merged tree equals an already graded tree.

## 6. What did not work (do not repeat)

- Container restarts (00:18, 04:02 UTC) kill every background job and workflow; only committed and
  pushed work survives. Workflow resume reuses only the longest unchanged prefix of agent calls;
  editing an early prompt reruns everything after it.
- Whole-pipeline runs in Octave (24–94 min) die in restarts; the Octave `fmincon` stand-in fails the
  v1.0 modular Stage 2 (exitflag −2). Do not compare Octave with MATLAB v1.0 numbers (owner).
- Spec errors made by the orchestrator: the 40 mm shoulder radius (sampling artefact), the 80 %
  height band, the same file assigned to two tasks, KG used for Z_CG, an incomplete quick fix of the
  amendment. Verify claims on the code and the exact deck before writing them.
- Deleting remote branches returns HTTP 403 (policy): do not retry.
- Reading a workflow's result from the notification text can mislead; read its `journal.jsonl`.
- Hand edits of the dense contract by the orchestrator were rejected twice: route every contract
  change through an author agent and the grader.
- Owner: figures are judged by metrics on the section data (every drawn boundary has a face of the
  solid behind it, drawn areas equal the body's section areas, ballast level within contract I3), not
  by pictures of stand-ins. The owner reviews real XZ (y = 0) and YZ (x = 0) cross sections of the C1
  bodies, from J1 on (real kernel), with the metrics printed beside them.

## 7. Next steps

1. When T2a is accepted: launch T3 and T2b in parallel from `task/T2a` (lane script, mode parallel).
2. As lanes return: merge T0d and T4a when accepted; T4b after J1; T5–T10 wait for J2.
3. J1, owner checkpoint after T3 (measured numbers), then J2 with the first whole-pipeline runs,
   the baseline regenerated once, owner checkpoints after T6 and T10; then G.
4. Session resources: the 7-day rate limit was at warning level at 12:20 UTC (reset 21:00 UTC);
   context 688k of 1M used.
