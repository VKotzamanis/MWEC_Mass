# HANDOFF — MWEC_Mass exact-geometry refactor

Updated 2026-10-07 16:40 UTC on `claude/lucid-cray-7o9442`. Process file for agents: the
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

Main `claude/lucid-cray-7o9442` = `f52e232`: general-kernel amendment, AGENTS rule 12 (judge by metrics), SK2, T0d; suite 27 passed, 1 skipped (476 s). Every task branch `task/<id>` is on GitHub.

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
| T5 | UHPC Stage 3 a: split, build, check, store | U | accepted (round 5); waits for J2 | 7, 8, 8, 8, 9 | `31f470f` |
| T6 | UHPC Stage 3 b: optimisation, spill, closest fail | U | implementing (`wf_e905b818-430`, base task/T5) | — | task/T6 |
| T7 | thin-shell rebuild | S | resumed after restart (`wf_1ef512c2-5e3`), rounds 6–7 | 7, 7, 7, 8, 8 | task/T7 `4ae74bb` |
| T8 | figures from the realised solid | O | restarted (`wf_69d3ee54-efc`), Opus fix of the round-8 verdict, rounds 9–10 | 8, 8, 8, 8, 6, 6, 6, 6 | task/T8 `b70b1ff` |
| T10 | Stage-3 STEP exports | O | accepted, waits for J2 | 9 (3) | `07856a2` |
| T0d | geometry cache of `build_config` (C1: fresh build 421 s CPU, reload 0.44 s) | P | merged `f52e232` | 9 (2); merge 10 | `158d77d` |
| T4a | Stage-2 changes (delete `c_mono`, `c_mass_min`; bottom-filled start; bounds) | P | implementing, restarted (`wf_6bfb8540-04d`, then T4b) | — | task/T4a `7ad75d5` |
| T4b | Stage-2 floors from the kernel (merges after J1) | P | queued after T4a | — | — |
| SK2 | stand-ins updated to spec2 | SK2 | merged `1bb06d6` | 9 (3); merge 10 | `753a11e` |
| T2a | exact-path offset, fold trim, adaptive fit | K1 | implementing, restarted (`wf_6f9b1ce1-590`): F3, F3b/F4, F1 and part of F2 committed | — | task/T2a `c825835` |
| T2b | general path (refits, flat regions, mirrors) | K | after T2a, ∥ T3 | — | — |
| T3 | bodies and exact properties | K | after T2a | — | — |
| J1 | merge T2a, T2b, T3; owner checkpoint | — | pending | — | — |
| J2 | merge T5, T6, T7, T8, T10 (in order); first whole-pipeline runs; owner checkpoints after T6, T10 | — | pending | — | — |
| G | T11 cleanup → T12 docs → T13 final review | — | pending | — | — |

## 3. How the work is run

- **Metrics decide (AGENTS rule 12).** Every acceptance, rejection and report to the owner cites
  numbers computed from the code's outputs. Never send the owner a picture as evidence; send the
  failing metric. Figures for the owner: real C1 XZ (y = 0) and YZ (x = 0) sections from J1 on, with
  the metrics printed beside them.
- The orchestrator delegates; implementers and graders are subagents run by the Workflow tool.
  Scripts: `docs/plans/orchestration/lane.js` (with rule 15 and grader item 7, judge by metrics) and
  `merge.js`. New lane launches use the session copy `workflows/scripts/mwec-lane-v2.js` (a task with
  `resume: {last, r0}` starts from its last verdict: fix, then two grading rounds); resume a
  workflow started before 13:50 UTC only with its own script `mwec-lane-wf_d681e461-e66.js`
  (resume reuses cached agent calls only while the prompts are unchanged).
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
- **Owner (16:35 UTC): one Octave process per core, never shared, never spread across CPUs.**
  `docs/plans/orchestration/octave_one_core.sh` is installed as `/usr/bin/octave` (original moved to
  `/usr/bin/octave.real`): each run takes a free core (lock `/tmp/octave-core-<k>.lock`, held until it
  exits), is pinned there with `taskset`, runs single-threaded BLAS, and waits while all 4 cores are
  busy. Check it after a container restart (`taskset -cp` on a running `octave-cli`). Lane rule 16.
  Bitwise tests recorded under multi-threaded BLAS: if one fails only now, check the thread count first.

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
- Orchestrator decisions of 14:40 UTC (contract errata to apply at J2 through an author agent and
  grader): F10 `pass` = the four metrics and flotation; the GM equality residual is reported beside
  the check (AGENTS item 27 and §5 item 11 clarified on main from the owner's items 4.4 and 32; the
  SK stand-in F10 still counts the GM row and must be aligned). F11 `out_dir` = the type folder
  (`export_stage3` adds `step`). Thin-shell floors live in `config.per_strip_density_lb` (T4b).
  t_max,i per UHPC module = F2b on the module's hollow z-range − ε/2 (contract "Bounds" line).
  Kernel errors (VoidClosed, JointNotNested, FitNotConverged) are failed evaluations in T6. F6
  integral caching is internal to T3 (F5/F6 signatures unchanged). A module whose t_min shell
  exceeds its split is reported in the Stage-3 log (no S8 field).
- SK2 spec issues for J1/T3 (contract errata with the batch above): whether a void end (constant-z
  end piece of an inner set) counts as a flat part for F6b `side` (stand-in: yes); name the error
  for a bad `side` (stand-in `mwecmass:solid:BadSide`); a row cut at a plane takes that height
  bitwise (F3b, F5; the stand-in does); J1 adapts `test_sk_outer_nurbs` (asserts that
  `force_general` is ignored). The SK stand-in F10 still counts the GM row in `pass` (align at J1).
  T5 edited `tests/standin_kit/test_sk_not_analytic.m` and `test_sk_realised.m` (merge with SK2 clean).
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

- Container restarts (00:18, 04:02, 15:54 UTC) kill every background job and workflow; only committed and
  pushed work survives. After a restart: recover each run's launch input from the session transcript,
  check every worktree (delete `octave-workspace` crash dumps), and relaunch each lane from its last
  verdict with a generated `mwec-lane-<lane>-restart.js` (args embedded) and a note to continue from
  the committed state. Workflow resume reuses only the longest unchanged prefix of agent calls;
  editing an early prompt reruns everything after it.
- Whole-pipeline runs in Octave (24–94 min) die in restarts; the Octave `fmincon` stand-in fails the
  v1.0 modular Stage 2 (exitflag −2). Do not compare Octave with MATLAB v1.0 numbers (owner).
- Spec errors made by the orchestrator: the 40 mm shoulder radius (sampling artefact), the 80 %
  height band, the same file assigned to two tasks, KG used for Z_CG, an incomplete quick fix of the
  amendment. Verify claims on the code and the exact deck before writing them.
- Deleting remote branches returns HTTP 403 (policy): do not retry.
- 15:12–15:21 UTC: every `git push` and GitHub API write returned HTTP 500 (reads worked,
  githubstatus.com reported no incident); it cleared by itself. `scratchpad/push_retry.sh` retries
  main and fast-forwards lagging task branches. A grader that requires origin = HEAD may reject only
  for an outage like this; discount such rejections.
- Reading a workflow's result from the notification text can mislead; read its `journal.jsonl`.
- Hand edits of the dense contract by the orchestrator were rejected twice: route every contract
  change through an author agent and the grader.
- T8 round 4: I sent the owner rendered pictures of stand-in figures instead of the metric the
  grader had already measured: dashed void boundaries with no face of the solid behind them at the
  3 interior module edges of the thin-shell cylinder stand-in. The owner had asked from the start
  for judgement by metrics; the handoff did not state it as a rule. It is now AGENTS rule 12.

## 7. Next steps

1. When T2a is accepted: launch T3 and T2b in parallel from `task/T2a` (lane script, mode parallel).
2. As lanes return: merge T0d and T4a when accepted; T4b after J1; T5–T10 wait for J2.
3. **Owner (14:50 UTC): "When you get the first REAL cross sections for UHPC, show me the figure and
   pause your work."** The first is the J1 checkpoint: a C1 UHPC body built by the real kernel
   (T2a, T2b, T3), XZ (y = 0) and YZ (x = 0) sections with the rule-12 metrics printed beside them.
   Then stop every running workflow (work is pushed per commit), launch and merge nothing, and wait
   for the owner's reply.
   J1, owner checkpoint after T3 (measured numbers), then J2 with the first whole-pipeline runs,
   the baseline regenerated once, owner checkpoints after T6 and T10; then G.
4. Session resources: the 7-day rate limit was at warning level at 12:20 UTC (reset 21:00 UTC);
   context 688k of 1M used.
