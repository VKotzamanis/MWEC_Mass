# Owner decisions of 2026-10-07 (apply to AGENTS.md and the plan after Phase A merges)

1. Delete c_mono: stage2_constraints.m (c_mono, its fallback size in the catch block) and the copy
   c_monotonic in optim/solve_2d_surrogate.m (lines ~226-269, fallback size). Owner's reason it
   existed: the optimiser did not put material low first and returned un-optimised results.
   Safeguard announced to the owner: Stage 2 gets a second start from a bottom-filled point at the
   Stage-1 draft (every module at its floor, then modules filled solid from the keel up until
   mass = rho_w * V_sub); keep the start with the lower objective; log both.
   Stability stays enforced by GM >= gm_min (Stage 2) and GM = GM_Stage2 (Stage 3).
2. Stage-3 objective (both modes, consistency): come closer to the Stage-2 converged solution.
   Minimise sum of ((X3 - X2)/X2)^2 over Z_CG = CG_total(3), coupled T_heave, coupled T_pitch.
   Equalities: flotation (exact, solver tolerance) and GM = GM_Stage2. The 10 % check
   (mass_acceptable_pct) decides pass/fail after every solve. Supersedes OD1 / section 3 item 15
   (range penalties) for Stage 3. Stage 1 and 2 objectives unchanged.
3. UHPC: t of the hollow shell above the ballast in the ballast module k* is a variable >= t_min
   (rule 4), start value t_min in the split. Unknowns at fixed draft: ballast level, t_k*, t of
   each hollow module above k*. Wording: "wall module" = the solid module only; the UHPC layer
   around air in a hollow module is the "shell".
4. Thin shell t_max = half the thickness of the slender wall (neck), where the two offset plates
   meet; from the exact geometry; 0.10 m for C1. Delete the 10 %-90 % height band
   (thin_shell/solve.m:77-80) and the plan's "middle 80 %" sentence (plan line 248).
5. Delete the UHPC pre-check with the 0.95/1.05 factors (modular_precast/solve_and_extract.m:95-106):
   invented tolerance and abort instead of closest fail. Add to the deletion register (T5).
6. Delete c_mass_min (implied by the floor bounds) in T4; told to owner, no objection.
7. Draft released only when mass balance cannot be met at the Stage-2 draft (owner's rule).
8. OD10 by owner's wording ("check if mass balance is satisfied", "Z_CG, GM and periods ...
   within the acceptable bound"): flotation exact; mass_acceptable_pct for Z_CG, GM, periods.
9. Dropped (owner): table-vs-exact hydrostatics note. Do not raise again.
10. Owner approved (2026-10-07): push every task branch (task/T*, task/spec) to GitHub as a backup;
    delete the remote task branch once it is merged into claude/lucid-cray-7o9442. Future workflow
    RULES: implementers and fixers push their own task branch after each commit (never any other
    branch). Phase A: background loop re-pushes task/T0, T1, T9 every 5 min (stop it after merge).
11. Spec branch task/spec (f51f760, accepted 9/10) merges right after Phase A integration.
12. Owner (2026-10-07): parser resolves each entity once (T0c) and build_config geometry products
    are saved/reloaded with a fingerprint (T0d). Waves revised: B = T0c, T0d, T2; C = T0b, T3;
    D = T4a then T4b. Full-pipeline regression opt-in via MWEC_REGRESSION=1 (in T0d).
    Spec branch task/spec now at the commit after 64f53ac (pushed).
13. Owner (04:3x): drop T0e ("You can't expect to match my matlab results when we're running
    different geometry and initilization"); no whole-pipeline run before Waves B-E are
    implemented; run many subagents in parallel. Plan: interface contract first (body = closed
    B-rep per module, T9 struct), then lanes in parallel. Grader improvements: deferred-findings
    ledger, grader-written exact test, contract conformance, brief check.
    Workflow facts: concurrency per workflow = min(16, CPUs-2) = 2 here; resume reuses only the
    longest unchanged prefix of agent() calls.
14. STATE ~08:xx: Phase A merged (d93b52e) + README fix (25b7706); spec + contract merged (7fcd02d, pushed).
    Remote branch deletion is blocked by policy (403): do not retry; owner informed.
    Phase B started: workflow war84khjb (SK opus-high, T0b sonnet-high, parallel) and w2vcxr2ml (T0c).
    Generic lane script: workflows/scripts/mwec-lane-wf_d681e461-e66.js (args: lane, mode chain|parallel,
    tasks[{id,wt,branch,base,model,effort,brief}]). Workflows never merge; I merge one at a time with an
    integration agent + grader. Lane plan (contract section 6): after SK+T0b merged -> K (T2->T3, opus high),
    U (T5->T6), S (T7), O (T8 || T10); after T0c merged too -> P (T0d->T4a->T4b; T4b merges after J1).
    J1 = merge T2+T3 (owner checkpoint). J2 = merge T5,T6,T7,T8,T10 chain after J1 and lane P, then the
    first whole-pipeline runs, owner checkpoints after T6 and T10. G: T11->T12->T13.
    Open for owner: contract section 9 item 1 (kernel only for the C1 class, explicit errors otherwise) -
    asked, no answer yet.
15. Owner: "No it needs to be generalized." Kernel general (AGENTS 5.9 rule: exact NURBS where possible,
    else fitted with the same metrics). Amendment workflow w60dwvclg on task/spec2 (worktree wt/spec2):
    general path = z-parametrised refits through exact points, split at z-turning rows and C0 seams,
    untrimmed faces; only number eps = 0.01 t_min; scope one closed loop per section; T2 gets a non-C1
    test hull. Merge task/spec2 before lane K (T2) starts; SK-visible changes must be listed.
16. T0c accepted 9 (80ef938): 4.3x faster fixed workload; bit-identical over 1,589 calls. Merge workflow
    wuy9k0c77 running. Carry-overs: T0d brief must print build_config time (fresh build, new parser) as
    the 'after' value (T0 'before': 2,067 s under profiler); METHODS_ENGINE.md 205-208 and _graph/CODE_MAP
    regeneration -> T12. Contract amendment must require an explicit error for .ms2 entity types the
    parser cannot evaluate (today silently dropped) and assign the MS2Parser change to a task (file map:
    MS2Parser.m was T0c only) - check when w60dwvclg completes.
17. T0c merged + verified 10 (a935491). General amendment run 1 ended unaccepted (6): flat regions,
    stepped exact patches, F1 deck check too strict, SK-visible changes understated, mirror vertex union,
    OD3 stale, u/v wording. Run 2 = wdsp7ry0w (rounds 4-6) with O1 OD3 edit, O2 F1 swaps u/v,
    O3 T2 -> T2a (exact path, critical) + T2b (general path, parallel with T3), O4 new task SK2 (stand-in
    update) after the amendment merges. Parser facts: skips unknown types silently; unknown mirror plane
    read as y=0 with a warning (T1 outer_rows raises UnsupportedMirror).
18. SK merged + verified 9.5 (8eb4552). General amendment run 2 ended at 7 (round 6): one break (column
    vertices on all-exact hulls). Orchestrator applied the 3 round-6 fixes (31fafbf on task/spec2);
    focused verification grader wt5s1q1fn running. T0b still in round-2 grading (workflow war84khjb).
19. T0b accepted 9 (c64a4d0, 4 rounds), merge wb060iyun running. TODO after merge (before lanes):
    restore AGENTS rule 10 + section 3 item 13 to state the rename explicitly (T0b reworded them to pass
    its grep); exemption from the old-name grep: tools/, tests/baseline/*.json, Output/, docs/plans.
    Contract errata (after spec2 merges, with SK2): F3 exempt from NotAnalytic; F9 stand-in by
    geo.analytic or hull_name, F10 by props.analytic; file-map row for tests/standin_kit/* (J1/J2);
    thin-shell z_ballast == inner z_lo cuts the inner set (no zero-thickness layer); S4 precast
    joint-over-inner precedence; tests/README Folders table -> T12.
    Lane models: T5/T6/T7 opus high; T8/T10/T0d/T4a/T4b sonnet high; T2a/T2b/T3 opus high.
    Shared-file edits whose row has an unmerged earlier task are deferred to J2 (report exact edit).
20. T0b merged + pushed b4523aa (merge grader wb060iyun still verifying). Lanes launched from b4523aa:
    U wnd5r4qmw (T5 -> T6, opus), S weybyc2pd (T7, opus), O w4bne3qif (T8 || T10, sonnet),
    P wyf5fwsix (T0d -> T4a -> T4b, sonnet). Spec2 verification wt5s1q1fn running. After the T0b merge
    grader finishes: restore AGENTS rule 10 / item 13 wording (commit on main). Then spec2 merge, contract
    errata, SK2, lane K (T2a -> T3, T2b).
21. T0b merge verified 10; AGENTS rename wording restored (1f9d19e, pushed). My spec2 column fix (31fafbf)
    rejected at 6: split seams not bitwise (parent parameter vs parser [0,1]). Third amendment run
    wqbo6pykj: O5 shared-edge equality = same degree, bitwise control points and weights, knots equal
    after affine map to [0,1] (rounding bound justified); F5 writes one edge curve per shared boundary;
    O6 corners only, split once at each corner height after constant-z splits, row ends set to corners;
    O7 split_wall_box: T2a asserts F1 decisions only, full-body checks move to J1 test_join_general.
22. Owner: keep a progress tracker and a handoff, updated every 100k tokens. HANDOFF.md at repo root
    (tracker table + handoff), orchestration scripts in docs/plans/orchestration/. Update at every
    milestone and when get_session context_usage grows by 100k (688k at 12:20 UTC). 7-day rate limit
    at warning level (reset 21:00 UTC).
23. spec2 accepted 9 (c9e2231, round 8). My errata commit 55671ef rejected 6 (equality cut vs S4 precast,
    pole case, SK2 must cover it, height bound unstated, T12/tests README row, stand-in identification,
    I2 empty weights). Lesson: do not hand-edit the dense contract; delegate to author + grader.
    Workflow wq4kyyvwb: Opus author fixes, xhigh grader (2 rounds), then merge task/spec2 + merge grader.
    HANDOFF.md must be updated after this merge (main busy until then).
24. Owner (13:35, 13:38): "you and any agent judges from metrics" - stated from the start (Stage-3 numeric
    checks 10-06 22:15, per-metric closest-fail report 22:50, fitting metrics 23:04). Now AGENTS rule 12,
    HANDOFF section 3, lane.js rule 15 + grader item 7 (1ec249b, pushed). New launches: workflows/scripts/
    mwec-lane-v2.js; resume pre-13:50 workflows only with mwec-lane-wf_d681e461-e66.js. Never send the
    owner pictures as evidence; real C1 XZ/YZ sections from J1 on with metrics printed beside them.
25. Owner (14:50): "When you get the first REAL cross sections for UHPC, show me the figure and pause your
    work." -> at J1: real-kernel C1 UHPC XZ/YZ sections + metrics; show; TaskStop all workflows; no launches
    or merges until the owner replies. Recorded in HANDOFF section 7 item 3.
