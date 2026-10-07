export const meta = {
  name: 'mwec-lane',
  description: 'Run one MWEC_Mass Phase-B lane: each task implemented and graded by an Opus xhigh grader until accepted (chain or parallel); branches pushed, no merging',
  phases: [
    { title: 'Implement', detail: 'implementer per task, in its own worktree' },
    { title: 'Grade', detail: 'Opus xhigh grader, accepts only verified work scoring >= 9' },
    { title: 'Fix', detail: 'implementer addresses the grader fixes' },
  ],
}

const REPO = '/home/user/MWEC_Mass'
const MAIN = 'claude/lucid-cray-7o9442'
const SCRATCH = '/tmp/claude-0/-home-user-MWEC-Mass/4eed5a18-bc85-5f1a-b911-b09cefd00ef6/scratchpad'
const PLAN = 'docs/plans/2026-10-06-exact-geometry-stage3-step.md'
const CONTRACT = 'docs/plans/2026-10-07-interfaces.md'
const TRAILER = 'Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01K6BHrg7DUndGoJ3QxpTeRr'

const RULES = (t) => `You are an implementation subagent on the MWEC_Mass refactor (MATLAB code, tested with GNU Octave 8.4).
Worktree: ${t.wt} (branch ${t.branch}). If the worktree does not exist, create it first with: git -C ${REPO} worktree add ${t.wt} -b ${t.branch} ${t.base} . Work ONLY inside it; never edit ${REPO} or another worktree. Commit on your branch and push it after every commit with: git -C ${t.wt} push -q origin ${t.branch} (never push any other branch). Every commit message ends with these two lines:
${TRAILER}

Before writing code read, in your worktree: AGENTS.md (all), the contract ${CONTRACT} (all; it wins over the plan where they differ) and your task section of ${PLAN}. Hard rules:
1. The code is the source of truth: open and grep the .m files; never take behaviour from docs/.
2. Do not rename or move existing files. Create new files only where your task or the contract says.
3. Never assume a shape (AGENTS section 1 rules 1-3): geometry only from the MS2Parser parametric definition or splines fitted to points computed from it. A test may use known exact answers (contract C1 oracles, stand-in closed forms) as an independent oracle and must say so in one line.
4. Names (rule 10): z_ballast, rho_air, rho_ballast, rho_uhpc; no variable carries the name of one material and the value of another.
5. Comments only where they explain a non-obvious reason; no narration, history or task tags. Delete code and comments you make stale.
6. No invented tolerances (rule 5): tests assert exact-by-construction properties at a machine-precision bound you justify in one line, or gates named in AGENTS.md; approximation errors are printed, not gated.
7. Count every volume once (rule 11).
8. Owner instruction: no whole-pipeline run. Never set MWEC_REGRESSION. Test components: use the stand-ins (tests/standins, added with addpath(...,'-end')) for every producer that is not merged into your base yet, as the contract section 3 describes.
9. Tests live in tests/<area>/test_<name>.m (function files, no inputs or outputs, error() on failure, measured values printed). Run the suite with: octave --no-gui --quiet ${t.wt}/tests/run_tests.m
10. Production code stays MATLAB R2022b compatible (no Octave-only syntax or functions; shims cover Octave gaps in tests only).
11. Shared files: obey the contract's file ownership map (section 4). Edit a shared file only when every earlier task in its row is merged into your base; such edits are your last commit.
12. No side quests; do only your task.
13. Long runs: up to 10 minutes in the foreground; longer ones in the background, and wait for them before reporting. Leave no process running.
14. Commit after each completed deliverable.
15. Judge by metrics (AGENTS section 1 rule 12): every claim in your report is a number you computed from the code's outputs, with the command that produced it. Never judge a result, a figure included, by how it looks; a figure is judged by its section data.
16. Owner rule: one Octave process per core. /usr/bin/octave gives each run its own core (pinned, single-threaded BLAS) and waits while every core is busy. Run one Octave process at a time; never start parallel Octave runs.
Your final answer goes to the orchestrator: fill the schema factually. List every test command you ran with its key printed numbers, and every minor finding you defer, with the reason (nothing is dropped silently). Never claim anything you did not run or read.`

const IMPL_SCHEMA = { type: 'object', properties: {
  commit: { type: 'string' }, files_changed: { type: 'array', items: { type: 'string' } },
  tests: { type: 'array', items: { type: 'string' } }, acceptance: { type: 'array', items: { type: 'string' } },
  deferred: { type: 'array', items: { type: 'string' } }, open_issues: { type: 'array', items: { type: 'string' } },
  report: { type: 'string' } }, required: ['commit', 'report'] }

const GRADE_SCHEMA = { type: 'object', properties: {
  accepted: { type: 'boolean' }, score: { type: 'number' },
  scores: { type: 'object', properties: { correctness: { type: 'number' }, rules: { type: 'number' }, evidence: { type: 'number' }, quality: { type: 'number' }, scope: { type: 'number' } } },
  violations: { type: 'array', items: { type: 'string' } },
  required_fixes: { type: 'array', items: { type: 'object', properties: { file: { type: 'string' }, line: { type: 'string' }, issue: { type: 'string' }, fix: { type: 'string' }, verify: { type: 'string' }, severity: { type: 'string', enum: ['blocking', 'minor'] } }, required: ['issue', 'fix'] } },
  minor_only: { type: 'boolean' },
  spec_issues: { type: 'array', items: { type: 'string' } },
  own_test: { type: 'string' }, verified: { type: 'array', items: { type: 'string' } }, summary: { type: 'string' } },
  required: ['accepted', 'score', 'summary'] }

const gradePrompt = (t, impl, round) => `You are the GRADER of the MWEC_Mass refactor. Review one subagent's work independently and adversarially against AGENTS.md, the contract ${CONTRACT}, the task section of ${PLAN} and the brief. Do not modify tracked files; scratch files go under ${SCRATCH}/grader/${t.id}/.
Task ${t.id}, worktree ${t.wt}, branch ${t.branch}, base ${t.base}. Review round ${round}.
Brief given to the implementer:
<<<
${t.brief}
>>>
Implementer's report:
<<<
${JSON.stringify(impl)}
>>>
Procedure:
1. git -C ${t.wt} status must be clean, and origin/${t.branch} must equal HEAD. Inspect git log and git diff --name-status from the merge base with ${t.base}: renames (R) and changes outside the task's scope are violations; a shared file (contract section 4) edited out of turn is a violation.
2. Read every added or changed file completely. Check every rule: no shape assumption (grep the changed files for sqrt(.*pi), pi*r, compute_rmin, R_eq, cos_alpha, L_k, max_slope, homothetic), names, unneeded or stale comments, dead code, invented tolerances, volume counted once, MATLAB R2022b compatibility, scope. Check conformance with every contract item the task implements or consumes (signatures, struct fields, invariants, error identifiers).
3. Rerun the task's tests and the whole suite yourself (octave --no-gui --quiet tests/run_tests.m; never set MWEC_REGRESSION; never run the whole pipeline). Independently re-derive at least two key numbers or properties, and write at least one test of your own with a known exact answer in your scratch folder, run it, and report it in own_test.
4. For every acceptance item of the plan section, the contract and the brief: met with evidence you reproduced, or not met.
5. Minor findings: each must be fixed, or listed by the implementer as deferred with a reason; nothing dropped silently.
6. Report in spec_issues any requirement of the brief, plan or contract that contradicts AGENTS.md, the owner decisions or another contract item; do not enforce such a requirement.
7. Judge only by metrics (AGENTS section 1 rule 12): every violation, required fix and acceptance item cites a number you computed from the code's outputs (a test result, a measured value, an exact reference). Judge a figure by its section data (every drawn boundary has a face of the solid behind it, the drawn areas equal the body's section areas, the drawn ballast level equals z_ballast), never by looking at a rendering.
8. Mark every required fix with severity. blocking: a wrong or unverified result, a failing test, a rule violation (AGENTS section 1 or the owner decisions), a contract mismatch, missing evidence. minor: comment or docstring wording, report text (deferred-list entries, line numbers in them), a stray untracked or committed junk file, an unused variable or dead branch, a test message. Set minor_only = true when every required fix is minor and the work would be accepted once they are done; a separate quick check then verifies only those items, so make each minor fix's 'verify' a concrete command or read.
9. Score 0-10: correctness, rules, evidence, quality, scope; overall = the minimum. accepted = true ONLY if there is no violation, every acceptance item is met with evidence you reproduced, all tests pass, and the overall score >= 9. Otherwise list every required fix precisely (file, line, issue, fix, how to verify). Be strict; never accept what you could not verify.`

const fixPrompt = (t, g, round) => `${RULES(t)}

${t.brief}

The GRADER rejected round ${round} of your branch (score ${g.score}). Fix every item below in ${t.wt}, rerun the affected tests and the suite, commit, push, and report. Grader verdict:
<<<
${JSON.stringify({ violations: g.violations, required_fixes: g.required_fixes, spec_issues: g.spec_issues, summary: g.summary })}
>>>`

const QUICK_SCHEMA = { type: 'object', properties: { pass: { type: 'boolean' }, failed: { type: 'array', items: { type: 'string' } }, report: { type: 'string' } }, required: ['pass', 'report'] }

const quickPrompt = (t, g, impl) => `You are the QUICK CHECKER of the MWEC_Mass refactor (owner-approved: minor findings are confirmed by a quick check, not a full review). Task ${t.id}, worktree ${t.wt}, branch ${t.branch}. Do not modify tracked files.
Check exactly these items of the last full verdict, each by its 'verify' instruction, and nothing else:
<<<
${JSON.stringify(g.required_fixes || [])}
>>>
Also check: git -C ${t.wt} status is clean, origin/${t.branch} equals HEAD, and the commits since the verdict touch only what these items need (git diff from the commit the verdict reviewed: ${JSON.stringify(g.reviewed_commit || '')} or the previous HEAD named in the fixer report). Run only the test files those commits change or that cover the changed lines (one Octave process at a time). Fixer report:
<<<
${JSON.stringify(impl)}
>>>
pass = true only if every item is done as specified, the tree is clean and pushed, the changes stay within the items, and the tests you ran pass. List each failed item.`

async function minorPath(t, g, round) {
  // fixer corrects the minor items, a quick check confirms them (at most two tries)
  let impl = null
  for (let k = 1; k <= 2; k++) {
    const fixed = await agent(fixPrompt(t, g, round), { label: `fix-minor:${t.id}#${round}.${k}`, phase: 'Fix', model: t.model, effort: t.effort, schema: IMPL_SCHEMA })
    if (fixed) impl = fixed
    const q = await agent(quickPrompt(t, g, impl), { label: `quick:${t.id}#${round}.${k}`, phase: 'Grade', model: 'sonnet', effort: 'medium', schema: QUICK_SCHEMA })
    log(`${t.id} quick check ${round}.${k}: pass ${q && q.pass}`)
    if (q && q.pass) return { impl, accepted: true, quick: q }
    if (q) g = { ...g, required_fixes: (g.required_fixes || []).filter(f => (q.failed || []).some(x => x.includes(f.issue.slice(0, 40)) || x.includes(f.file || '#'))), summary: `Quick check failed: ${q.report}` }
  }
  return { impl, accepted: false }
}

async function runTask(t) {
  // t.resume = { last: grader verdict, r0: first grading round } continues a task from its last verdict
  let impl, r0 = 1, rounds = 4
  if (t.resume && t.resume.grade_first) {
    // continue with grading of committed work: t.resume.impl is the implementer's last report
    r0 = t.resume.r0
    rounds = t.resume.rounds || 2
    impl = t.resume.impl
  } else if (t.resume && t.resume.last && t.resume.last.minor_only) {
    const m = await minorPath(t, t.resume.last, t.resume.r0 - 1)
    if (m.accepted) return { id: t.id, branch: t.branch, accepted: true, score: t.resume.last.score, rounds: t.resume.r0 - 1, minor_fixes_quick_checked: true, commit: m.impl.commit, deferred: m.impl.deferred || [], grader: t.resume.last.summary, quick: m.quick.report }
    r0 = t.resume.r0
    rounds = 2
    impl = m.impl
  } else if (t.resume) {
    r0 = t.resume.r0
    rounds = 2
    impl = await agent(fixPrompt(t, t.resume.last, r0 - 1), { label: `fix:${t.id}#${r0 - 1}`, phase: 'Fix', model: t.model, effort: t.effort, schema: IMPL_SCHEMA })
  } else {
    impl = await agent(`${RULES(t)}\n\n${t.brief}`, { label: `impl:${t.id}`, phase: 'Implement', model: t.model, effort: t.effort, schema: IMPL_SCHEMA })
  }
  if (!impl) return { id: t.id, accepted: false, error: 'implementer died' }
  for (let round = r0; round < r0 + rounds; round++) {
    const g = await agent(gradePrompt(t, impl, round), { label: `grade:${t.id}#${round}`, phase: 'Grade', model: 'opus', effort: 'xhigh', schema: GRADE_SCHEMA })
    if (!g) return { id: t.id, accepted: false, error: 'grader died', rounds: round }
    log(`${t.id} round ${round}: score ${g.score}, accepted ${g.accepted}`)
    if (g.accepted && g.score >= 9) return { id: t.id, branch: t.branch, accepted: true, score: g.score, rounds: round, commit: impl.commit, deferred: impl.deferred || [], spec_issues: g.spec_issues || [], grader: g.summary }
    if (g.minor_only) {
      const m = await minorPath(t, g, round)
      if (m.impl) impl = m.impl
      if (m.accepted) return { id: t.id, branch: t.branch, accepted: true, score: g.score, rounds: round, minor_fixes_quick_checked: true, commit: impl.commit, deferred: impl.deferred || [], spec_issues: g.spec_issues || [], grader: g.summary, quick: m.quick.report }
      if (round === r0 + rounds - 1) return { id: t.id, branch: t.branch, accepted: false, score: g.score, rounds: round, violations: g.violations || [], required_fixes: g.required_fixes || [], spec_issues: g.spec_issues || [], grader: g.summary }
      continue
    }
    if (round === r0 + rounds - 1) return { id: t.id, branch: t.branch, accepted: false, score: g.score, rounds: round, violations: g.violations || [], required_fixes: g.required_fixes || [], spec_issues: g.spec_issues || [], grader: g.summary }
    const fixed = await agent(fixPrompt(t, g, round), { label: `fix:${t.id}#${round}`, phase: 'Fix', model: t.model, effort: t.effort, schema: IMPL_SCHEMA })
    if (fixed) impl = fixed
  }
}

const L = args
let results = []
if (L.mode === 'chain') {
  for (const t of L.tasks) {
    const r = await runTask(t)
    results.push(r || { id: t.id, accepted: false, error: 'task failed' })
    if (!r || !r.accepted) { log(`${t.id} not accepted: lane ${L.lane} stops here`); break }
  }
} else {
  results = (await parallel(L.tasks.map(t => () => runTask(t)))).map((r, i) => r || { id: L.tasks[i].id, accepted: false, error: 'task failed' })
}
return { lane: L.lane, results }
