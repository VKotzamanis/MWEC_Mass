export const meta = {
  name: 'mwec-merge',
  description: 'Merge accepted MWEC_Mass task branches into claude/lucid-cray-7o9442, run the suite, push, and have an Opus xhigh grader verify the merge',
  phases: [
    { title: 'Merge', detail: 'merge, test, push' },
    { title: 'Verify', detail: 'Opus xhigh grader checks the merge' },
  ],
}
const REPO = '/home/user/MWEC_Mass'
const MAIN = 'claude/lucid-cray-7o9442'
const TRAILER = 'Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01K6BHrg7DUndGoJ3QxpTeRr'
const branches = args.branches
const SCHEMA = { type: 'object', properties: { pre_head: { type: 'string' }, merged_head: { type: 'string' }, pushed: { type: 'boolean' }, conflicts: { type: 'array', items: { type: 'string' } }, tests_summary: { type: 'string' }, report: { type: 'string' } }, required: ['pre_head', 'merged_head', 'pushed', 'report'] }
const GRADE = { type: 'object', properties: { accepted: { type: 'boolean' }, score: { type: 'number' }, violations: { type: 'array', items: { type: 'string' } }, summary: { type: 'string' } }, required: ['accepted', 'score', 'summary'] }
const m = await agent(`Merge step of the MWEC_Mass refactor. In ${REPO} (branch ${MAIN}; the working tree must be clean, check first), record the current HEAD as pre_head, then merge these accepted task branches in this order with "git merge --no-ff <branch>": ${branches.join(', ')}. If a conflict appears, resolve it so that both sides' intended content is kept, and describe each resolution. Every merge commit message ends with:
${TRAILER}
Then run the whole suite: octave --no-gui --quiet ${REPO}/tests/run_tests.m (never set MWEC_REGRESSION). If a test fails, do not push; report it with its output. If all pass, push with: git -C ${REPO} push -u origin ${MAIN} (on a network error retry after 2, 4, 8, 16 s; never force-push). Modify no file except to resolve conflicts. Report pre_head, the merged HEAD SHA, whether the push succeeded, conflicts and the suite summary.`, { label: 'merge', phase: 'Merge', model: 'sonnet', effort: 'high', schema: SCHEMA })
if (!m) return { merged: false, error: 'merge agent died' }
const g = await agent(`You are the GRADER of the MWEC_Mass refactor. Verify a merge in ${REPO} on branch ${MAIN}: the branches ${branches.join(', ')} were merged on top of ${m.pre_head}. Merge report: ${JSON.stringify(m)}. Check: (1) git diff ${m.pre_head}..HEAD equals the union of the branches' own changes (for each branch, git diff $(git merge-base ${m.pre_head} <branch>)..<branch> appears unchanged in the result), except documented conflict resolutions, which you judge; (2) no other file changed; (3) rerun the whole suite yourself (octave --no-gui --quiet tests/run_tests.m; never set MWEC_REGRESSION) and report its summary; (4) git ls-remote origin ${MAIN} equals HEAD if the report says pushed. accepted = true only if all four hold; score 0-10.`, { label: 'verify', phase: 'Verify', model: 'opus', effort: 'xhigh', schema: GRADE })
return { merged: m.pushed, merge: m, verify: g }
