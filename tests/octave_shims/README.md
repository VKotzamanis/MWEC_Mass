# Octave test shims

Test-only replacements for MATLAB functions that GNU Octave 8.4 lacks or implements differently.
Put this folder on the Octave path **before** `src/` when running tests. Never add it to the
production path; MATLAB runs must not see these files.

| Shim | Why |
|---|---|
| `startsWith.m`, `endsWith.m` | Octave 8 strips a single-space pattern to empty and errors (`MS2Parser.m:81`). |

Further shims (`contains`, `discretize`, `datetime`, `optimoptions`, an `fmincon` wrapper over
Octave's core `sqp` that accepts infeasible starts) are added by task T0 of
`docs/plans/2026-10-06-exact-geometry-stage3-step.md`.
