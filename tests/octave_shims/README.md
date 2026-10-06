# Octave test shims

Test-only replacements for MATLAB functions that GNU Octave 8.4 lacks or implements differently.
`tests/run_tests.m` puts this folder on the Octave path **before** `src/`. Never add it to the
production path; MATLAB runs must not see these files.

| Shim | Why |
|---|---|
| `startsWith.m`, `endsWith.m` | Octave 8 strips a single-space pattern to empty and errors (`MS2Parser.m:81`). |
| `contains.m` | Missing in Octave 8 (`MS2Parser`, `build_config`, HAMS route). Char/cellstr, `'IgnoreCase'`. |
| `discretize.m` | Missing (`properties_2d`). Numeric edge vector form: bin index, NaN outside. |
| `datetime.m` | Missing. Only `datetime('now', 'Format', fmt)`, returned as char (`optim/run.m`, `export_figure.m`). |
| `issorted.m` | Octave lacks the `'strictascend'` mode (`load_hydro_cache.m`); other forms go to the built-in. |
| `optimoptions.m` | Missing. Returns a struct of the name/value pairs; the solver name goes in `SolverName`. |
| `fmincon.m` | Missing (Octave's optim-package version rejects infeasible starts). Built on core `sqp`; the help text states the problem mapping, the `exitflag` mapping, and what is not available. |

The shims implement the call forms the pipeline uses, not the whole MATLAB interface.
