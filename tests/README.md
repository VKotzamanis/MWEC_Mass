# Tests

GNU Octave 8.4 runs the tests; production code stays MATLAB R2022b compatible. Install the toolchain
once per container with `bash tools/install_toolchain.sh` (Octave, gnuplot, octave-optim, Python gmsh
and h5py).

## Run

```bash
octave --no-gui --quiet tests/run_tests.m                       # every tests/**/test_*.m
TESTS_FILTER=geometry octave --no-gui --quiet tests/run_tests.m  # only paths containing "geometry"
```

`run_tests.m` puts `tests/octave_shims` before `src/` on the path, runs each test file in a
`try/catch`, prints `PASS`/`FAIL` with the elapsed time, prints a summary and exits with status 1 if
any test failed. Test file names must be unique across folders.

## Conventions

- One file per test: `tests/<area>/test_<name>.m`, a function with no inputs and no outputs.
- A test calls `error()` to fail and prints its measured values with `fprintf`.
- Assert only what is exact by construction (symmetry, closure, identity of two code paths) at a
  machine-precision bound stated and justified in a comment, or a gate named in `AGENTS.md`.
  Approximation errors are measured and printed, never gated.
- C1-specific knowledge may serve as an independent oracle only; the test says so in one line.
- Shims are test-only and never go on the production path.

## Folders

| Folder | Holds |
|---|---|
| `octave_shims/` | Replacements for MATLAB functions Octave 8.4 lacks (`fmincon` over `sqp`, `optimoptions`, `contains`, `discretize`, `datetime`, `issorted`, `startsWith`, `endsWith`). See its README. |
| `shim/` | Tests of the shims (`test_fmincon_shim.m`). |
| `geometry/` | Geometry tests against the independent reference (`test_c1_reference_sections.m`). |
| `reference/` | `c1_reference.py`: evaluator of `Input/C1.ms2` written from the deck text; run on request by the geometry test. |
| `regression/` | `test_pipeline_baseline.m`: reruns the Octave pipeline and requires the numbers in `baseline/octave_v1_baseline.json`. |
| `baseline/` | `octave_v1_baseline.json` (this repository's Octave pipeline, v1.0 code) and `matlab_v1_reference.json` (the MATLAB v1.0 results read from `Output/C1_*_results.mat`). |
| `fixtures/` | `C1_wamit_cache_v5.mat`: the BEM cache converted to a format Octave can read. |

`tools/` holds the scripts that make the baseline and fixture files: `convert_v73_to_v5.py`,
`make_matlab_v1_reference.py`, `baseline_run.m`, `write_octave_baseline.m`.

## Pipeline under Octave

`baseline_run.m` runs `mwecmass.driver.run` on a temporary copy of `src/`, the deck and the converted
cache (the pipeline derives its output folder from the location of `src/`), with every figure, log
and diagnostic switch off, so the repository's `Output/` is never written. Octave's `fmincon` shim
calls the `OutputFcn` only at `init` and `done`, and `sqp` differs from MATLAB's `sqp`, so Octave and
MATLAB numbers differ; the baseline pins the Octave numbers for regression, and the MATLAB values are
printed next to them for reference.
