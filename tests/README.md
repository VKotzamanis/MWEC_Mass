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
| `octave_shims/` | Replacements for MATLAB functions Octave 8.4 lacks (`fmincon` over `sqp`, `optimoptions`, `contains`, `discretize`, `datetime`, `issorted`, `java` (SHA-256 only), `double`/`logical` (`.empty` form), `startsWith`, `endsWith`). See its README. |
| `shim/` | Tests of the shims (`test_fmincon_shim.m`, `test_misc_shims.m`). |
| `geometry/` | Geometry tests against the independent reference (`test_c1_reference_sections.m`). |
| `reference/` | `c1_reference.py`: evaluator of `Input/C1.ms2` written from the deck text; run on request by the geometry test. |
| `regression/` | `test_pipeline_baseline.m`: reruns the Octave pipeline and requires the numbers in `baseline/octave_v1_baseline.json` (see Pipeline under Octave). |
| `baseline/` | `octave_v1_baseline.json` (this repository's Octave pipeline, v1.0 code) and `matlab_v1_reference.json` (the MATLAB v1.0 results read from `Output/C1_*_results.mat`). |
| `fixtures/` | `C1_wamit_cache_v5.mat`: the BEM cache converted to a format Octave can read. |

`tools/` holds the scripts that make the baseline and fixture files: `convert_v73_to_v5.py`,
`make_matlab_v1_reference.py`, `baseline_run.m`, `write_octave_baseline.m`.

## Pipeline under Octave

`baseline_run.m` runs `mwecmass.driver.run` on a temporary copy of `src/`, the deck and the converted
cache (the pipeline derives its output folder from the location of `src/`), with every figure, log
and diagnostic switch off except the results MAT-file, which is written into the temporary folder and
read back. The repository's `Output/` is never written and no production file is changed. Octave
cannot read the v7.3 BEM cache, hence the converted copy in `fixtures/`. The `fmincon` shim calls the
`OutputFcn` only at `init` and `done`, and Octave's `sqp` is not MATLAB's, so Octave and MATLAB
numbers differ (Stage 1 and Stage 2 end in different local optima); the baselines pin the Octave
numbers for regression and `matlab_v1_reference.json` holds the MATLAB v1.0 values for comparison.

Two presets, one baseline file each, one JSON member per line:

| Preset | Inputs | Baseline file |
|---|---|---|
| `full` | author inputs of `WEC_User_Input.m`, except `modular_precast.n_sub = 101` (the default 100 makes `build_config` read a roundoff-sized polygon area in the circle-based floor code) | `baseline/octave_v1_baseline.json` |
| `fast` | `full` plus `geometry.n_z_levels = 40`, `thin_shell.n_z_grid = 60`, `modular_precast.n_z_grid = 60` | `baseline/octave_v1_baseline_fast.json` |

`test_pipeline_baseline.m` reruns the `fast` preset by default and requires the recorded text
`jsonencode(summary)` to be identical (every double bit-for-bit). Environment variables:
`TESTS_BASELINE_PRESET=full` for the full case, `TESTS_BASELINE_MODES=thin_shell` (comma separated)
to limit the modes. Regenerate a baseline, one process per mode if wanted:

```bash
octave --no-gui --quiet --eval "addpath('tools'); write_octave_baseline({'thin_shell'}, 'fast')"
```
