# Tests

GNU Octave 8.4 runs the tests; production code stays MATLAB R2022b compatible. Install the toolchain
once per container with `bash tools/install_toolchain.sh` (Octave, gnuplot, octave-optim, Python gmsh,
h5py, numpy and scipy).

## Run

```bash
octave --no-gui --quiet tests/run_tests.m                       # every tests/**/test_*.m
TESTS_FILTER=geometry octave --no-gui --quiet tests/run_tests.m  # only paths containing "geometry"
```

`run_tests.m` puts `tests/octave_shims` before `src/` on the path, runs each test file in a
`try/catch`, prints `PASS`/`FAIL` with the elapsed time, prints a summary and exits with status 1 if
any test failed. Test file names must be unique across folders.

The pipeline regression in `regression/` takes 24 to 94 minutes per mode, so the default run lists
it as `SKIP`. Run it with `MWEC_REGRESSION=1`, for example the modular-precast `fast` preset
(about 24 minutes):

```bash
MWEC_REGRESSION=1 TESTS_FILTER=regression octave --no-gui --quiet tests/run_tests.m
```

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
| `geometry/` | Geometry tests against the independent reference (`test_c1_reference_sections.m`, `test_c1_reference_profile.m`) and `test_ms2_parser_identity.m`, which replays every public `MS2Parser` evaluation on `Input/C1.ms2` and on `fixtures/ms2_all_entity_types.ms2` (all entity types, references to undefined entities, entities with an empty brace list, and each public evaluator named after an entity type called with the struct of `model.entities`) and requires bit-identical outputs and errors against `fixtures/ms2_parser_reference.mat`. `ms2_parser_cases.m` and `ms2_parser_record.m` hold the decks, grids and calls. |
| `reference/` | `c1_reference.py`: evaluator of `Input/C1.ms2` written from the deck text: section area and inertias per height, the profile with outward normals, and the normal-offset half-width per height and wall thickness. |
| `regression/` | `test_pipeline_baseline.m`: reruns the Octave pipeline and requires the numbers recorded in `baseline/` (see Pipeline under Octave). Runs only with `MWEC_REGRESSION=1`. |
| `baseline/` | `octave_v1_baseline.json` (full preset, both modes), `octave_v1_baseline_fast.json` (fast preset, modular precast) and `matlab_v1_reference.json` (the MATLAB v1.0 results read from `Output/C1_*_results.mat`); `test_baseline_vs_matlab.m` prints the recorded Octave values next to the MATLAB ones from these files, without a pipeline run. The MATLAB file also holds the Stage-2 design (`results.Final3D`) per mode. |
| `step/`, `step_check.py` | STEP writer tests, `test_step_export_stage3.m` (Stage-3 file set on the fixture bodies) and the gmsh/OpenCASCADE import check of STEP files. |
| `fixtures/` | `C1_wamit_cache_v5.mat`: the BEM cache converted to a format Octave can read. `ms2_parser_reference.mat`: the outputs of `MS2Parser` before it resolved each entity once (written by `tools/make_ms2_parser_reference.m`; not to be regenerated from a changed parser). `ms2_all_entity_types.ms2`: a deck with every entity type. |

`tools/` holds the scripts that make the baseline and fixture files: `convert_v73_to_v5.py`,
`make_matlab_v1_reference.py`, `baseline_run.m`, `write_octave_baseline.m`, `make_ms2_parser_reference.m`.

## Pipeline under Octave

`baseline_run.m` runs `mwecmass.driver.run` on a temporary copy of `src/`, the deck and the converted
cache (the pipeline derives its output folder from the location of `src/`), with every figure, log
and diagnostic switch off except the results MAT-file, which is written into the temporary folder and
read back. The repository's `Output/` is never written and no production file is changed. Octave
cannot read the v7.3 BEM cache, hence the converted copy in `fixtures/`. The `fmincon` shim calls the
`OutputFcn` only at `init` and `done`, and Octave's `sqp` is not MATLAB's, so Octave and MATLAB
numbers differ (Stage 1 and Stage 2 end in different local optima); the baselines pin the Octave
numbers for regression and `matlab_v1_reference.json` holds the MATLAB v1.0 values for comparison.

Known difference: for modular precast Octave's Stage 1 picks another draft node (`vertical_shift`
0.7786, MATLAB 1.075) and Stage 2 stops with exitflag -2 at 0.409 (MATLAB 0.983, exitflag 1); the final
numbers agree with MATLAB only because today's Stage 3 re-optimises from scratch. Octave results are not
compared with MATLAB v1.0 numbers (owner, 2026-10-07): Stage-2 and Stage-3 tests check the formulation
directly, and the owner's MATLAB run is the final check.

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
