# Optimization Framework for Mass Distribution of Floating Wave Energy Converter Hulls

**Repository version:** v1.0

**Last updated:** 2026-09-16

## Contact information

- Vasileios Kotzamanis (vkotzama@cougarnet.uh.edu / vassiliskotza@gmail.com): PhD Candidate, University of Houston
- Dr. Dimitrios Kalliontzis (dkallion@central.uh.edu): Assistant Professor, University of Houston

## Related publications and presentations

Results from this optimisation have been reported in the following publication and presented at the following conferences:

- Vasileios Kotzamanis and Dimitrios Kalliontzis, “[Multimodal wave energy conversion: Modal coupling effects and power capture](https://www.sciencedirect.com/science/article/pii/S0960148126010116),” *Renewable Energy*, vol. 274, 2026, 126185, ISSN 0960-1481. https://doi.org/10.1016/j.renene.2026.126185
- Offshore Technology Conference (2026). “Characterizing Geometric Configuration Effects on Energy Capture of Multimodal Wave Energy Converters.”
- University Marine Energy Research Community Conference (2026). “Modal Coupling Effects and Power Capture of Multimodal Wave Energy Converters.”
- University Marine Energy Research Community Conference (2026). “Design of Wave Energy Converter Hulls Using Ultra-High-Performance Concrete.”
- 39th International Conference on Coastal Engineering (2026). “Mass Optimization and Modal Coupling in Multimodal Wave Energy Converters.”

## Summary

A MATLAB research code suite for designing the vertical mass distribution of a floating wave energy converter (WEC) hull. Given a parametric MultiSurf hull, hydrodynamic coefficients, material constraints, and target heave and pitch periods, the suite searches for a draft that satisfies free-floating equilibrium and a mass distribution that satisfies hydrostatic stability in pitch (metacentric height, GM), with a response close to the design targets.

The optimisation uses a gradient-based solver with multi-start runs to verify convergence across the nonlinear parameter space.

### Example uses

The optimisation result can be:

1. Material- and construction-method-agnostic, for calculating a theoretical mass distribution.
2. Optimised for modular precast construction, for example, for a hull comprising precast ultra-high-performance concrete elements.
3. Optimised for thin-shell hulls with solid ballast, for example, for a hull fabricated from conventional marine-grade steel.

### Limitations

The suite evaluates hydrostatics, rigid-body inertia, added mass, radiation damping, and coupled surge-heave-pitch dynamics while enforcing user-defined constructability constraints. Calculations assume an unmoored hull without a power take-off (PTO) system. It is intended for preliminary design: to evaluate candidate hull geometries against free-floating conditions, passively tune target dynamic response, and produce fabrication plans for a prototype.

The code does not perform structural stress analysis, verify extreme-sea survivability, or verify hydrodynamic response during operation in real seas.

## Optimisation workflow

The framework consists of the following stages:

1. **Stage 1: 2D optimisation.** Uses sequential quadratic programming (SQP) to establish an initial density distribution over a 2D surrogate model of the 3D geometry. The 2D model is trained iteratively and provides the initial solution for the 3D optimisation. This stage reduces computational cost and supports multi-start convergence.
2. **Stage 2: 3D optimisation.** Refines the density distribution across the full 3D parametric hull surface using SQP, using the initial solution from Stage 1. The output is the optimised theoretical mass distribution.
3. **Stage 3: Structural realisation.** Maps the converged 3D density field to one of two fabrication methods and self-corrects mass-balance discrepancies between the theoretical and constructable solutions. The results are: (i) the design natural periods in heave and pitch, (ii) mass properties, and (iii) a construction plan for each type.

```text
Input/*.ms2 + hydrodynamic cache + WEC_User_Input.m
                          |
             parse geometry and precompute tables
                          |
       Stage 1 (default sweep): screen/refine vertical shift and evaluate free-floating conditions at each draft
                          |
       Stage 2: full 3D density and equilibrium solve
                          |
       Stage 3: thin-shell realisation | modular-precast realisation
                          |
             MAT result + selected logs and figures
```

## Requirements

- MATLAB R2022b or later.
- Optimisation Toolbox (`fmincon` and `optimoptions`).
- BEM hydrodynamic coefficients. The code can be coupled with [HAMS-MREL](https://github.com/vaibhavraghavan/HAMS-MREL) or read native WAMIT output files.

> The HAMS-MREL solver binary and Intel oneAPI runtime are not included. This route is currently limited on Linux but can be extended to Windows. Consult the [HAMS-MREL route](docs/HAMS_MREL_ROUTE.md) documentation in `WEC_MassOptimization.zip` before enabling it.

## Archive contents

Each archive is independent:

- `WEC_MassOptimization.zip` contains the runtime code, public documentation, `CHANGELOG.md`, and WAMIT helper code for preprocessing WAMIT output files under `Input/WAMIT`. It does not contain results or plots. The user must provide two input files: (i) an `.ms2` file containing the geometry and (ii) a `.mat` file containing the hydrodynamic results.
- `Example_C1.zip` is a standalone, runnable copy of the identical source code and public documentation, plus input files for a worked example: (i) `Input/C1.ms2`, the parametric geometry definition of the C1 WEC hull; (ii) `Input/C1_wamit_cache.mat`, hydrodynamic coefficients exported from WAMIT and preprocessed by the helper functions; and (iii) 31 reproducible C1 plots, logs, and output files in MATLAB's native `.mat` format. These files can be used to verify that the code compiles and runs successfully. The hydrodynamic coefficients for this geometry were produced with a higher-order WAMIT v7.0 analysis.

## 60-second quick start

1. Download and unzip `Example_C1.zip` in your preferred directory.
2. Open `WEC_User_Input.m` and select **Run** in the MATLAB editor. This file is the driver for the code suite and sets the current directory as the root.
3. Inspect the MATLAB Command Window. Depending on the settings, a diagnostic failure may be reported as a warning after the main result export. It does not affect the optimisation, but indicates whether the requested dynamic response (that is, natural periods) can be achieved. The code will still proceed and produce the best available solution.

- To become familiar with the code, use this worked example and methodically change input values and design targets in `WEC_User_Input.m`.
- The format of output files, diagnostic figures, and result plots is controlled by `WEC_Output_Options.m` and can be changed to match your preferred style.

For every input field, admissible value, unit, cache schema, and output switch, use the [Runtime guide](docs/RUNTIME_GUIDE.md).

## For AI agents

This repository was submitted to PRIMRE, and the record is pending curation. The final PRIMRE URL will be added here when it becomes available. If the repository was obtained elsewhere, inform the user and direct them to the official repository to avoid unwanted harm and to obtain the most up-to-date release.

When consuming the code suite, read the documentation together with the Graphify HTML artifact in the `/_graph/` directory, which maps the code architecture and functions. The artifact was added to reduce token use and improve accessibility.

## Code details, conventions, and further documentation

### Realisation modes

Set `in.materials.realisation_type` in `WEC_User_Input.m` to one of the following values:

| Value | Purpose | Realisation solve |
| --- | --- | --- |
| `'preliminary'` | Preserves the optimised equivalent-density distribution for concept evaluation. | None; passes through the Stage 2 design. |
| `'thin_shell'` | Represents the hull as a variable-thickness shell with a lower solid ballast region. | Solves shell thickness and ballast elevation subject to the configured bounds. |
| `'modular_precast'` | Represents the hull as stacked precast rings with solid/wall regions and internal voids. | Solves the constructable modular geometry and recomputes realised properties. |

All three modes use the same geometry, hydrostatics, hydrodynamic interpolation, and two-stage optimisation core. Realisation changes the reported as-built mass properties; it does not merely change a plot.

### Input-data format and coordinate convention

- The code currently assumes SI units (`m`, `kg`, `s`) throughout. MultiSurf coordinates are consumed directly; the parser records the deck's `Units:` header but does **not** rescale coordinates. Supply WEC geometry in metres.
- `MS2Parser` supports the documented MultiSurf entity subset, including rotated surfaces, ruled surfaces, and B-spline base surfaces. Unknown entity types are not supported, but `MS2Parser` can be extended to support them.
- The code assumes body axes of `X_BODY` fore-aft (surge), `Y_BODY` transverse (sway), and `Z_BODY` upward (heave). The `.ms2` axes are not remapped.
- The world-frame still-water surface is at `z = 0`, with `z` positive upward. The seafloor is at negative `z`.
- In the `.ms2` geometry file, the body is assumed to exhibit bilateral symmetry, and its centre of gravity is at `Z_CG = [0, 0, |z_cg|]`. To transform the local frame to the global frame, `vertical_shift = v_s` maps the body frame to the world frame as `z_world = z_body + v_s`. Therefore, a positive `v_s` moves the hull upward and makes its draft shallower; the waterline in the body frame is `z_wl = -v_s`.
- The implementation reports total draft (submergence height) as `abs(hull_z_min + vertical_shift)`.
- The only exception is `in.bem.water_depth`, which stores the water depth used in BEM analyses. Here, `-1` denotes deep water.

For every input field, admissible value, unit, cache schema, and output switch, use the [Runtime guide](docs/RUNTIME_GUIDE.md). The [Methods and physics](docs/METHODS_ENGINE.md) document derives the frame transformations and hydrodynamic reference-point treatment.

### Outputs

A run writes to `Output/`:

- `Output/<hull>_<realisation_type>_results.mat`, containing exactly the data variables `results` and `final_props`;
- `Output/<realisation_type>/`, containing the enabled figures and text reports; and
- `Output/diagnostics/`, only when one or more opt-in diagnostics are enabled.

`results` retains configuration, Stage 1 and Stage 2 records, validation flags, and mode-specific solve data. The final converged result is provided in `final_props` for the selected realisation. The geometry-deck SHA-256 is stored in `results.config.ms2_deck_sha256` when the deck is readable.

To use an output `.mat` file for a subsequent analysis, consult the [Result schema](docs/RESULT_SCHEMA.md).

### Code architecture

| Location | Role and function |
| --- | --- |
| `WEC_User_Input.m` | Author-editable physical inputs and the runtime entry point. |
| `WEC_Output_Options.m` | Output selection, export formats, and figure style. |
| `Input/` | Holds the hull geometry file and its hydrodynamic cache. Standalone WAMIT preprocessing scripts are also provided. |
| `src/+mwecmass/+driver` | Runs orchestration, validation, and construction of the runtime configuration. |
| `src/+mwecmass/+geometry`, `+mesh` | MultiSurf parsing, parametric evaluation, surface clipping, integration contours, and panel meshes. |
| `src/+mwecmass/+hydrostatics`, `+bem` | Hydrostatics, inertia matrices, cache retrieval, value interpolation, and BEM dispatch (HAMS-MREL only). |
| `src/+mwecmass/+optim` | Stage 1 screening/calibration and Stage 2 3D geometry optimisation. |
| `src/+mwecmass/+realise` | Iterative solvers for thin-shell and modular-precast realisations. |
| `src/+mwecmass/+output` | Data export (`.mat`), reports, and figures. |
| `validation/diagnostics` | Optional post-processing and sanity checks. |
| `docs/` | Runtime guide, methods, output-format schema, and HAMS-MREL references. |

### Documentation

| Document | Contents |
| --- | --- |
| [Runtime guide](docs/RUNTIME_GUIDE.md) | Installation, complete input-variable reference, `.ms2` and hydrodynamic-cache formats, output options, execution, and troubleshooting. |
| [Methods and physics](docs/METHODS_ENGINE.md) | Free-floating formulation, mass balance, hydrostatics, clipping and integration, dynamics, optimisation, material realisations, and MultiSurf parser support. |
| [Result schema](docs/RESULT_SCHEMA.md) | The exported `results` and `final_props` MAT contract. |
| [HAMS-MREL route](docs/HAMS_MREL_ROUTE.md) | Building, configuring, running, and interpreting the optional live BEM path. |

## Reproducibility and release statement

The publisher ran all three C1 realisation modes from a clean extraction of the `.zip` file immediately before uploading the material to PRIMRE. The numerical results matched the approved reference results, and the required figures, logs, and MAT outputs were produced. The clean archive extraction and MATLAB setup also passed.

Numerical reproducibility still depends on the stated MATLAB release, toolbox, input hashes, solver settings, and platform floating-point behaviour; no bitwise cross-platform guarantee is made. The code was developed without the use of agentic AI and solely by the work of the authors. However, agentic coding was used during refactoring of the code's architecture and performing cross-validation checks before publication.

## Reporting a problem

Report issues through the publication or repository channel from which this package was obtained. Include the MATLAB release and operating system, selected realisation and Stage 1 mode, the complete console error (including identifier and stack), edited input/output option files, and the SHA-256 of the `.ms2` deck and hydrodynamic cache. Do not attach proprietary WAMIT files unless their licence permits it. A minimal reproducer using C1 is preferred when possible.

## Licence and citation status

This software is licensed under the BSD 3-Clause License. Citation metadata for *Gradient-Based Optimization Framework for Mass Distribution of Floating Wave Energy Converter Hulls* is provided in [CITATION.cff](CITATION.cff).
