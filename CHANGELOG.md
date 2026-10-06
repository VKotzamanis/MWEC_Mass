# Changelog

##First public release (2026-09-16)

This is the first public release of the Gradient-Based Optimization Framework for Mass Distribution of Floating Wave Energy Converter Hulls. The release has been tested for compatibility in both MATLAB R2022b and R2026a.

### Included

- Runnable MATLAB source code for the complete mass-distribution workflow.
- Preliminary, thin-shell, and modular-precast realisation types.
- The cached C1 hydrodynamic workflow, including the supplied input deck and cache data.
- Corrected publication figures and checked C1 result artefacts in MATLAB, PNG, PDF, and log formats. MATLAB `.fig` files are included where generated.
- Public documentation for setup and execution: `README.md` and `docs/RUNTIME_GUIDE.md`; computational methods: `docs/METHODS_ENGINE.md`; result schema: `docs/RESULT_SCHEMA.md`; and the optional HAMS-MREL route: `docs/HAMS_MREL_ROUTE.md`.
- Verification across all three realisation types, with the corresponding C1 result records and expected output files.
- Opt-in diagnostics that operate on explicitly supplied result data.
- Two release archives: a source-only package without C1 input or output data, and a standalone C1 example package containing the source code, documentation, example inputs, and checked outputs.
- BSD 3-Clause licensing in `LICENSE` and citation metadata in `CITATION.cff`.

### Known caveats

- The external HAMS executable is not bundled. The HAMS-MREL route requires a separately installed solver and user configuration.
- The supplied hydrodynamic cache is tied to its C1 period grid and run conditions. Changes to those conditions require rebuilding or replacing the cache.
