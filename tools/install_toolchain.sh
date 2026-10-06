#!/usr/bin/env bash
# Idempotent test-toolchain install (Ubuntu 24.04, root). Commands below are the ones that worked.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
# GNU Octave 8.4.0 (no octave-cli package on noble; 'octave' is the CLI+GUI meta pkg),
# gnuplot-nox for headless print(), octave-optim 1.6.2 (provides fmincon/nonlin_min; pulls octave-statistics)
apt-get install -y octave gnuplot-nox octave-optim
# STEP verification: pip wheel bundles OpenCASCADE; no extra system libs were needed.
python3 -m pip install gmsh
# sanity checks
octave --no-gui --quiet --eval "disp(version())"
python3 -c "import gmsh; print('gmsh', gmsh.__version__)"
