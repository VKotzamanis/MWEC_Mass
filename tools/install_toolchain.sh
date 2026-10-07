#!/usr/bin/env bash
# Idempotent test-toolchain install (Ubuntu 24.04, root). Commands below are the ones that worked.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
# GNU Octave 8.4.0 (no octave-cli package on noble; 'octave' is the CLI+GUI meta pkg),
# gnuplot-nox for headless print(), octave-optim 1.6.2 (provides fmincon/nonlin_min; pulls octave-statistics)
apt-get install -y octave gnuplot-nox octave-optim
# gmsh (STEP verification): the wheel bundles OpenCASCADE. h5py: reads MATLAB v7.3 files for tools/convert_v73_to_v5.py.
# numpy, scipy: tests/reference/c1_reference.py, tools/make_matlab_v1_reference.py and tools/convert_v73_to_v5.py.
python3 -m pip install gmsh h5py numpy scipy
# sanity checks
octave --no-gui --quiet --eval "disp(version())"
python3 -c "import gmsh; print('gmsh', gmsh.__version__)"
python3 -c "import numpy, scipy, h5py; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'h5py', h5py.__version__)"
