#!/bin/bash
# Installed as /usr/bin/octave (the original moved to /usr/bin/octave.real).
# Owner rule: each Octave process runs on its own core, single-threaded, pinned there.
# A call takes a free core (one lock file per core, held until Octave exits), else waits for one.
N=$(nproc)
while true; do
  for ((k = 0; k < N; k++)); do
    exec 9>"/tmp/octave-core-$k.lock"
    if flock -n 9; then
      export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GOTO_NUM_THREADS=1
      exec taskset -c "$k" /usr/bin/octave.real "$@"
    fi
    exec 9>&-
  done
  sleep 2
done
