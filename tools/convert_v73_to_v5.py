"""Convert a MATLAB v7.3 (HDF5) .mat file into a v5 .mat file that GNU Octave 8.4 can load.

Octave cannot read the cell arrays of a v7.3 file (object references). Test-only tool; the
production Input/ file is never modified.

Usage: python3 tools/convert_v73_to_v5.py <in_v73.mat> <out_v5.mat> [variable ...]
Without variable names every top-level variable is converted.
Requires h5py, numpy and scipy.
"""
import sys

import h5py
import numpy as np
from scipy.io import savemat


def to_matlab_order(a):
    # HDF5 stores MATLAB arrays with reversed dimensions.
    return a.T if a.ndim >= 2 else a


def read_dataset(f, ds):
    cls = ds.attrs.get("MATLAB_class", b"").decode()
    raw = ds[()]
    if cls == "cell":
        out = np.empty(raw.shape, dtype=object)
        for idx in np.ndindex(raw.shape):
            out[idx] = read_node(f, f[raw[idx]])
        return to_matlab_order(out)
    if cls == "char":
        arr = np.asarray(raw)
        if arr.size == 0:
            return ""
        return "".join(chr(c) for c in arr.T.ravel())
    if cls == "canonical empty":
        return np.zeros((0, 0))
    if raw.dtype.names and "real" in raw.dtype.names:
        raw = raw["real"] + 1j * raw["imag"]
    if cls == "logical":
        raw = raw.astype(bool)
    return to_matlab_order(np.asarray(raw))


def read_node(f, node):
    if isinstance(node, h5py.Group):
        return {k: read_node(f, v) for k, v in node.items() if k != "#refs#"}
    return read_dataset(f, node)


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    src, dst, names = argv[1], argv[2], argv[3:]
    with h5py.File(src, "r") as f:
        keys = names or [k for k in f.keys() if k != "#refs#"]
        data = {k: read_node(f, f[k]) for k in keys}
    savemat(dst, data, format="5", do_compression=False)


if __name__ == "__main__":
    main(sys.argv)
