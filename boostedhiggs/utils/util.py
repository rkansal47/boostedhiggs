# coding: utf-8

"""
Some helper functions.
"""

import os
import subprocess
import numpy as np
import json
import law
import resource
from math import inf
from law.util import colored
import scipy


def ulimit(*, lower=False, **kwargs):
    for key, value in kwargs.items():
        res = getattr(resource, "RLIMIT_%s" % key.upper())
        soft, hard = [inf if v == resource.RLIM_INFINITY else v for v in resource.getrlimit(res)]
        if value in ("hard", max):
            value = hard
        else:
            value = min(value, hard)
        if not lower:
            value = max(value, soft)
        resource.setrlimit(
            res, tuple(resource.RLIM_INFINITY if v == inf else v for v in (value, hard))
        )


def calc_checksum(*paths, **kwargs):
    exclude = law.util.make_list(kwargs.get("exclude", ["*.pyc", "*.git*"]))
    exclude = " ".join("! -path '{}'".format(p) for p in exclude)

    sums = []
    for path in paths:
        path = os.path.expandvars(os.path.expanduser(path))
        if os.path.isfile(path):
            cmd = 'sha1sum "{}"'.format(path)
        elif os.path.isdir(path):
            cmd = (
                'files="$( find "{}" -type f {} -print | sort -z )"; '
                "(for f in $files; do sha1sum $f; done) | sha1sum".format(path, exclude)
            )
        else:
            raise IOError("file or directory '{}' does not exist".format(path))

        code, out, _ = law.util.interruptable_popen(
            cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        if code != 0:
            raise Exception("checksum calculation failed")

        sums.append(out.strip().split(" ")[0])

    if len(sums) == 1:
        return sums[0]
    else:
        cmd = 'echo "{}" | sha1sum'.format(",".join(sums))
        code, out, _ = law.util.interruptable_popen(
            cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        if code != 0:
            raise Exception("checksum combination failed")

        return out.strip().split(" ")[0]


class DotDict(dict):
    def __init__(self, *args, **kwargs):
        super(DotDict, self).__init__(*args, **kwargs)

    def __getattr__(self, attr):
        return self.get(attr)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self

    def __setstate__(self, state):
        self.update(state)
        self.__dict__ = self

    def __dir__(self):
        return sorted(dict.__dict__.keys() + self.keys())


def iter_chunks(iterable, chunksize):
    from itertools import zip_longest

    assert 0 < chunksize
    assert chunksize == int(chunksize)
    null = object()
    for chunk in zip_longest(*([iter(iterable)] * chunksize), fillvalue=null):
        yield [v for v in chunk if v is not null]

def parametrized(dec):
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)

        return repl

    return layer


def linsplit(n, l):
    import math

    lengths = []
    avg = float(n) / l
    rest = 0.0
    for i in range(l):
        length = avg + rest
        length_int = int(math.floor(length))
        rest = length - length_int
        lengths.append(length_int)
    return lengths


def round_significant(number, sig=1):
    from math import log10, floor

    return round(number, sig - int(floor(log10(abs(number)))) - 1)


def optimize_binning(
    full_edges,
    s_vals,
    b_vals,
    s_errs,
    b_errs,
    n_start_bins,
    n_min_bins,
    y_low,
    y_high,
    x_min=None,
    x_max=None,
    callback=None,
    silent=False,
):

    # defaults
    if s_errs is None:
        s_errs = np.zeros_like(s_vals)
    if b_errs is None:
        b_errs = np.zeros_like(b_vals)

    # some input checks
    assert s_vals.sum() > 0
    assert b_vals.sum() > 0
    assert len(s_vals) == len(full_edges) - 1
    assert len(b_vals) == len(full_edges) - 1
    assert s_errs.shape == s_vals.shape
    assert b_errs.shape == b_vals.shape
    assert n_start_bins >= n_min_bins
    assert y_low <= y_high
    _mini_bin_widths = full_edges[1:] - full_edges[:-1]
    assert _mini_bin_widths.max() - _mini_bin_widths.min() < 1e-6

    # helpers
    def raise_min_bins():
        if silent:
            return None
        else:
            raise Exception(
                "bin contents insufficient for n_min_bins {} and y_low {}: {}".format(
                    n_min_bins, y_low, s_vals.sum() + b_vals.sum()
                )
            )

    def select_vals(vals, start, stop):
        vals = np.array(vals)
        vals[start] = vals[: start + 1].sum()
        vals[stop - 1] = vals[stop - 1 :].sum()
        return vals[start:stop]

    def select_errs(errs, start, stop):
        errs = np.array(errs)
        errs[start] = (errs[: start + 1] ** 2.0).sum() ** 0.5
        errs[stop - 1] = (errs[stop - 1 :] ** 2.0).sum() ** 0.5
        return errs[start:stop]

    def sizes_to_edges(sizes):
        return full_edges[[0] + np.cumsum(sizes).tolist()]

    # when x_min or x_max are "auto", auto detect the centrally populated range
    vals = s_vals + b_vals
    if x_min == "auto":
        x_min = full_edges[np.argwhere(vals > 0).reshape(-1)[0]]
    if x_max == "auto":
        x_max = full_edges[np.argwhere(vals > 0).reshape(-1)[-1] + 1]

    # x_min and x_max define the approximate range of optimized edges to return, so when they are
    # set, find the outer most approximate edges and limit all arrays
    start, stop = 0, len(s_vals)
    if x_min is not None:
        start = int(np.argwhere(full_edges <= x_min).reshape(-1)[-1])
    if x_max is not None:
        stop = int(np.argwhere(full_edges >= x_max).reshape(-1)[0])
    full_edges = full_edges[start : stop + 1]
    s_vals, s_errs = select_vals(s_vals, start, stop), select_errs(s_errs, start, stop)
    b_vals, b_errs = select_vals(b_vals, start, stop), select_errs(b_errs, start, stop)

    # recompute things
    vals = s_vals + b_vals
    itg = vals.sum()
    # errs = (s_errs**2. + b_errs**2.)**0.5

    # detect early when the bin contents are insufficient to fill n_min_bins with at least y_low
    if itg < n_min_bins * y_low:
        return raise_min_bins()

    # start with the requested number of bins and an even binning
    # for easier handling, keep track of bin widths ("sizes" below) in units of bins defined by
    # full_edges ("mini bins"), e.g. if bin 0 has a size 5, it combines the first 5 mini bins
    n_bins = n_start_bins
    sizes = None
    while True:
        if sizes is None:
            sizes = linsplit(len(s_vals), n_bins)
            print("start from even binning with {} bins".format(colored(n_bins, "green")))
            if callable(callback):
                callback(sizes_to_edges(sizes))

        # get bin contents and errors
        split_points = np.cumsum(sizes)[:-1]
        binned_vals = np.array([sum(s) for s in np.split(vals, split_points)])
        # binned_errs = np.array([sum(s**2.)**0.5 for s in np.split(errs, split_points)])
        # binned_rels = binned_errs / binned_vals
        # binned_rels[np.isnan(binned_rels)] = 0.

        # identify bins that are below y_low / above y_high
        low_bins = np.argwhere(binned_vals < y_low).reshape(-1)
        high_bins = np.argwhere(binned_vals >= y_high).reshape(-1)

        # stop when there are no low bins
        if not len(low_bins):
            break

        # when there are no high bins with size > 1 to extract bin contents from
        # reduce the number of bins and start over
        high_bin_sizes = np.array([sizes[b] for b in high_bins])
        if not len(high_bins) or (high_bin_sizes == 1).all():
            n_bins -= 1
            if n_bins >= n_min_bins:
                print("reducing n_bins to {}".format(n_bins))
                sizes = None
                continue
            else:
                return raise_min_bins()

        # find the low bin with the smallest content, select the outermost in case of multiple ones
        smallest_low_bins = np.argwhere(binned_vals == binned_vals.min()).reshape(-1)
        dst_bin = sorted(smallest_low_bins, key=lambda i: abs(i - 0.5 * (n_bins - 1)))[-1]

        # find the widest high bin, select the one with the largest content in case of multiple ones
        widest_high_bins = high_bins[high_bin_sizes == high_bin_sizes.max()]
        src_bin = sorted(widest_high_bins, key=lambda i: binned_vals[i])[-1]

        # reduce the size of the widest high bin and increase that of the smallest low bin
        sizes[src_bin] -= 1
        sizes[dst_bin] += 1

    # convert sizes back into optimized edges
    edges = sizes_to_edges(sizes)

    # call the callback one last time
    if callable(callback):
        callback(edges)

    return edges


class NumpyEncoder(json.JSONEncoder):
    """ Custom encoder for numpy data types """

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):

            return int(obj)

        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)

        elif isinstance(obj, (np.complex_, np.complex64, np.complex128)):
            return {"real": obj.real, "imag": obj.imag}

        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()

        elif isinstance(obj, (np.bool_)):
            return bool(obj)

        elif isinstance(obj, (np.void)):
            return None

        return json.JSONEncoder.default(self, obj)


def wquant(x, w, q):
    """
    x: bin centers
    w: bin heights (bin content)
    q: quantiles
    """
    assert x.shape == w.shape
    assert x.ndim == 1
    assert q.ndim == 1
    assert np.all((0 <= q) & (q <= 1))
    i = np.argsort(x)
    x = x[i]
    w = w[i]
    c = np.cumsum(w)
    inter = scipy.interpolate.interp1d(c, x, kind="nearest")
    return inter(q[1:-1] * c[-1])
