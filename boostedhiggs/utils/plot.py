import hist
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mplhep as hep
from coffea.hist.plot import poisson_interval, clopper_pearson_interval, normal_interval
import numpy as np
from functools import reduce
from operator import add


def plot1d(
    h,
    ax=None,
    clear=True,
    overlay=None,
    stack=False,
    overflow=False,
    line_opts=None,
    fill_opts=None,
    error_opts=None,
    legend_opts={},
    density=False,
    binwnorm=None,
    order=None,
    legend_labels={},
    legend_title=None,
    log=False,
    yrange=[],  # [min, max] -> [1e-3, 1e6] (only min: [1e-3, None])
):
    if ax is None:
        ax = plt.gca()
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()

    if line_opts is None and fill_opts is None and error_opts is None:
        if stack:
            fill_opts = {}
        else:
            line_opts = {}
            error_opts = {}

    if overlay is not None:
        assert len(h.axes) == 2
        overlay = h.axes[overlay]
        assert isinstance(overlay, hist.axis.StrCategory)
        (axis,) = tuple(
            ax for ax in h.axes if isinstance(ax, (hist.axis.Regular, hist.axis.Variable))
        )

    else:
        assert len(h.axes) == 1
        (axis,) = h.axes

    if not isinstance(axis, (hist.axis.Regular, hist.axis.Variable)):
        raise NotImplementedError("Plot a sparse axis (e.g. bar chart)")

    else:
        ax.set_xlabel(axis.label)
        if binwnorm:
            ax.set_ylabel("Entries / Bin width")
        elif h.metadata:
            ax.set_ylabel(h.metadata.get("y_title", "Entries"))
        else:
            ax.set_ylabel("Entries")

        edges = axis.edges
        if order is None:
            identifiers = [*overlay] if overlay is not None else [None]
        else:
            identifiers = order

        plot_info = {
            "identifier": identifiers,
            "label": [legend_labels[i] for i in identifiers] if legend_labels else identifiers,
            "sumw": [],
            "sumw2": [],
        }
        for i, identifier in enumerate(identifiers):
            if identifier is None:
                sumw = h.view(flow=overflow)["value"]
                sumw2 = h.view(flow=overflow)["variance"]
            else:
                sumw = h[identifier, ...].view(flow=overflow)["value"]
                sumw2 = h[identifier, ...].view(flow=overflow)["variance"]

            plot_info["sumw"].append(sumw)
            plot_info["sumw2"].append(sumw2)

        def w2err(sumw, sumw2):
            err = []
            for a, b in zip(sumw, sumw2):
                err.append(np.abs(poisson_interval(a, b) - a))
            return err

        kwargs = None
        if line_opts is not None and error_opts is None:
            _error = None
        else:
            _error = w2err(plot_info["sumw"], plot_info["sumw2"])
        if fill_opts is not None:
            histtype = "fill"
            kwargs = fill_opts
        elif error_opts is not None and line_opts is None:
            histtype = "errorbar"
            kwargs = error_opts
        else:
            histtype = "step"
            kwargs = line_opts
        if kwargs is None:
            kwargs = {}

        hep.histplot(
            plot_info["sumw"],
            edges,
            label=plot_info["label"],
            yerr=_error,
            histtype=histtype,
            ax=ax,
            density=density,
            binwnorm=binwnorm,
            stack=stack,
            **kwargs,
        )
        if stack and error_opts is not None:
            stack_sumw = np.sum(plot_info["sumw"], axis=0)
            stack_sumw2 = np.sum(plot_info["sumw2"], axis=0)
            err = poisson_interval(stack_sumw, stack_sumw2)
            if binwnorm is not None:
                err *= binwnorm / np.diff(edges)[None, :]
            opts = {
                "step": "post",
                "label": "Sum unc.",
                "hatch": "///",
                "facecolor": "none",
                "edgecolor": (0, 0, 0, 0.5),
                "linewidth": 0,
            }
            opts.update(error_opts)
            ax.fill_between(
                x=edges, y1=np.r_[err[0, :], err[0, -1]], y2=np.r_[err[1, :], err[1, -1]], **opts
            )

        if legend_opts is not None:
            _label = overlay.label if overlay is not None else ""
            ax.legend(title=_label, **legend_opts)
        else:
            ax.legend(title=_label)
        if log:
            ax.set_yscale("log")
        ax.autoscale(axis="x", tight=True)
        ax.autoscale(axis="y", tight=False)
        if not log:
            ax.set_ylim(0, None)
        if yrange:
            ax.set_ylim(*yrange)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        title=legend_title,
        loc="upper left",
        bbox_to_anchor=(1.04, 1),
        borderaxespad=0,
    )

    return ax


def plotratio(
    num,
    denom,
    ax=None,
    clear=True,
    overflow=False,
    error_opts=None,
    denom_fill_opts=None,
    guide_opts=None,
    unc="clopper-pearson",
    label=None,
    ratio_yticks=[],  # [start, stop, step] -> [0.5, 1.5, 0.1]
):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()

    if error_opts is None and denom_fill_opts is None and guide_opts is None:
        error_opts = {}
        denom_fill_opts = {}

    (naxis,) = num.axes
    (daxis,) = denom.axes
    assert isinstance(naxis, (hist.axis.Regular, hist.axis.Variable))
    assert isinstance(daxis, (hist.axis.Regular, hist.axis.Variable))
    assert all(naxis.edges == daxis.edges)

    ax.set_xlabel(naxis.label)
    ax.set_ylabel("Ratio")
    edges = naxis.edges
    centers = naxis.centers

    sumw_num = num.view(flow=overflow)["value"]
    sumw2_num = num.view(flow=overflow)["variance"]
    sumw_denom = denom.view(flow=overflow)["value"]
    sumw2_denom = denom.view(flow=overflow)["variance"]

    rsumw = sumw_num / sumw_denom
    if unc == "clopper-pearson":
        rsumw_err = np.abs(clopper_pearson_interval(sumw_num, sumw_denom) - rsumw)
    elif unc == "poisson-ratio":
        # poisson ratio n/m is equivalent to binomial n/(n+m)
        rsumw_err = np.abs(clopper_pearson_interval(sumw_num, sumw_num + sumw_denom) - rsumw)
    elif unc == "num":
        rsumw_err = np.abs(poisson_interval(rsumw, sumw2_num / sumw_denom ** 2) - rsumw)
    elif unc == "normal":
        rsumw_err = np.abs(normal_interval(sumw_num, sumw_denom, sumw2_num, sumw2_denom))
    else:
        raise ValueError("Unrecognized uncertainty option: %r" % unc)

    if error_opts is not None:
        opts = {"label": label, "linestyle": "none"}
        opts.update(error_opts)
        emarker = opts.pop("emarker", "")
        errbar = ax.errorbar(x=centers, y=rsumw, yerr=rsumw_err, **opts)
        plt.setp(errbar[1], "marker", emarker)
    if denom_fill_opts is not None:
        unity = np.ones_like(sumw_denom)
        denom_unc = poisson_interval(unity, sumw2_denom / sumw_denom ** 2)
        opts = {"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
        opts.update(denom_fill_opts)
        ax.fill_between(
            edges,
            np.r_[denom_unc[0], denom_unc[0, -1]],
            np.r_[denom_unc[1], denom_unc[1, -1]],
            **opts,
        )
    if guide_opts is not None:
        opts = {"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1}
        opts.update(guide_opts)
        ax.axhline(1.0, **opts)

    if ratio_yticks:
        start, stop, step = ratio_yticks
        ax.set_yticks(np.arange(start, stop + step, step).tolist())
        for label in ax.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        ax.set_ylim(start, stop)

    if clear:
        ax.autoscale(axis="x", tight=True)
        ax.set_ylim(0, None)

    return ax


def blind(h_sig, h_bg, thresh=1e-3, h_sig_std=+1, h_bg_std=-1, flow=True):
    bg_sum = h_bg[::sum, ...]
    sig_sum = h_sig[::sum, ...]
    idx = (
        bg_sum.view(flow=flow)["value"] + h_bg_std * np.sqrt(bg_sum.view(flow=flow)["variance"])
    ) <= thresh * (
        sig_sum.view(flow=flow)["value"] + h_sig_std * sig_sum.view(flow=flow)["variance"]
    )
    return idx
