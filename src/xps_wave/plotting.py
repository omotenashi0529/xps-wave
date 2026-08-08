"""Plotting helpers for AnalysisResult: raw data, background, fit, separated peaks."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from xps_wave.pipeline import AnalysisResult


def plot_analysis(result: AnalysisResult, ax: Axes | None = None, title: str | None = None) -> Axes:
    x = result.spectrum.energy
    y = result.spectrum.intensity
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    ax.plot(x, y, "k.", ms=3, label="measured", zorder=3)
    ax.plot(x, result.background, "--", color="gray", lw=1, label="background", zorder=2)
    fitted_total = result.background + result.fit_result.best_fit
    ax.plot(x, fitted_total, "r-", lw=1.5, label="fit total", zorder=4)

    components = result.fit_result.eval_components(x=x)
    for comp in components.values():
        ax.fill_between(x, result.background, result.background + comp, alpha=0.25, zorder=1)

    for peak in result.peaks:
        base = float(np.interp(peak.center, x, result.background))
        ax.annotate(
            f"{peak.label}\n{peak.center:.1f} eV",
            (peak.center, base + peak.height),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=7,
        )

    ax.set_xlabel("Binding Energy (eV)")
    ax.set_ylabel("Intensity (counts)")
    ax.invert_xaxis()
    ax.legend(fontsize=8, loc="best")
    if title:
        ax.set_title(title)
    return ax
