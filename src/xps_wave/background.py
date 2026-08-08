"""Background models for XPS spectra: linear and the standard iterative Shirley background.

Convention: x is binding energy in ascending order. In a real XPS core-level
scan, the low-BE side of a region sits "before" the photoemission peak (fewer
inelastically-scattered electrons) and the high-BE side sits "after" it (more
scattered electrons contribute a step up in the background). Both functions
below rely on that convention; Spectrum already guarantees ascending energy.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import minimum_filter1d, uniform_filter1d


def linear_background(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Straight line between the first and last (x, y) points of the region."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x[-1] == x[0]:
        return np.full_like(y, y[0])
    slope = (y[-1] - y[0]) / (x[-1] - x[0])
    return y[0] + slope * (x - x[0])


def shirley_background(x: np.ndarray, y: np.ndarray, tol: float = 1e-6, max_iter: int = 200) -> np.ndarray:
    """Classic iterative Shirley (1972) background, active-electron-count formulation.

    B(x_i) = y_high + (y_low - y_high) * S(i) / S(0)

    where S(i) is the area of (y - B) from x_i to the high-BE end of the region,
    y_low = y[0] (low-BE endpoint) and y_high = y[-1] (high-BE endpoint). B is
    refined until it stops changing by more than `tol` (relative to the peak
    height) or `max_iter` is reached.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    y_low, y_high = y[0], y[-1]
    peak_scale = max(np.max(y) - min(y_low, y_high), 1e-12)

    background = np.linspace(y_low, y_high, x.size)
    for _ in range(max_iter):
        signal = y - background
        signal = np.clip(signal, 0.0, None)
        cum_from_low = cumulative_trapezoid(signal, x, initial=0.0)
        total_area = cum_from_low[-1]
        if total_area <= 1e-12:
            break
        area_i_to_high = total_area - cum_from_low
        new_background = y_high + (y_low - y_high) * area_i_to_high / total_area
        if np.max(np.abs(new_background - background)) < tol * peak_scale:
            background = new_background
            break
        background = new_background
    return background


def rolling_minimum_background(x: np.ndarray, y: np.ndarray, window_ev: float = 20.0) -> np.ndarray:
    """A rolling-minimum ("opening") baseline, for wide survey scans.

    A single linear or Shirley background is not meaningful across an entire
    survey (hundreds of eV, many unrelated core-level steps); this instead
    tracks the lower envelope of the data with a `window_ev`-wide rolling
    minimum, smooths it, and clips it to never exceed the data. It only
    needs to be good enough to expose candidate peaks to `scipy.find_peaks`
    (see xps_wave.survey), not to be quantitatively accurate - each detected
    peak is re-fit in its own narrow window with a local linear/Shirley
    background afterwards.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    step = float(np.median(np.diff(x))) if x.size > 1 else 1.0
    window = max(int(round(window_ev / step)), 3)
    baseline = minimum_filter1d(y, size=window, mode="nearest")
    baseline = uniform_filter1d(baseline, size=max(window // 2, 1), mode="nearest")
    return np.minimum(baseline, y)


def subtract_background(x: np.ndarray, y: np.ndarray, kind: str = "shirley", **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Return (background, signal_above_background) for the requested background kind."""
    if kind == "shirley":
        background = shirley_background(x, y, **kwargs)
    elif kind == "linear":
        background = linear_background(x, y)
    elif kind == "rolling_min":
        background = rolling_minimum_background(x, y, **kwargs)
    else:
        raise ValueError(f"unknown background kind: {kind!r} (expected 'shirley', 'linear', or 'rolling_min')")
    signal = y - background
    return background, signal
