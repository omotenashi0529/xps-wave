"""Wide-spectrum (survey scan) analysis: locate every significant peak across a
full-range XPS spectrum, fit each locally, and identify the most likely
element / chemical state for each one.

This is the coarser sibling of the Ti-2p-specific doublet fit in
``peakfit.fit_ti_2p_doublets``: it trades chemical-state precision for
breadth, so that a full survey scan (0-1200 eV or whatever the instrument
covers) can be scanned automatically instead of region-by-region by hand.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks

from xps_wave.background import subtract_background
from xps_wave.identify import Identification, identify_peaks
from xps_wave.peakfit import PeakResult, fit_peaks
from xps_wave.spectrum import Spectrum


def find_survey_peaks(
    spectrum: Spectrum,
    prominence_fraction: float = 0.05,
    min_distance_ev: float = 3.0,
    background_kind: str = "rolling_min",
    smooth_window_ev: float = 1.5,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Background-subtract the whole spectrum and return candidate peak indices.

    A single linear or Shirley background across an entire survey scan
    (hundreds of eV, many unrelated core-level steps) is not physically
    meaningful, so the default here is a rolling-minimum envelope
    (`background.rolling_minimum_background`) just to expose candidate peak
    locations; each candidate is then re-fit in its own narrow window with a
    proper local background by `analyze_wide_spectrum`.

    Peak *detection* runs on a lightly smoothed copy of the signal
    (`smooth_window_ev`) so that counting-statistics shot noise - narrow,
    single-point spikes - doesn't get picked up as spurious peaks; the
    returned `signal` itself is unsmoothed, since that's what gets fit.
    """
    x, y = spectrum.energy, spectrum.intensity
    _, signal = subtract_background(x, y, kind=background_kind)
    step = float(np.median(np.diff(x))) if x.size > 1 else 1.0
    distance = max(int(round(min_distance_ev / step)), 1)
    smooth_window = max(int(round(smooth_window_ev / step)), 1)
    smoothed = uniform_filter1d(signal, size=smooth_window)
    prominence = prominence_fraction * max(smoothed.max() - smoothed.min(), 1e-9)
    idx, _ = find_peaks(smoothed, prominence=max(prominence, 1e-9), distance=distance)
    return x, signal, list(idx)


def analyze_wide_spectrum(
    spectrum: Spectrum,
    window_half_width_ev: float = 4.0,
    prominence_fraction: float = 0.05,
    min_distance_ev: float = 3.0,
    background_kind: str = "rolling_min",
    identify_tolerance: float = 1.0,
) -> list[tuple[PeakResult, Identification]]:
    """Locate every significant peak in a wide/survey spectrum, fit each locally
    with a single PseudoVoigt peak, and identify the most likely element/state.

    Returns a list of (PeakResult, Identification) pairs sorted by binding
    energy - the peak table plus its most likely element assignment that the
    user asked for.
    """
    x, signal, idx = find_survey_peaks(spectrum, prominence_fraction, min_distance_ev, background_kind)
    results: list[tuple[PeakResult, Identification]] = []
    for i in idx:
        center_guess = x[i]
        lo, hi = center_guess - window_half_width_ev, center_guess + window_half_width_ev
        mask = (x >= lo) & (x <= hi)
        if mask.sum() < 5:
            continue
        # Fine local re-baseline: the whole-spectrum rolling-min pass above
        # only exposes candidate peaks, so each one gets a proper local
        # Shirley background (appropriate for a single, roughly-isolated
        # core-level peak) before the single-peak fit.
        _, local_signal = subtract_background(x[mask], signal[mask], kind="shirley")
        peaks, _ = fit_peaks(x[mask], local_signal, n_peaks=1)
        peak = peaks[0]
        peak.label = f"BE={peak.center:.1f}eV"
        ident = identify_peaks([peak], tolerance=identify_tolerance)[0]
        results.append((peak, ident))
    results.sort(key=lambda pair: pair[0].center)
    return results
