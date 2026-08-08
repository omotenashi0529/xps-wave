"""Peak deconvolution ("waveform separation") for XPS regions using lmfit.

Two entry points:

* ``fit_peaks`` - generic multi-peak deconvolution for an arbitrary region.
  Peak count is auto-detected from local maxima unless given explicitly.
* ``fit_ti_2p_doublets`` - physics-constrained fit for the Ti 2p region: each
  chemical state is modelled as a spin-orbit doublet (2p3/2 + 2p1/2) with the
  literature splitting and the fixed 2:1 (4:2 sub-shell degeneracy) area
  ratio built in, so the fit cannot converge on a spin-orbit-violating
  solution. This is the model used for the Ti-focused analysis and for the
  accuracy evaluation against the real reference spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from lmfit.model import ModelResult
from lmfit.models import PseudoVoigtModel
from scipy.signal import find_peaks

from xps_wave.reference import TI_2P_STATES, TiState


@dataclass
class PeakResult:
    label: str
    center: float  # eV, binding energy - the peak's X-axis value
    height: float  # Y-axis intensity (counts) of the separated waveform at `center`
    area: float
    fwhm: float
    fraction: float  # 0 = pure Gaussian, 1 = pure Lorentzian


def _build_composite(prefixes: list[str]) -> PseudoVoigtModel:
    models = [PseudoVoigtModel(prefix=p) for p in prefixes]
    composite = models[0]
    for m in models[1:]:
        composite = composite + m
    return composite


def _peak_result(result: ModelResult, prefix: str, label: str) -> PeakResult:
    p = result.params
    return PeakResult(
        label=label,
        center=p[f"{prefix}center"].value,
        height=p[f"{prefix}height"].value,
        area=p[f"{prefix}amplitude"].value,
        fwhm=p[f"{prefix}fwhm"].value,
        fraction=p[f"{prefix}fraction"].value,
    )


def fit_peaks(
    x: np.ndarray,
    y_above_background: np.ndarray,
    n_peaks: int | None = None,
    min_prominence: float | None = None,
) -> tuple[list[PeakResult], ModelResult]:
    """Fit a background-subtracted region as a sum of PseudoVoigt peaks."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y_above_background, dtype=float)
    step = float(np.median(np.diff(x))) if x.size > 1 else 1.0

    prominence = min_prominence if min_prominence is not None else 0.05 * max(y.max() - min(y.min(), 0.0), 1e-9)
    peak_idx, _ = find_peaks(y, prominence=max(prominence, 1e-9))

    if n_peaks is not None:
        if len(peak_idx) >= n_peaks:
            peak_idx = np.sort(peak_idx[np.argsort(y[peak_idx])[::-1][:n_peaks]])
        else:
            peak_idx = np.linspace(0, len(x) - 1, n_peaks + 2)[1:-1].astype(int)
    if len(peak_idx) == 0:
        peak_idx = np.array([int(np.argmax(y))])

    prefixes = [f"p{i}_" for i in range(len(peak_idx))]
    composite = _build_composite(prefixes)
    params = composite.make_params()
    for prefix, idx in zip(prefixes, peak_idx):
        params[f"{prefix}center"].set(value=x[idx], min=x.min(), max=x.max())
        params[f"{prefix}sigma"].set(value=max(5 * step, 0.3), min=step / 2, max=(x.max() - x.min()))
        params[f"{prefix}amplitude"].set(value=max(y[idx], 1e-6) * 2.0, min=0)
        params[f"{prefix}fraction"].set(value=0.5, min=0, max=1)

    result = composite.fit(y, params, x=x)
    peaks = [_peak_result(result, prefix, f"peak{i+1}") for i, prefix in enumerate(prefixes)]
    peaks.sort(key=lambda p: p.center)
    return peaks, result


def fit_ti_2p_doublets(
    x: np.ndarray,
    y_above_background: np.ndarray,
    states: list[TiState] | None = None,
    splitting_tolerance: float = 0.15,
) -> tuple[list[PeakResult], ModelResult]:
    """Fit the Ti 2p region as a sum of literature-constrained spin-orbit doublets.

    Each state in `states` (default: all 4 in TI_2P_STATES - metal, Ti(II),
    Ti(III), Ti(IV)) contributes two PseudoVoigt components (2p3/2, 2p1/2)
    whose separation is tied to the literature spin-orbit splitting (allowed
    to refine within +/- `splitting_tolerance` eV to absorb small calibration
    offsets) and whose amplitude ratio is fixed at the 2:1 sub-shell
    degeneracy ratio - not a free parameter, since that ratio is set by atomic
    physics, not by the sample.
    """
    states = states or TI_2P_STATES
    x = np.asarray(x, dtype=float)
    y = np.asarray(y_above_background, dtype=float)
    step = float(np.median(np.diff(x))) if x.size > 1 else 1.0

    prefixes = []
    for i in range(len(states)):
        prefixes += [f"s{i}a_", f"s{i}b_"]
    composite = _build_composite(prefixes)
    params = composite.make_params()

    y_span = max(y.max() - y.min(), 1e-6)
    for i, state in enumerate(states):
        a, b = f"s{i}a_", f"s{i}b_"

        # Sigma is bounded tightly around the literature FWHM for this state
        # rather than left free over the whole region: with several strongly
        # overlapping chemical states only ~1-2 eV apart, an unconstrained
        # width lets one component's amplitude collapse to zero while a
        # neighbour absorbs its area (a well-known degeneracy in multi-state
        # Ti 2p fitting) - the same reason real CasaXPS fits constrain
        # component widths from single-phase reference standards.
        sigma_a_guess = max(state.fwhm_2p32 / 2, step)
        sigma_b_guess = max(state.fwhm_2p12 / 2, step)

        # Center is bounded near the literature value for this state (not the
        # whole region): with 4 states only ~1.2-1.8 eV apart, a wide-open
        # center lets two components collapse onto the same real feature,
        # leaving another real feature completely unmodelled.
        params[f"{a}center"].set(value=state.be_2p32, min=state.be_2p32 - 0.8, max=state.be_2p32 + 0.8)
        params[f"{a}sigma"].set(value=sigma_a_guess, min=sigma_a_guess * 0.5, max=sigma_a_guess * 1.8)
        params[f"{a}amplitude"].set(value=y_span, min=0)
        params[f"{a}fraction"].set(value=0.5, min=0, max=1)

        params.add(f"split{i}", value=state.splitting,
                    min=state.splitting - splitting_tolerance,
                    max=state.splitting + splitting_tolerance)
        # The lower-j spin-orbit partner (2p1/2) sits at HIGHER binding energy
        # than 2p3/2 - e.g. Ti(IV) 2p3/2 = 458.6 eV, 2p1/2 = 464.3 eV.
        params[f"{b}center"].set(expr=f"{a}center + split{i}")
        params[f"{b}sigma"].set(value=sigma_b_guess, min=sigma_b_guess * 0.5, max=sigma_b_guess * 1.8)
        params[f"{b}amplitude"].set(expr=f"{a}amplitude / {state.area_ratio_32_to_12}")
        params[f"{b}fraction"].set(expr=f"{a}fraction")

    result = composite.fit(y, params, x=x)

    peaks = []
    for i, state in enumerate(states):
        peaks.append(_peak_result(result, f"s{i}a_", f"{state.state} 2p3/2"))
        peaks.append(_peak_result(result, f"s{i}b_", f"{state.state} 2p1/2"))
    peaks.sort(key=lambda p: p.center)
    return peaks, result
