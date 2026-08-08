"""End-to-end XPS analysis pipeline: load -> background-subtract -> peak-fit
(waveform separation) -> identify -> tabulate."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from lmfit.model import ModelResult

from xps_wave.background import subtract_background
from xps_wave.identify import Identification, identify_peaks, identify_ti_doublet_peaks
from xps_wave.peakfit import PeakResult, fit_peaks, fit_ti_2p_doublets
from xps_wave.reference import TiState
from xps_wave.spectrum import Spectrum


@dataclass
class AnalysisResult:
    spectrum: Spectrum
    background: np.ndarray
    signal: np.ndarray
    peaks: list[PeakResult]
    identifications: list[Identification]
    fit_result: ModelResult

    def table(self) -> pd.DataFrame:
        """Peak table: BE (X-axis), intensity (Y-axis), and the identified species."""
        rows = []
        for peak, ident in zip(self.peaks, self.identifications):
            rows.append(
                {
                    "label": peak.label,
                    "BE_eV": round(peak.center, 3),
                    "intensity": round(peak.height, 1),
                    "area": round(peak.area, 1),
                    "fwhm_eV": round(peak.fwhm, 3),
                    "element": ident.matched_element,
                    "orbital": ident.matched_orbital,
                    "state": ident.matched_state,
                    "ref_BE_eV": ident.reference_be,
                    "delta_eV": None if ident.delta_ev is None else round(ident.delta_ev, 3),
                    "confidence": ident.confidence,
                }
            )
        return pd.DataFrame(rows)


class XPSAnalyzer:
    """Convenience wrapper: background subtraction + peak deconvolution + identification."""

    def __init__(self, spectrum: Spectrum):
        self.spectrum = spectrum

    def analyze_region(
        self,
        low: float | None = None,
        high: float | None = None,
        background_kind: str = "shirley",
        n_peaks: int | None = None,
        identify_tolerance: float = 1.0,
    ) -> AnalysisResult:
        """Generic region fit: auto-detected peak count, unconstrained PseudoVoigt peaks."""
        region = self.spectrum.window(low, high) if low is not None and high is not None else self.spectrum
        background, signal = subtract_background(region.energy, region.intensity, kind=background_kind)
        peaks, fit_result = fit_peaks(region.energy, signal, n_peaks=n_peaks)
        identifications = identify_peaks(peaks, tolerance=identify_tolerance)
        return AnalysisResult(region, background, signal, peaks, identifications, fit_result)

    def analyze_ti_2p(
        self,
        low: float = 450.0,
        high: float = 470.0,
        background_kind: str = "shirley",
        states: list[TiState] | None = None,
        splitting_tolerance: float = 0.15,
    ) -> AnalysisResult:
        """Ti 2p specific analysis: literature spin-orbit-doublet-constrained fit.

        Each candidate chemical state (metal / Ti(II) / Ti(III) / Ti(IV)) is
        fit as a 2p3/2+2p1/2 pair with the known splitting and 2:1 area ratio
        built into the model (see peakfit.fit_ti_2p_doublets), which is what
        makes this far more robust than a generic N-peak fit for this element.
        """
        region = self.spectrum.window(low, high, name=f"{self.spectrum.name} Ti 2p")
        background, signal = subtract_background(region.energy, region.intensity, kind=background_kind)
        peaks, fit_result = fit_ti_2p_doublets(
            region.energy, signal, states=states, splitting_tolerance=splitting_tolerance
        )
        identifications = identify_ti_doublet_peaks(peaks, states=states)
        return AnalysisResult(region, background, signal, peaks, identifications, fit_result)
