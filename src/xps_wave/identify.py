"""Element / chemical-state identification: match fitted peak positions against
the reference binding-energy database (xps_wave.reference)."""

from __future__ import annotations

from dataclasses import dataclass, field

from xps_wave.peakfit import PeakResult
from xps_wave.reference import TI_2P_STATES, SURVEY_REFERENCE_TABLE, ReferenceLine, TiState


@dataclass
class Identification:
    peak_label: str
    peak_center: float
    matched_element: str | None
    matched_orbital: str | None
    matched_state: str | None
    reference_be: float | None
    delta_ev: float | None
    confidence: str  # "high" / "medium" / "low" / "unmatched"
    candidates: list[tuple[str, str, str, float, float]] = field(default_factory=list)
    # candidates: (element, orbital, state, reference_be, delta_ev), best-first


def _confidence(abs_delta: float, tolerance: float) -> str:
    if abs_delta <= 0.3:
        return "high"
    if abs_delta <= 0.7:
        return "medium"
    if abs_delta <= tolerance:
        return "low"
    return "unmatched"


def find_candidates(
    center: float,
    tolerance: float = 1.0,
    table: list[ReferenceLine] | None = None,
    top_n: int = 3,
) -> list[tuple[ReferenceLine, float]]:
    """Return up to `top_n` reference lines within `tolerance` eV of `center`, closest first."""
    table = table if table is not None else SURVEY_REFERENCE_TABLE
    scored = [(line, abs(line.be - center)) for line in table]
    scored = [s for s in scored if s[1] <= tolerance]
    scored.sort(key=lambda s: s[1])
    return scored[:top_n]


def identify_peaks(
    peaks: list[PeakResult],
    tolerance: float = 1.0,
    table: list[ReferenceLine] | None = None,
) -> list[Identification]:
    """Identify the most likely element/chemical-state for each fitted peak."""
    results = []
    for peak in peaks:
        matches = find_candidates(peak.center, tolerance=tolerance, table=table)
        if not matches:
            results.append(
                Identification(
                    peak_label=peak.label,
                    peak_center=peak.center,
                    matched_element=None,
                    matched_orbital=None,
                    matched_state=None,
                    reference_be=None,
                    delta_ev=None,
                    confidence="unmatched",
                )
            )
            continue
        best_line, best_delta = matches[0]
        results.append(
            Identification(
                peak_label=peak.label,
                peak_center=peak.center,
                matched_element=best_line.element,
                matched_orbital=best_line.orbital,
                matched_state=best_line.state,
                reference_be=best_line.be,
                delta_ev=peak.center - best_line.be,
                confidence=_confidence(best_delta, tolerance),
                candidates=[(l.element, l.orbital, l.state, l.be, peak.center - l.be) for l, _ in matches],
            )
        )
    return results


def identify_ti_doublet_peaks(
    peaks: list[PeakResult],
    states: list[TiState] | None = None,
) -> list[Identification]:
    """Identification for peaks produced by `peakfit.fit_ti_2p_doublets`.

    Those peaks already carry a structurally-certain (state, orbital)
    assignment from the constrained doublet model itself (encoded in
    `peak.label`, e.g. "TiO2 (rutile/anatase) 2p3/2") - this just looks up
    each one's literature binding energy for the delta-vs-reference column,
    instead of re-discovering the assignment via the generic nearest-neighbour
    search in `identify_peaks` (which only tabulates 2p3/2 lines and would
    have nothing sensible to match 2p1/2 peaks against).
    """
    states = states or TI_2P_STATES
    by_state_name = {s.state: s for s in states}
    results = []
    for peak in peaks:
        state_name, _, orbital = peak.label.rpartition(" ")
        state = by_state_name.get(state_name)
        if state is None:
            results.append(
                Identification(peak.label, peak.center, None, None, None, None, None, "unmatched")
            )
            continue
        reference_be = state.be_2p32 if orbital == "2p3/2" else state.be_2p32 + state.splitting
        delta = peak.center - reference_be
        results.append(
            Identification(
                peak_label=peak.label,
                peak_center=peak.center,
                matched_element="Ti",
                matched_orbital=orbital,
                matched_state=state.state,
                reference_be=reference_be,
                delta_ev=delta,
                confidence=_confidence(abs(delta), tolerance=1.5),
            )
        )
    return results
