from xps_wave.identify import find_candidates, identify_peaks, identify_ti_doublet_peaks
from xps_wave.peakfit import PeakResult
from xps_wave.reference import TI_2P_STATES


def _peak(label, center, height=100.0, area=100.0, fwhm=1.0, fraction=0.5):
    return PeakResult(label=label, center=center, height=height, area=area, fwhm=fwhm, fraction=fraction)


def test_find_candidates_matches_carbon_1s():
    matches = find_candidates(284.85, tolerance=0.5)
    assert matches
    best_line, delta = matches[0]
    assert best_line.element == "C"
    assert delta < 0.1


def test_identify_peaks_flags_unmatched_when_far_from_any_reference():
    result = identify_peaks([_peak("mystery", 999.0)], tolerance=0.5)
    assert result[0].confidence == "unmatched"
    assert result[0].matched_element is None


def test_identify_peaks_confidence_tiers():
    high = identify_peaks([_peak("p", 284.81)], tolerance=1.0)[0]
    medium = identify_peaks([_peak("p", 285.3)], tolerance=1.0)[0]
    assert high.confidence == "high"
    assert medium.confidence in ("medium", "low")


def test_identify_ti_doublet_peaks_uses_structural_labels():
    state = TI_2P_STATES[3]  # TiO2
    peaks = [
        _peak(f"{state.state} 2p3/2", state.be_2p32),
        _peak(f"{state.state} 2p1/2", state.be_2p32 + state.splitting),
    ]
    idents = identify_ti_doublet_peaks(peaks)
    assert all(i.matched_element == "Ti" for i in idents)
    assert {i.matched_orbital for i in idents} == {"2p3/2", "2p1/2"}
    assert all(abs(i.delta_ev) < 1e-6 for i in idents)
