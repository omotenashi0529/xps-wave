import numpy as np

from xps_wave.background import linear_background, shirley_background, subtract_background


def _synthetic_region(step_be: bool = False):
    x = np.linspace(450.0, 462.0, 200)
    step = 8.0 * (x - x[0]) / (x[-1] - x[0])  # rising step, low->high BE
    peak = 100.0 * np.exp(-0.5 * ((x - 456.0) / 0.8) ** 2)
    baseline = 10.0
    return x, baseline + step + peak


def test_linear_background_matches_endpoints():
    x, y = _synthetic_region()
    bg = linear_background(x, y)
    assert np.isclose(bg[0], y[0])
    assert np.isclose(bg[-1], y[-1])


def test_shirley_background_matches_endpoints_and_is_monotonic_step():
    x, y = _synthetic_region()
    bg = shirley_background(x, y)
    assert np.isclose(bg[0], y[0], atol=1e-3)
    assert np.isclose(bg[-1], y[-1], atol=1e-3)
    # Shirley background should rise (low-BE side lower, high-BE side higher)
    assert bg[-1] > bg[0]


def test_shirley_background_is_near_zero_signal_away_from_the_peak():
    x, y = _synthetic_region()
    _, signal = subtract_background(x, y, kind="shirley")
    far_from_peak = np.abs(x - 456.0) > 3.0
    assert np.max(np.abs(signal[far_from_peak])) < 2.0


def test_subtract_background_recovers_positive_signal_under_peak():
    x, y = _synthetic_region()
    background, signal = subtract_background(x, y, kind="shirley")
    peak_idx = np.argmin(np.abs(x - 456.0))
    assert signal[peak_idx] > 50.0  # most of the Gaussian peak height survives
    assert signal.min() > -5.0  # no large negative excursions
