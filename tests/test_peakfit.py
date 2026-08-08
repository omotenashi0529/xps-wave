import numpy as np

from xps_wave.peakfit import fit_peaks, fit_ti_2p_doublets
from xps_wave.reference import TI_2P_STATES


def test_fit_peaks_recovers_single_gaussian_center_and_area():
    x = np.linspace(280.0, 290.0, 300)
    true_center, true_sigma, true_amp = 284.8, 0.6, 500.0
    y = true_amp / (true_sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - true_center) / true_sigma) ** 2)

    peaks, result = fit_peaks(x, y, n_peaks=1)

    assert len(peaks) == 1
    assert abs(peaks[0].center - true_center) < 0.05
    assert peaks[0].height > 0
    assert result.success


def test_fit_peaks_separates_two_overlapping_gaussians():
    x = np.linspace(280.0, 290.0, 400)

    def gauss(c, s, a):
        return a * np.exp(-0.5 * ((x - c) / s) ** 2)

    y = gauss(284.8, 0.5, 100) + gauss(286.5, 0.5, 60)
    peaks, _ = fit_peaks(x, y, n_peaks=2)

    assert len(peaks) == 2
    centers = sorted(p.center for p in peaks)
    assert abs(centers[0] - 284.8) < 0.15
    assert abs(centers[1] - 286.5) < 0.15


def test_fit_ti_2p_doublets_respects_spin_orbit_area_ratio():
    x = np.linspace(450.0, 470.0, 400)
    state = TI_2P_STATES[3]  # TiO2, Ti(IV)

    def pseudo_gauss(c, fwhm, a):
        sigma = fwhm / 2.355
        return a * np.exp(-0.5 * ((x - c) / sigma) ** 2)

    y = (
        pseudo_gauss(state.be_2p32, state.fwhm_2p32, 1000.0)
        + pseudo_gauss(state.be_2p32 + state.splitting, state.fwhm_2p12, 500.0)
    )

    peaks, result = fit_ti_2p_doublets(x, y, states=[state])
    assert len(peaks) == 2
    p32, p12 = sorted(peaks, key=lambda p: p.center)
    assert abs(p32.center - state.be_2p32) < 0.1
    assert abs((p12.center - p32.center) - state.splitting) < 0.05
    # 2:1 degeneracy ratio must hold by construction (it is a fixed constraint,
    # not a free parameter of the fit)
    assert abs(p32.area / p12.area - state.area_ratio_32_to_12) < 1e-6
