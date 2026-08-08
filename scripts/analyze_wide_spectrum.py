"""Demo of the wide-spectrum (survey scan) analysis path.

This uses a SYNTHETIC survey spectrum, not a downloaded measurement - unlike
the Ti 2p accuracy evaluation (scripts/evaluate_accuracy.py), no open,
freely-redistributable *raw* wide-survey XPS spectrum with a clear license
was found while building this project (see README.md, "Data provenance").
The synthetic spectrum is still physically grounded: peak positions come
from xps_wave.reference.SURVEY_REFERENCE_TABLE (real literature binding
energies), each peak has a realistic Ti-2p-like FWHM, and a cumulative-step
background plus Poisson noise are added to mimic what a real Al-Kalpha
survey looks like. Its only job here is to exercise the peak-search +
per-peak fit + identify code path across a wide energy range; it is not a
substitute for the accuracy evaluation.

Run: uv run python scripts/analyze_wide_spectrum.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xps_wave.spectrum import Spectrum  # noqa: E402
from xps_wave.survey import analyze_wide_spectrum  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

# (element label, BE eV, relative height, FWHM eV) - a plausible surface
# composition for a titanium-oxide sample with typical adventitious
# contamination, spanning most of a standard Al Kalpha survey window.
SYNTHETIC_PEAKS = [
    ("O 1s", 530.0, 900.0, 1.6),
    ("Ti 2p3/2", 458.6, 700.0, 1.3),
    ("Ti 2p1/2", 464.3, 350.0, 1.8),
    ("C 1s", 284.8, 300.0, 1.5),
    ("Si 2p", 103.3, 120.0, 1.7),
    ("Na 1s", 1071.7, 80.0, 2.0),
]


def make_synthetic_survey(rng_seed: int = 0) -> Spectrum:
    x = np.arange(0.0, 1150.0, 0.5)
    rng = np.random.default_rng(rng_seed)

    background = 400.0 + 0.9 * x  # inelastic-scattering background rises with BE
    for _, be, height, _ in SYNTHETIC_PEAKS:
        background += 0.6 * height / (1.0 + np.exp((x - be) / 3.0))  # step after each peak

    signal = np.zeros_like(x)
    for _, be, height, fwhm in SYNTHETIC_PEAKS:
        sigma = fwhm / 2.355
        signal += height * np.exp(-0.5 * ((x - be) / sigma) ** 2)

    intensity = background + signal
    intensity = rng.poisson(np.clip(intensity, 0, None)).astype(float)
    return Spectrum(energy=x, intensity=intensity, name="synthetic wide survey (demo)")


def main() -> None:
    spectrum = make_synthetic_survey()
    # prominence_fraction/min_distance_ev are tuned higher than the library
    # defaults for this specific noise level and peak spacing - real usage
    # should tune these per dataset, same as any peak-picking algorithm.
    hits = analyze_wide_spectrum(
        spectrum, window_half_width_ev=3.5, prominence_fraction=0.1, min_distance_ev=4.0
    )

    print(f"{'BE (eV)':>10} {'intensity':>10} {'fwhm':>6}  element  orbital  state             delta(eV)  confidence")
    for peak, ident in hits:
        print(
            f"{peak.center:10.1f} {peak.height:10.1f} {peak.fwhm:6.2f}  "
            f"{ident.matched_element or '?':7} {ident.matched_orbital or '-':7}  "
            f"{ident.matched_state or '-':16} "
            f"{'' if ident.delta_ev is None else f'{ident.delta_ev:+.2f}':>10}  {ident.confidence}"
        )

    true_positions = ", ".join(f"{label}={be}eV" for label, be, _, _ in SYNTHETIC_PEAKS)
    print(f"\n(ground truth used to build this synthetic spectrum: {true_positions})")

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(spectrum.energy, spectrum.intensity, "k-", lw=0.6)
    for peak, ident in hits:
        ax.axvline(peak.center, color="C1", lw=0.8, ls="--")
        ax.annotate(
            ident.matched_element or "?",
            (peak.center, peak.height * 1.05 + np.interp(peak.center, spectrum.energy, spectrum.intensity)),
            fontsize=8,
            ha="center",
        )
    ax.set_xlabel("Binding Energy (eV)")
    ax.set_ylabel("Intensity (counts)")
    ax.invert_xaxis()
    ax.set_title("Wide-spectrum peak search + identification (synthetic demo spectrum)")
    fig.tight_layout()
    RESULTS_DIR.mkdir(exist_ok=True)
    fig.savefig(RESULTS_DIR / "wide_spectrum_demo.png", dpi=150)
    print(f"\nsaved: {RESULTS_DIR / 'wide_spectrum_demo.png'}")


if __name__ == "__main__":
    main()
