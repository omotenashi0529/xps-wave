"""Accuracy evaluation: run this project's own Ti 2p peak-deconvolution pipeline
against a real measured spectrum and score it against an independent,
expert-fitted ground truth.

Ground truth source: the "Mixed Titanium Sample" VAMAS file (see
scripts/fetch_reference_data.py) embeds the CasaXPS component fit that
M.C. Biesinger (Surface Science Western) published as the reference fit for
this exact spectrum - 8 peaks (4 chemical states x 2p3/2/2p1/2 spin-orbit
doublet), each with a binding energy, area, and width, independently derived
by an XPS expert using commercial software. This script re-extracts that
embedded fit as ground truth and compares it to what `xps_wave` recovers on
its own, using only the raw counts and the literature constraints in
`xps_wave.reference`.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xps_wave.io import read_vamas  # noqa: E402
from xps_wave.pipeline import XPSAnalyzer  # noqa: E402
from xps_wave.plotting import plot_analysis  # noqa: E402

from fetch_reference_data import DEST, fetch  # noqa: E402

_CASA_COMP_RE = re.compile(
    r"CASA comp \(\*([^*]+)\*\).*?Area\s+([\d.]+).*?Position\s+([\d.]+)"
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def parse_ground_truth(block_comment: str, photon_energy_eV: float) -> pd.DataFrame:
    rows = []
    for m in _CASA_COMP_RE.finditer(block_comment):
        label, area, ke = m.group(1), float(m.group(2)), float(m.group(3))
        rows.append({"label": label, "ground_truth_BE_eV": photon_energy_eV - ke, "ground_truth_area": area})
    df = pd.DataFrame(rows).sort_values("ground_truth_BE_eV").reset_index(drop=True)
    if len(df) != 8:
        raise RuntimeError(f"expected 8 CASA components in the reference file, parsed {len(df)}")
    return df


def match_nearest(ground_truth: pd.DataFrame, fitted_centers: list[float]) -> pd.DataFrame:
    """Greedy nearest-energy pairing between ground truth and fitted peaks (no reuse)."""
    remaining = list(enumerate(fitted_centers))
    matches = []
    for _, gt_row in ground_truth.iterrows():
        gt_be = gt_row["ground_truth_BE_eV"]
        best_j, (best_idx, best_be) = min(enumerate(remaining), key=lambda kv: abs(kv[1][1] - gt_be))
        remaining.pop(best_j)
        matches.append(best_be)
    result = ground_truth.copy()
    result["fitted_BE_eV"] = matches
    result["error_eV"] = result["fitted_BE_eV"] - result["ground_truth_BE_eV"]
    return result


def main() -> None:
    fetch()
    spectrum = read_vamas(DEST)
    print(f"loaded '{spectrum.name}' - {len(spectrum)} points, "
          f"calibration shift applied: {spectrum.metadata['calibration_shift_eV']:+.3f} eV\n")

    ground_truth = parse_ground_truth(spectrum.metadata["block_comment"], spectrum.metadata["photon_energy_eV"])

    result = XPSAnalyzer(spectrum).analyze_ti_2p()
    fitted_centers = [p.center for p in result.peaks]

    comparison = match_nearest(ground_truth, fitted_centers)
    comparison["abs_error_eV"] = comparison["error_eV"].abs()

    pd.set_option("display.width", 120)
    print("=== Peak position accuracy vs. expert CasaXPS reference fit (same spectrum) ===")
    print(comparison[["label", "ground_truth_BE_eV", "fitted_BE_eV", "error_eV"]].round(3).to_string(index=False))

    mae = comparison["abs_error_eV"].mean()
    rmse = np.sqrt((comparison["error_eV"] ** 2).mean())
    max_err = comparison["abs_error_eV"].max()
    worst = comparison.loc[comparison["abs_error_eV"].idxmax(), "label"]
    print(f"\nMAE  = {mae:.3f} eV")
    print(f"RMSE = {rmse:.3f} eV")
    print(f"max  = {max_err:.3f} eV  (worst: {worst})")

    # The two dominant chemical states (large area in the reference fit) are
    # the ones that matter most for "which element/state is this": report
    # their accuracy separately from the two minor, heavily-overlapping
    # intermediate oxidation states (Ti(II)/Ti(III)), which are intrinsically
    # harder to resolve from 2p data alone - see peakfit.fit_ti_2p_doublets.
    major = comparison.nlargest(4, "ground_truth_area")
    minor = comparison.nsmallest(4, "ground_truth_area")
    print(f"\nDominant species (Ti metal, TiO2) MAE = {major['abs_error_eV'].mean():.3f} eV")
    print(f"Minor/overlapping species (TiO, Ti2O3) MAE = {minor['abs_error_eV'].mean():.3f} eV")

    print("\n=== Full peak table with element/state identification ===")
    print(result.table().to_string(index=False))

    RESULTS_DIR.mkdir(exist_ok=True)
    comparison.round(3).to_csv(RESULTS_DIR / "ti2p_accuracy_report.csv", index=False)
    ax = plot_analysis(result, title=f"Ti 2p deconvolution vs. expert reference - MAE={mae:.2f} eV")
    ax.figure.tight_layout()
    ax.figure.savefig(RESULTS_DIR / "ti2p_accuracy.png", dpi=150)
    print(f"\nsaved: {RESULTS_DIR / 'ti2p_accuracy_report.csv'}")
    print(f"saved: {RESULTS_DIR / 'ti2p_accuracy.png'}")


if __name__ == "__main__":
    main()
