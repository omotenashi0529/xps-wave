"""Generic CLI entry point: analyze a user-supplied XPS file end-to-end.

Background subtraction -> peak deconvolution -> element/chemical-state
identification -> peak table (CSV) + annotated plot (PNG). This is what the
`xps-analysis` skill (.claude/skills/xps-analysis/SKILL.md) calls when the
user hands over a new measurement file, as opposed to the two evaluation-
specific scripts (fetch_reference_data.py / evaluate_accuracy.py) which are
pinned to the bundled Ti 2p reference spectrum.

Examples:
    uv run python scripts/run_analysis.py data/raw/mixed_ti_2p.vms --mode ti2p
    uv run python scripts/run_analysis.py my_region.csv --mode region --low 280 --high 292 --n-peaks 2
    uv run python scripts/run_analysis.py my_survey.csv --mode wide --kinetic-energy --photon-energy 1486.6
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xps_wave.io import read_two_column, read_vamas  # noqa: E402
from xps_wave.pipeline import XPSAnalyzer  # noqa: E402
from xps_wave.plotting import plot_analysis  # noqa: E402
from xps_wave.survey import analyze_wide_spectrum  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("file", help="Path to a VAMAS (.vms) file or a 2-column text/CSV file")
    parser.add_argument("--format", choices=["auto", "vamas", "csv"], default="auto")
    parser.add_argument("--mode", choices=["ti2p", "region", "wide"], default="ti2p")
    parser.add_argument("--low", type=float, default=None, help="region/ti2p low BE (eV)")
    parser.add_argument("--high", type=float, default=None, help="region/ti2p high BE (eV)")
    parser.add_argument("--n-peaks", type=int, default=None, help="mode=region: fixed peak count (default: auto-detect)")
    parser.add_argument("--kinetic-energy", action="store_true", help="CSV energy column is kinetic energy, not binding energy")
    parser.add_argument("--photon-energy", type=float, default=1486.6, help="source photon energy (eV) for --kinetic-energy; default Al Kalpha")
    parser.add_argument("--delimiter", default=None, help="CSV delimiter (default: whitespace)")
    parser.add_argument("--skip-header", type=int, default=1)
    parser.add_argument("--background", choices=["shirley", "linear"], default="shirley")
    parser.add_argument("--out", default="results/custom_analysis", help="output directory for the peak table + plot")
    return parser.parse_args()


def load_spectrum(args: argparse.Namespace):
    path = Path(args.file)
    fmt = args.format
    if fmt == "auto":
        fmt = "vamas" if path.suffix.lower() == ".vms" else "csv"
    if fmt == "vamas":
        return read_vamas(path)
    return read_two_column(
        path,
        delimiter=args.delimiter,
        skip_header=args.skip_header,
        energy_is_kinetic=args.kinetic_energy,
        photon_energy=args.photon_energy if args.kinetic_energy else None,
    )


def main() -> None:
    args = parse_args()
    spectrum = load_spectrum(args)
    print(
        f"loaded '{spectrum.name}': {len(spectrum)} points, "
        f"BE range {spectrum.energy.min():.1f}-{spectrum.energy.max():.1f} eV"
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "wide":
        hits = analyze_wide_spectrum(spectrum)
        rows = [
            {
                "BE_eV": round(peak.center, 3),
                "intensity": round(peak.height, 1),
                "fwhm_eV": round(peak.fwhm, 3),
                "element": ident.matched_element,
                "orbital": ident.matched_orbital,
                "state": ident.matched_state,
                "confidence": ident.confidence,
            }
            for peak, ident in hits
        ]
        table = pd.DataFrame(rows)
        print(table.to_string(index=False))
        table.to_csv(out_dir / "wide_spectrum_peaks.csv", index=False)
        print(f"\nsaved: {out_dir / 'wide_spectrum_peaks.csv'}")
        return

    analyzer = XPSAnalyzer(spectrum)
    if args.mode == "ti2p":
        result = analyzer.analyze_ti_2p(low=args.low or 450.0, high=args.high or 470.0, background_kind=args.background)
    else:
        result = analyzer.analyze_region(low=args.low, high=args.high, n_peaks=args.n_peaks, background_kind=args.background)

    table = result.table()
    print(table.to_string(index=False))
    table.to_csv(out_dir / "peak_table.csv", index=False)

    ax = plot_analysis(result, title=Path(args.file).name)
    ax.figure.tight_layout()
    ax.figure.savefig(out_dir / "analysis.png", dpi=150)
    print(f"\nsaved: {out_dir / 'peak_table.csv'}")
    print(f"saved: {out_dir / 'analysis.png'}")


if __name__ == "__main__":
    main()
