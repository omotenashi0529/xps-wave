---
name: xps-analysis
description: "Use this skill whenever the user provides or references XPS (X-ray Photoelectron Spectroscopy) measurement data and wants peak deconvolution/waveform separation, binding-energy (X-axis) and intensity (Y-axis) extraction, or element/chemical-state identification. Trigger on: XPS, 光電子分光, X線光電子分光, 波形分離, ピーク分離, ピークフィッティング, デコンボリューション, 結合エネルギー, binding energy, VAMAS, .vms files, survey/wide spectrum analysis, Ti 2p or other core-level regions, Shirley background. Also use it when the user asks to extend this project to a new element/orbital, or to re-run the accuracy evaluation against the reference spectrum."
license: Project-internal skill for xps_wave — no external license constraints.
---

# XPS spectrum analysis (xps_wave)

Wraps the `xps_wave` Python package in this repo: background subtraction ->
peak deconvolution -> element/chemical-state identification, for both a
single core-level region (e.g. Ti 2p) and a full wide/survey scan. Built and
validated 2026-08-08/09 against a real measured Ti 2p spectrum (see
`README.md` "データ出所・精度評価" — 0.05 eV MAE on the dominant chemical
states).

## Setup (once per machine)

```bash
uv sync   # installs numpy, scipy, lmfit, pandas, matplotlib, vamas, requests
```

Everything below assumes the current directory is the repo root (where
`pyproject.toml` lives) and commands run via `uv run`.

## Analyzing a user-supplied file

`scripts/run_analysis.py` is the generic entry point — point it at whatever
the user handed over:

| Situation | Command |
|---|---|
| VAMAS (`.vms`) file, Ti 2p region | `uv run python scripts/run_analysis.py <file>.vms --mode ti2p` |
| VAMAS file, some other single region | `uv run python scripts/run_analysis.py <file>.vms --mode region --low <eV> --high <eV> [--n-peaks N]` |
| 2-column text/CSV, binding-energy already | `uv run python scripts/run_analysis.py <file>.csv --mode region --delimiter , --low <eV> --high <eV>` |
| 2-column text/CSV in **kinetic** energy | add `--kinetic-energy --photon-energy 1486.6` (Al Kalpha; use 1253.6 for Mg Kalpha) |
| Full survey/wide scan, unknown composition | `uv run python scripts/run_analysis.py <file> --mode wide` |

Every mode prints a peak table (BE eV, intensity, FWHM, identified
element/orbital/state, confidence) and writes it + an annotated plot to
`--out` (default `results/custom_analysis/`). Read the table back to the
user rather than re-deriving it — it's already in the right units and has
the identification attached.

**If `--mode ti2p` gives a bad fit** (peaks collapsing onto each other,
`fit_result.success == False` with a real error rather than the benign
"Could not estimate error-bars" message): the sample's actual composition
may not have all 4 chemical states, or the region window (`--low`/`--high`,
default 450-470 eV) needs adjusting to fully bracket both the 2p3/2 and
2p1/2 doublets. Pass `states=[...]` (a subset of `xps_wave.reference.TI_2P_STATES`)
programmatically via `XPSAnalyzer.analyze_ti_2p(states=...)` if only some
species are expected to be present — this removes the extra, competing
components entirely rather than just down-weighting them.

## Re-running the accuracy evaluation

```bash
uv run python scripts/evaluate_accuracy.py   # downloads the reference file if missing, needs network
```

Compares this project's own Ti 2p fit against the expert (CasaXPS)
component fit embedded in the reference VAMAS file's text comments. Re-run
this after changing anything in `background.py` or `peakfit.py` to confirm
accuracy hasn't regressed — the current baseline is 0.05 eV MAE on Ti
metal/TiO2, 0.55 eV MAE on the overlapping TiO/Ti2O3 minor states (see
`results/ti2p_accuracy_report.csv` for the exact numbers).

## Extending to a new element

1. **New chemical states of an existing spin-orbit-doublet element** (e.g. a
   5th Ti state): add a `TiState(...)` entry to `TI_2P_STATES` in
   `src/xps_wave/reference.py` — everything downstream (`fit_ti_2p_doublets`,
   `identify_ti_doublet_peaks`, `XPSAnalyzer.analyze_ti_2p`) picks it up
   automatically via the `states=` parameter.
2. **A different doublet element entirely** (e.g. Fe 2p, Cu 2p): copy the
   pattern in `peakfit.fit_ti_2p_doublets` — it is not Ti-specific in its
   mechanics, only in which dataclass (`TiState`) and reference list it
   reads. The physical constraints that make it robust (fixed sub-shell
   degeneracy area ratio, splitting tied within a small tolerance, sigma
   bounded near a literature FWHM rather than left free) are the reusable
   part; generalize `TiState` to a generic `DoubletState` if this happens
   more than once.
3. **A new element/orbital for wide-spectrum identification only** (no
   doublet fitting, just "what's near this BE"): add a `ReferenceLine(...)`
   entry to `SURVEY_REFERENCE_TABLE` in the same file, with a source
   citation — see the existing entries for the expected level of rigor
   (real literature/NIST values, not guesses).

## Known limitations (tell the user if relevant, don't silently paper over them)

- `SURVEY_REFERENCE_TABLE` values are representative literature values
  (NIST XPS Database / Moulder Handbook), not a live database query — for
  quantitative work point the user at https://srdata.nist.gov/xps/.
- No automatic charge-referencing/calibration for arbitrary input files.
  `read_vamas` only auto-corrects the specific CasaXPS `Calib M=.. A=..`
  comment convention if present; a user's own uncalibrated CSV will need a
  manual BE shift (`Spectrum.energy += shift`) before fitting, same as any
  real XPS workflow requires (e.g. referencing adventitious C 1s to 284.8 eV).
- `--mode wide` peak search defaults (`prominence_fraction=0.05`,
  `min_distance_ev=3.0`) are tuned for a "typical" survey SNR — noisy real
  data will need `analyze_wide_spectrum(..., prominence_fraction=...)`
  tuned higher (see `scripts/analyze_wide_spectrum.py` for a worked example
  going from very noisy to clean).
- The wide-spectrum demo (`scripts/analyze_wide_spectrum.py`) uses a
  synthetic spectrum, not real data — no open, redistributable raw survey
  spectrum with Ti was found when this project was built. If the user later
  finds/provides one, prefer it for any future accuracy claims about
  wide-spectrum mode specifically.

## Where things live

Full architecture in `README.md` ("アーキテクチャ"). In short:
`src/xps_wave/{spectrum,io,background,peakfit,reference,identify,survey,pipeline,plotting}.py`,
CLI scripts in `scripts/`, tests in `tests/` (`uv run pytest -q`).
