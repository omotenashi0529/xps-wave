"""Download the real Ti 2p reference spectrum used for accuracy evaluation.

Source: "Mixed Titanium Sample, Ti 2p" VAMAS file, an expert-fitted (CasaXPS)
real measured spectrum published on the XPS Reference Pages for Titanium
(http://www.xpsfitting.com/2008/09/titanium.html), maintained by
M.C. Biesinger, Surface Science Western, Western University. The file
carries the CasaXPS component fit (chemical-state assignment, binding
energies, spin-orbit splitting, area ratios) as embedded text, which
`scripts/evaluate_accuracy.py` uses as independently expert-derived ground
truth to score this project's own peak-fitting pipeline against.

This is downloaded at run time (not vendored in the repo) and is intended
for research/educational, non-commercial use consistent with the source
page's stated purpose (a community reference resource for XPS peak
fitting) - re-check http://www.xpsfitting.com for current terms before any
other use.
"""

from __future__ import annotations

from pathlib import Path

import requests

GOOGLE_DRIVE_FILE_ID = "1AA-XCdgku2k1jhCXSTSfr1JgC9kJQtHy"
SOURCE_PAGE = "http://www.xpsfitting.com/2008/09/titanium.html"
DEST = Path(__file__).resolve().parent.parent / "data" / "raw" / "mixed_ti_2p.vms"


def fetch(dest: Path = DEST, force: bool = False) -> Path:
    if dest.exists() and not force:
        print(f"already present: {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://drive.google.com/uc?export=download&id={GOOGLE_DRIVE_FILE_ID}"
    resp = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
    resp.raise_for_status()
    if not resp.content.startswith(b"VAMAS"):
        raise RuntimeError(
            "downloaded content does not look like a VAMAS file - "
            f"the Google Drive link may have changed; check {SOURCE_PAGE}"
        )
    dest.write_bytes(resp.content)
    print(f"downloaded {len(resp.content)} bytes -> {dest}")
    print(f"source: {SOURCE_PAGE}")
    return dest


if __name__ == "__main__":
    fetch()
