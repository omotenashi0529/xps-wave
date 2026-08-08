"""Loaders for XPS measurement files: generic 2-column text/CSV and VAMAS (.vms, ISO 14976)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from xps_wave.spectrum import Spectrum, ke_to_be

_CASA_CALIB_RE = re.compile(r"Calib\s+M\s*=\s*([\d.+-]+)\s+A\s*=\s*([\d.+-]+)\s+BE\s+ADD", re.IGNORECASE)


def read_two_column(
    path: str | Path,
    energy_col: int = 0,
    intensity_col: int = 1,
    delimiter: str | None = None,
    skip_header: int = 0,
    energy_is_kinetic: bool = False,
    photon_energy: float | None = None,
    name: str | None = None,
) -> Spectrum:
    """Read a generic two-column energy/intensity text or CSV file.

    Most XPS software (CasaXPS, Avantage, Multipak, XPSPEAK) can export a region or
    survey scan as a plain two-column ASCII/CSV file; this is the lowest common
    denominator format so it is the default input path for real measurement data.
    """
    path = Path(path)
    data = np.genfromtxt(path, delimiter=delimiter, skip_header=skip_header)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"expected at least 2 numeric columns in {path}, got shape {data.shape}")
    energy = data[:, energy_col]
    intensity = data[:, intensity_col]
    valid = ~(np.isnan(energy) | np.isnan(intensity))
    energy, intensity = energy[valid], intensity[valid]

    if energy_is_kinetic:
        if photon_energy is None:
            raise ValueError("photon_energy is required when energy_is_kinetic=True")
        energy = ke_to_be(energy, photon_energy)

    return Spectrum(energy=energy, intensity=intensity, name=name or path.stem, source=str(path))


def read_vamas(path: str | Path, block_index: int = 0, name: str | None = None) -> Spectrum:
    """Read a VAMAS (ISO 14976) .vms file, as exported by CasaXPS/Kratos/PHI tools.

    VAMAS blocks store the abscissa as kinetic energy by convention; this loader
    converts to binding energy using the block's own photon energy
    (`analysis_source_characteristic_energy`) so all Spectrum objects share the
    same BE convention regardless of source format.
    """
    from vamas import Vamas

    path = Path(path)
    vms = Vamas(str(path))
    if not vms.blocks:
        raise ValueError(f"no data blocks found in VAMAS file {path}")
    block = vms.blocks[block_index]

    intensity_var = next(
        (cv for cv in block.corresponding_variables if "intens" in cv.label.lower()),
        block.corresponding_variables[0],
    )
    intensity = np.asarray(intensity_var.y_values, dtype=float)

    # VAMAS's "number of ordinate values" counts every corresponding variable's
    # values together (points * num_corresponding_variables), not just points.
    n = len(intensity)
    x_start = block.x_start
    x_step = block.x_step
    kinetic_energy = x_start + x_step * np.arange(n)

    hv = block.analysis_source_characteristic_energy
    x_units = (block.x_units or "").lower()
    if "kinetic" in (block.x_label or "").lower() or x_units in ("ev",):
        # VAMAS almost always stores kinetic energy; ke_to_be is a no-op-safe
        # conversion since it only needs hv, which every real block carries.
        energy = ke_to_be(kinetic_energy, hv) if hv else kinetic_energy
    else:
        energy = kinetic_energy

    # CasaXPS often stores the raw, as-measured KE axis in the VAMAS block
    # itself and keeps its charge-referencing calibration only as a text
    # annotation in block_comment ("Calib M = <measured BE> A = <accepted BE>
    # BE ADD"). Without applying that shift, this block's stored BE axis
    # disagrees with the CasaXPS-fitted component positions (also stored as
    # text) by exactly that offset. Apply it automatically when present so
    # the returned Spectrum is in the same, literature-calibrated BE frame
    # as any embedded fit results.
    calibration_shift_eV = 0.0
    calib_match = _CASA_CALIB_RE.search(block.block_comment or "")
    if calib_match:
        measured, accepted = float(calib_match.group(1)), float(calib_match.group(2))
        calibration_shift_eV = accepted - measured
        energy = energy + calibration_shift_eV

    metadata = {
        "sample_identifier": block.sample_identifier,
        "species_label": block.species_label,
        "transition_or_charge_state_label": block.transition_or_charge_state_label,
        "photon_energy_eV": hv,
        "block_comment": block.block_comment,
        "technique": block.technique,
        "calibration_shift_eV": calibration_shift_eV,
    }
    label = name or block.sample_identifier or path.stem
    return Spectrum(energy=energy, intensity=intensity, name=label, source=str(path), metadata=metadata)
