"""Spectrum container and energy-scale conversions shared by the whole package."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def ke_to_be(kinetic_energy: np.ndarray, photon_energy: float, work_function: float = 0.0) -> np.ndarray:
    """Convert kinetic energy (eV) to binding energy (eV): BE = hv - KE - work_function."""
    return photon_energy - np.asarray(kinetic_energy, dtype=float) - work_function


def be_to_ke(binding_energy: np.ndarray, photon_energy: float, work_function: float = 0.0) -> np.ndarray:
    """Convert binding energy (eV) to kinetic energy (eV): KE = hv - BE - work_function."""
    return photon_energy - np.asarray(binding_energy, dtype=float) - work_function


@dataclass
class Spectrum:
    """A single XPS energy/intensity trace plus provenance metadata.

    energy is always stored as binding energy (eV), ascending order, since that is
    the convention chemists read peak tables in and the convention the reference
    binding-energy database uses.
    """

    energy: np.ndarray
    intensity: np.ndarray
    name: str = "spectrum"
    source: str = ""
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.energy = np.asarray(self.energy, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
        if self.energy.shape != self.intensity.shape:
            raise ValueError("energy and intensity arrays must have the same shape")
        if self.energy.size >= 2 and self.energy[0] > self.energy[-1]:
            self.energy = self.energy[::-1].copy()
            self.intensity = self.intensity[::-1].copy()

    def window(self, low: float, high: float, name: str | None = None) -> "Spectrum":
        """Return a new Spectrum restricted to [low, high] eV binding energy."""
        lo, hi = min(low, high), max(low, high)
        mask = (self.energy >= lo) & (self.energy <= hi)
        if not np.any(mask):
            raise ValueError(f"no data points in window [{lo}, {hi}] eV for '{self.name}'")
        return Spectrum(
            energy=self.energy[mask],
            intensity=self.intensity[mask],
            name=name or f"{self.name} [{lo:g}-{hi:g} eV]",
            source=self.source,
            metadata=dict(self.metadata),
        )

    def __len__(self) -> int:
        return self.energy.size
