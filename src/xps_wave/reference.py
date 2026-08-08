"""Reference binding-energy database used for element / chemical-state identification.

Two tiers:

* ``TI_2P_STATES`` - high-precision Ti 2p chemical-state parameters (binding
  energy, spin-orbit splitting, 2p3/2:2p1/2 area ratio) taken from Biesinger et
  al., Appl. Surf. Sci. 257 (2010) 887-898 ("Resolving surface chemical states
  in XPS analysis of first row transition metals, oxides and hydroxides: Sc,
  Ti, V, Cu and Zn"), cross-checked against xpsfitting.com's Titanium
  reference page (M.C. Biesinger, Surface Science Western) and against the
  real, expert-fitted "Mixed Titanium Sample" VAMAS file bundled with this
  project's accuracy evaluation (see ``scripts/fetch_reference_data.py``).
  These are the values ``peakfit.fit_ti_2p`` uses to build constrained
  spin-orbit-doublet models.

* ``SURVEY_REFERENCE_TABLE`` - a broader, lower-precision table of principal
  photoelectron lines for ~30 common elements, used by the wide-spectrum
  (survey) identification pass. Values are the commonly tabulated ones found
  in the NIST XPS Database (https://srdata.nist.gov/xps/) and Moulder et al.,
  *Handbook of X-ray Photoelectron Spectroscopy* (Physical Electronics, 1995).
  They are typical literature values, not a substitute for the NIST database
  for rigorous quantitative work - see the accuracy-evaluation report for how
  well they perform against a real measured spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class TiState:
    state: str
    formal_charge: str
    be_2p32: float
    splitting: float
    area_ratio_32_to_12: float
    fwhm_2p32: float
    fwhm_2p12: float
    lineshape: str
    source: str


TI_2P_STATES: list[TiState] = [
    TiState("Ti metal", "Ti(0)", 453.8, 6.05, 2.0, 0.55, 0.75, "LA(1.1,5,7) (metallic, asymmetric)",
            "Biesinger et al. 2010; xpsfitting.com Titanium page"),
    TiState("TiO", "Ti(II)", 455.4, 5.72, 2.0, 1.6, 2.0, "GL(30)",
            "Biesinger et al. 2010; xpsfitting.com Titanium page"),
    TiState("Ti2O3", "Ti(III)", 457.2, 5.72, 2.0, 1.6, 2.0, "GL(30)",
            "Biesinger et al. 2010; xpsfitting.com Titanium page"),
    TiState("TiO2 (rutile/anatase)", "Ti(IV)", 458.6, 5.72, 2.0, 1.0, 1.8, "GL(30)",
            "Biesinger et al. 2010; NIST XPS Database; xpsfitting.com Titanium page"),
]


@dataclass(frozen=True)
class ReferenceLine:
    element: str
    orbital: str
    state: str
    be: float
    source: str = "NIST XPS Database / Moulder Handbook (typical literature value)"


SURVEY_REFERENCE_TABLE: list[ReferenceLine] = [
    ReferenceLine("Li", "1s", "metal/Li2O", 55.0),
    ReferenceLine("C", "1s", "adventitious/C-C", 284.8),
    ReferenceLine("C", "1s", "C-O/C=O", 286.5),
    ReferenceLine("N", "1s", "amine/nitride", 399.5),
    ReferenceLine("O", "1s", "lattice oxide", 530.0),
    ReferenceLine("O", "1s", "hydroxide/adsorbed", 532.0),
    ReferenceLine("F", "1s", "fluoride", 684.8),
    ReferenceLine("Na", "1s", "Na+", 1071.7),
    ReferenceLine("Mg", "2p", "metal", 49.8),
    ReferenceLine("Mg", "2p", "MgO", 50.3),
    ReferenceLine("Al", "2p", "metal", 72.8),
    ReferenceLine("Al", "2p", "Al2O3", 74.4),
    ReferenceLine("Si", "2p", "elemental", 99.3),
    ReferenceLine("Si", "2p", "SiO2", 103.3),
    ReferenceLine("P", "2p", "phosphate", 133.0),
    ReferenceLine("S", "2p3/2", "sulfide", 161.5),
    ReferenceLine("S", "2p3/2", "sulfate", 168.9),
    ReferenceLine("Cl", "2p3/2", "chloride", 198.7),
    ReferenceLine("K", "2p3/2", "K+", 293.6),
    ReferenceLine("Ca", "2p3/2", "Ca2+", 346.4),
    ReferenceLine("Ti", "2p3/2", "Ti metal", 453.8),
    ReferenceLine("Ti", "2p3/2", "TiO", 455.4),
    ReferenceLine("Ti", "2p3/2", "Ti2O3", 457.2),
    ReferenceLine("Ti", "2p3/2", "TiO2", 458.6),
    ReferenceLine("V", "2p3/2", "metal/oxide", 512.9),
    ReferenceLine("Cr", "2p3/2", "metal", 574.1),
    ReferenceLine("Cr", "2p3/2", "Cr2O3", 576.5),
    ReferenceLine("Mn", "2p3/2", "metal/oxide", 638.8),
    ReferenceLine("Fe", "2p3/2", "metal", 706.8),
    ReferenceLine("Fe", "2p3/2", "Fe2O3/FeOOH", 710.9),
    ReferenceLine("Co", "2p3/2", "metal/oxide", 778.1),
    ReferenceLine("Ni", "2p3/2", "metal", 852.6),
    ReferenceLine("Ni", "2p3/2", "NiO", 855.7),
    ReferenceLine("Cu", "2p3/2", "metal/Cu2O", 932.6),
    ReferenceLine("Cu", "2p3/2", "CuO", 933.6),
    ReferenceLine("Zn", "2p3/2", "ZnO/metal", 1021.8),
    ReferenceLine("Ga", "3d", "metal/oxide", 19.5),
    ReferenceLine("Ge", "3d", "metal/oxide", 29.5),
    ReferenceLine("As", "3d", "metal/oxide", 41.5),
    ReferenceLine("Se", "3d", "elemental", 55.0),
    ReferenceLine("Br", "3d", "bromide", 68.6),
    ReferenceLine("Zr", "3d5/2", "ZrO2", 182.5),
    ReferenceLine("Nb", "3d5/2", "Nb2O5", 207.2),
    ReferenceLine("Mo", "3d5/2", "metal", 227.9),
    ReferenceLine("Mo", "3d5/2", "MoO3", 232.6),
    ReferenceLine("Ag", "3d5/2", "metal", 368.3),
    ReferenceLine("Sn", "3d5/2", "SnO2", 486.6),
    ReferenceLine("Sb", "3d5/2", "oxide", 530.0),
    ReferenceLine("Ba", "3d5/2", "BaO/BaCO3", 780.0),
    ReferenceLine("La", "3d5/2", "La2O3", 835.4),
    ReferenceLine("W", "4f7/2", "metal", 31.4),
    ReferenceLine("W", "4f7/2", "WO3", 35.8),
    ReferenceLine("Pt", "4f7/2", "metal", 71.2),
    ReferenceLine("Au", "4f7/2", "metal", 84.0),
    ReferenceLine("Pb", "4f7/2", "PbO", 138.8),
    ReferenceLine("Bi", "4f7/2", "Bi2O3", 158.9),
]
