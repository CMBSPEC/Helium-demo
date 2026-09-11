"""Python interface for the Helium-demo C++ atom routines."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from ._core import (
    HeliumAtom,
    heI_triplet_two_photon_2s_1s_rate_per_s,
    heI_two_photon_2s_1s_rate_per_s,
)


Atom = HeliumAtom
PathLike = Union[str, Path]
HEI_TRIPLET_TWO_PHOTON_2S_1S_RATE_PER_S = heI_triplet_two_photon_2s_1s_rate_per_s
HEI_TWO_PHOTON_2S_1S_RATE_PER_S = heI_two_photon_2s_1s_rate_per_s


def default_data_path() -> Path:
    """Return the default Helium.vX data directory for this source checkout."""
    return Path(__file__).resolve().parents[3] / "Helium.vX" / "Helium.Data"


def helium(
    shells: int = 10,
    *,
    j_resolved_shells: int = 10,
    quadrupole_shells: int = 10,
    intercombination_shells: int = 10,
    data_path: Optional[PathLike] = None,
    initialize_photoionization: bool = True,
    message_level: int = -1,
) -> HeliumAtom:
    """Create a neutral helium atom with the defaults used by the C++ demo."""
    if data_path is None:
        data_path = default_data_path()

    return HeliumAtom(
        shells,
        j_resolved_shells,
        quadrupole_shells,
        intercombination_shells,
        str(data_path),
        initialize_photoionization,
        message_level,
    )


__all__ = [
    "Atom",
    "HEI_TRIPLET_TWO_PHOTON_2S_1S_RATE_PER_S",
    "HEI_TWO_PHOTON_2S_1S_RATE_PER_S",
    "HeliumAtom",
    "default_data_path",
    "heI_triplet_two_photon_2s_1s_rate_per_s",
    "heI_two_photon_2s_1s_rate_per_s",
    "helium",
]
