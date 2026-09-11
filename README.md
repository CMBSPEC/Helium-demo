# Helium-demo

Helium-demo provides a compact C++ example for the neutral helium atom setup in
`Helium.vX`, together with a Python interface for interactive work and
notebook-based demonstrations. The C++ code builds the `helium_atom_demo`
executable and exercises the helium routines together with the minimal shared
numerical utilities from `Tools`.

Hydrogenic.vX is included only because the helium module uses hydrogenic
bound-bound, photoionization, and recombination helpers internally. The examples
in this repository are intentionally helium-focused.

## Project layout

| Path | Purpose |
| --- | --- |
| `main.cpp` | Example program that constructs neutral helium and prints representative atomic data. |
| `Helium.vX/HeI_Atom.*` | Core neutral helium atom classes and level/transition access. |
| `Helium.vX/HeI_Ric_Topbase.*` | TOPbase photoionization cross-section and recombination integrals. |
| `Helium.vX/He4_Quantum_Defects.*` | Quantum-defect energy helper for high-lying He I levels. |
| `Helium.vX/Helium.Data/` | Transition tables and TOPbase photoionization data used by the helium setup. |
| `Hydrogenic.vX/` | Internal dependency used by the helium implementation. |
| `Tools/` | Shared constants, integration routines, line profiles, and simple numerical helpers. |
| `python/` | Installable Python package and Jupyter notebook examples. |
| `Makefile` and `Makefile.in` | Build entry points for the C++ demo executable. |

## What the C++ example does

The demo constructs a `Gas_of_HeI_Atoms` object with ten shells, j-resolved
triplet levels, intercombination lines, and quadrupole lines enabled. It reports:

- level degeneracies, ionization energies, and excitation energies for selected
  singlet and triplet states;
- selected He I resonance, intercombination, and quadrupole transitions;
- the special He I two-photon rates stored outside the regular transition
  table;
- photoionization and recombination rates at a representative radiation
  temperature;
- TOPbase photoionization cross sections for selected levels;
- simple rescaling behavior for the fine-structure constant and electron mass.

## Dependencies

The C++ build expects:

- a C++17 compiler with OpenMP support;
- GSL headers and libraries;
- Boost headers;
- GNU `make`.

The Python package additionally uses:

- Python 3.9 or newer;
- `pybind11` and `numpy` at build time;
- optional `jupyter` and `matplotlib` for notebooks and plots.

On macOS with Homebrew, the typical compiled dependency setup is:

```sh
brew install gcc gsl boost
```

`Makefile.in` currently defaults to:

```make
CXX = /opt/homebrew/bin/g++-16
DEF_INC_PATH = /opt/homebrew/include
DEF_LIB_PATH = /opt/homebrew/lib
```

If your compiler or Homebrew prefix differs, override these variables when
invoking `make` or `pip`.

## Build and run the C++ demo

From the repository root:

```sh
make
make run
```

List available Makefile targets:

```sh
make help
```

Clean generated build products:

```sh
make clean
```

Update the `Helium.vX` and `Hydrogenic.vX` submodules from their configured
remote branches:

```sh
make update-modules
```

The executable is written to:

```sh
./helium_atom_demo
```

Intermediate C++ object and dependency files are written under:

```sh
./build/
```

## Python package and Jupyter

The Makefile prefers `/opt/homebrew/bin/python3` when it exists, then falls back
to the first `python3` on `PATH`. You can override this with `PYTHON=...` if
your Python packages live in a different environment. Avoid using plain
`PYTHON=python` from the repository root, since this repository has a `python/`
folder and some shell setups can otherwise report `Permission denied`.

On macOS with Homebrew:

```sh
brew install gcc gsl boost python-setuptools pybind11 numpy matplotlib jupyterlab
```

Build the local Python extension:

```sh
make python-build
```

Open the demo notebook:

```sh
make run-jupyter
```

`make run-jupyter` runs `make python-build` first, then opens the notebook with
an available Jupyter command and sets `PYTHONPATH` to include `python/src`. This
lets the notebook import the local `helium` extension without installing it into
JupyterLab's private Python environment.

## Use from Python

```python
import numpy as np
from helium import helium

he = helium(shells=10)
print(he.summary())
print(he.level(2, 1, S=0, J=1))
print(he.transition(2, 1, 0, 1, 1, 0, 0, 0))
print(he.two_photon_rate(2, 0, 0, 0))

nu0 = he.level(1, 0)["ionization_frequency_hz"]
nu = np.geomspace(1.001, 5.0, 100)*nu0
data = he.topbase_photoionization(1, 0, 0, 0, nu)
print(data["sigma_cm2"])
```

## Scientific background

The helium atom setup was developed for detailed calculations of cosmological
recombination. This requires neutral-helium singlet and triplet levels,
j-resolved low-lying triplet states, intercombination lines, non-dipole
transitions, photoionization cross sections, and recombination coefficients.

The transition tables are based primarily on Drake and Morton, with additional
forbidden-transition data from Lach and Pachucki. Recombination approximations
for selected levels use Benjamin, Skillman, and Smits. Photoionization data are
taken from TOPbase where available, and the hydrogenic module provides the
fallback hydrogenic machinery used internally by the helium code.

## Related Repositories

- [Helium.vX](https://github.com/CMBSPEC/Helium.vX)
- [Hydrogenic.vX](https://github.com/CMBSPEC/Hydrogenic.vX)
- [Hydrogenic-demo](https://github.com/CMBSPEC/Hydrogenic-demo)

## Related Literature And Data Sources

The code and tables are connected to several atomic-physics and recombination
references, including:

- G. W. F. Drake and D. C. Morton, "A Multiplet Table for Neutral Helium
  (^4He I) with Transition Rates", ApJS 170, 251, 2007.
- R. A. Benjamin, E. D. Skillman, and D. P. Smits, "Improving Predictions for
  Helium Emission Lines", ApJ 514, 307, 1999.
- G. Lach and K. Pachucki, "Forbidden transitions in the helium atom",
  Phys. Rev. A 64, 042510, 2001.
- W. Cunto, C. Mendoza, F. Ochsenbein, and C. J. Zeippen, "TOPbase at the CDS",
  A&A 275, L5, 1993.
- J. A. Rubino-Martin, J. Chluba, and R. A. Sunyaev, "Lines in the cosmic
  microwave background spectrum from the epoch of cosmological helium
  recombination", A&A 485, 377, 2008.
- J. Chluba and R. M. Thomas, "Towards a complete treatment of the
  cosmological recombination problem", MNRAS 412, 748, 2011.

These references are listed to document the origin and scientific context of
the data and methods used here; this repository provides examples built around
the helium atom setup.

**Acknowledgements:** This repository was made available and documented with the help of Codex. The related release work was supported in part by a grant of access to OpenAI models through the ChatGPT for Academic Researchers program.
