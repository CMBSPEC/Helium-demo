# Helium-demo Python Interface

This package exposes a small Python interface to the C++ neutral helium atom
demo. It is intended for interactive checks, plots, and notebooks, not as a
fully packaged standalone distribution.

Build the extension from the repository root:

```sh
make python-build
```

Then use it with `PYTHONPATH=python/src`, or launch the notebook through:

```sh
make run-jupyter
```

Basic use:

```python
import numpy as np
from helium import helium

he = helium(shells=10)
print(he.summary())

line = he.transition(2, 1, 0, 1, 1, 0, 0, 0)
print(line)
print(he.two_photon_rate(2, 0, 0, 0))

nu0 = he.level(1, 0)["ionization_frequency_hz"]
nu = np.geomspace(1.001, 5.0, 100)*nu0
photo = he.topbase_photoionization(1, 0, 0, 0, nu)
```

The default factory resolves `Helium.vX/Helium.Data/` relative to this source
checkout. Pass `data_path=...` to use a different helium data directory.
