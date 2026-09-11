#============================================================================================
# Small wrapper for the Helium atom demo.
#============================================================================================

TOOLS_DIR = ./Tools
HYDROGENIC_DIR = ./Hydrogenic.vX
HELIUM_DIR = ./Helium.vX
PYTHON ?= $(shell command -v /opt/homebrew/bin/python3 2>/dev/null || command -v python3)
NOTEBOOK ?= python/notebooks/helium_python_demo.ipynb
PYTHONPATH := $(CURDIR)/python/src$(if $(PYTHONPATH),:$(PYTHONPATH))

SUBMAKE = $(MAKE) -f Makefile.in TOOLS_DIR=$(TOOLS_DIR) HYDROGENIC_DIR=$(HYDROGENIC_DIR) HELIUM_DIR=$(HELIUM_DIR)

.PHONY: all help h usage list targets --help -h run update-modules python-check python-build python py run-jupyter jupyter python-clean clean tidy

all:
	$(SUBMAKE) all

help:
	@echo " Helium-demo targets:"
	@echo "  make                 Build the C++ demo"
	@echo "  make help            Show this help message"
	@echo "  make h               Alias for help"
	@echo "  make usage           Alias for help"
	@echo "  make list            Alias for help"
	@echo "  make targets         Alias for help"
	@echo "  make run             Build and run the C++ demo"
	@echo "  make python-build    Build the Python extension in-place"
	@echo "  make python          Alias for python-build"
	@echo "  make py              Alias for python-build"
	@echo "  make run-jupyter     Build the Python extension and open the notebook"
	@echo "  make jupyter         Alias for run-jupyter"
	@echo "  make python-clean    Remove Python build products"
	@echo "  make clean           Remove C++ and Python build products"
	@echo "  make update-modules"
	@echo "                       Update Helium.vX and Hydrogenic.vX submodules"
	@echo ""
	@echo " Notes:"
	@echo "  make --help and make -h are GNU make options."
	@echo "  Use make help for project-specific help."

h usage list targets --help -h: help

run:
	$(SUBMAKE) run

update-modules:
	$(SUBMAKE) update-modules

python-check:
	@$(PYTHON) -c 'import importlib.util, sys; missing=[m for m in ("setuptools", "numpy", "pybind11") if importlib.util.find_spec(m) is None]; sys.exit("Missing Python build packages: " + ", ".join(missing) + ". Install with: brew install python-setuptools pybind11 numpy" if missing else 0)'

python-build: python-check
	cd python && $(PYTHON) setup.py build_ext --inplace

python py: python-build

run-jupyter: python-build
	@if command -v jupyter >/dev/null; then \
		PYTHONPATH="$(PYTHONPATH)" jupyter lab $(NOTEBOOK); \
	elif command -v jupyter-lab >/dev/null; then \
		PYTHONPATH="$(PYTHONPATH)" jupyter-lab $(NOTEBOOK); \
	elif command -v jupyter-notebook >/dev/null; then \
		PYTHONPATH="$(PYTHONPATH)" jupyter-notebook $(NOTEBOOK); \
	else \
		echo "Missing JupyterLab. Install with: brew install jupyterlab"; \
		exit 1; \
	fi

jupyter: run-jupyter

python-clean:
	rm -rf python/build python/src/helium/_core*.so python/src/helium_demo.egg-info python/src/helium/__pycache__ python/__pycache__ python/notebooks/.ipynb_checkpoints

clean: python-clean
	$(SUBMAKE) clean

tidy: clean

#============================================================================================
#============================================================================================
