from __future__ import annotations

import os
import pathlib
import shutil
import sys

import numpy
import pybind11
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


ROOT = pathlib.Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent


def _homebrew_prefix() -> pathlib.Path:
    env_prefix = os.environ.get("HOMEBREW_PREFIX")
    for prefix in [env_prefix, "/opt/homebrew", "/usr/local"]:
        if prefix and pathlib.Path(prefix).exists():
            return pathlib.Path(prefix)
    return pathlib.Path("/opt/homebrew")


def _configure_default_compiler() -> None:
    if os.environ.get("CXX"):
        return

    for compiler in [
        "/opt/homebrew/bin/g++-16",
        "/opt/homebrew/bin/g++-15",
        "/opt/homebrew/bin/g++-14",
        "/usr/local/bin/g++-16",
        "/usr/local/bin/g++-15",
        "/usr/local/bin/g++-14",
        shutil.which("g++"),
    ]:
        if compiler and pathlib.Path(compiler).exists():
            os.environ["CXX"] = compiler
            os.environ.setdefault("CC", compiler)
            return


_configure_default_compiler()

brew_prefix = _homebrew_prefix()
default_include = brew_prefix / "include"
default_lib = brew_prefix / "lib"

gsl_include = pathlib.Path(os.environ.get("GSL_INC_PATH", str(default_include)))
gsl_lib = pathlib.Path(os.environ.get("GSL_LIB_PATH", str(default_lib)))
boost_include = pathlib.Path(os.environ.get("BOOST_INC_PATH", str(default_include)))

include_dirs = [
    pybind11.get_include(),
    numpy.get_include(),
    str(REPO_ROOT / "Helium.vX"),
    str(REPO_ROOT / "Hydrogenic.vX"),
    str(REPO_ROOT / "Tools" / "Definitions"),
    str(REPO_ROOT / "Tools" / "Integration"),
    str(REPO_ROOT / "Tools" / "Line_profiles"),
    str(REPO_ROOT / "Tools" / "Simple_routines"),
    str(gsl_include),
    str(boost_include),
]

sources = [
    "src/helium/_core.cpp",
    "../Helium.vX/HeI_Atom.cpp",
    "../Helium.vX/HeI_Ric_Topbase.cpp",
    "../Helium.vX/He4_Quantum_Defects.cpp",
    "../Hydrogenic.vX/Atom.cpp",
    "../Hydrogenic.vX/Oscillator_strength.cpp",
    "../Hydrogenic.vX/Photoionization_cross_section.cpp",
    "../Hydrogenic.vX/Quadrupole_lines_HI.cpp",
    "../Hydrogenic.vX/Rec_Phot_BB.cpp",
    "../Hydrogenic.vX/Recombination_Integrals.cpp",
    "../Tools/Simple_routines/routines.cpp",
    "../Tools/Simple_routines/File.cpp",
    "../Tools/Integration/Chebyshev_Int.cpp",
    "../Tools/Integration/Integration_routines.GSL.cpp",
    "../Tools/Integration/Patterson.cpp",
    "../Tools/Line_profiles/Voigtprofiles.cpp",
]


class BuildExt(build_ext):
    def build_extensions(self) -> None:
        cxx = os.environ.get("CXX")
        if cxx and self.compiler.compiler_type == "unix":
            for attr in ["compiler_so", "compiler_cxx", "linker_so"]:
                value = getattr(self.compiler, attr, None)
                if isinstance(value, list) and value:
                    value[0] = cxx

        for ext in self.extensions:
            ext.extra_compile_args.extend([
                "-std=c++17",
                "-O2",
                "-Wall",
                "-pthread",
                "-fopenmp",
                "-DOPENMP_ACTIVATED",
            ])
            ext.extra_link_args.extend([
                "-pthread",
                "-fopenmp",
            ])

            if sys.platform == "darwin":
                ext.extra_compile_args.append("-Wno-unused-local-typedefs")

        super().build_extensions()


extension = Extension(
    "helium._core",
    sources=sources,
    include_dirs=include_dirs,
    define_macros=[("HEIDATADIR", '"./Helium.vX/"')],
    library_dirs=[str(gsl_lib)],
    libraries=["gsl", "gslcblas"],
    language="c++",
)


setup(
    name="helium-demo",
    version="0.1.0",
    description="Python bindings for the Helium-demo C++ atom routines",
    ext_modules=[extension],
    cmdclass={"build_ext": BuildExt},
)
