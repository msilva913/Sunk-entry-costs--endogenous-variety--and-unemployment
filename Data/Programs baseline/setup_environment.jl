# setup_environment.jl
# =============================================================================
# One-time setup of the project environment for the model programs.
# Creates Project.toml / Manifest.toml in this folder and installs every
# third-party package the model code loads. Run from this folder:
#     julia --startup-file=no setup_environment.jl
# On another machine with the Manifest committed, `Pkg.instantiate()` reproduces
# the same versions instead.
#
# SymPy needs a Python with sympy. PYTHON = "" makes PyCall use a private Conda
# distribution, so the setup does not depend on any system Python.
# =============================================================================
import Pkg
cd(@__DIR__)
Pkg.activate(".")

ENV["PYTHON"] = ""

pkgs = ["MKL", "DataFrames", "SymPy", "PyCall", "BenchmarkTools", "Parameters",
        "Distributions", "StatsBase", "CSV", "NLsolve", "Roots", "Optim",
        "LeastSquaresOptim", "PrettyPrinting", "LaTeXStrings", "Plots", "PGFPlotsX",
        "KernelDensity", "ShiftedArrays", "MappedArrays", "TexTables", "TypedTables", "GLM"]

Pkg.add(pkgs)
Pkg.build("PyCall")
Pkg.precompile()

# Load SymPy once so Conda installs the Python sympy module now, not mid-run.
@eval using SymPy
println("sympy version: ", SymPy.sympy.__version__)
println("Environment ready: ", Base.active_project())
