##################
### Here are the codes for installing the old version. You shall create a new environment in VScode for them.
### You can check the version logs at https://juliahub.com/ui/Home.
##################

using Pkg
Pkg.add(["DataFrames", "Optim"])
Pkg.add([
    Pkg.PackageSpec(name="Random"),  # Standard library
    Pkg.PackageSpec(name="Statistics"),  # Standard library
    Pkg.PackageSpec(name="SparseArrays"),  # Standard library
    Pkg.PackageSpec(name="LinearAlgebra"),  # Standard library
    Pkg.PackageSpec(name="SymPy", version="1.1.5"),
    Pkg.PackageSpec(name="BenchmarkTools", version="1.3.2"),
    Pkg.PackageSpec(name="Parameters", version="0.12.3"),
    Pkg.PackageSpec(name="Distributions", version="0.25.50"),
    Pkg.PackageSpec(name="StatsBase", version="0.33.15"),
    Pkg.PackageSpec(name="MKL", version="0.5.0"),
    Pkg.PackageSpec(name="Plots", version="1.38.0"),
    Pkg.PackageSpec(name="PyPlot", version="2.11.0"),
    Pkg.PackageSpec(name="CSV", version="0.10.5"),
    Pkg.PackageSpec(name="Roots", version="1.4.0"),
    Pkg.PackageSpec(name="LeastSquaresOptim", version="0.8.2"),
    Pkg.PackageSpec(name="PrettyPrinting", version="0.2.0"),
    Pkg.PackageSpec(name="LaTeXStrings", version="1.3.0"),
    Pkg.PackageSpec(name="KernelDensity", version="0.6.3"),
    Pkg.PackageSpec(name="ShiftedArrays", version="1.0.0"),
    Pkg.PackageSpec(name="QuantEcon", version="0.16.4"),
    Pkg.PackageSpec(name="MappedArrays", version="0.4.1"),
    Pkg.PackageSpec(name="TexTables", version="0.2.6"),
    Pkg.PackageSpec(name="TypedTables", version="1.4.1"),
    Pkg.PackageSpec(name="GLM", version="1.5.1"),
    Pkg.PackageSpec(name="Pandas", version="1.5.0"),
])
