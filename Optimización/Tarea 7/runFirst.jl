# Run this first to ensure all requirements are satisfied. You will need to accept the download of MNIST dataset and it might take a minute.

using Pkg
Pkg.update()
Pkg.add("MLDatasets")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")

using MLDatasets

MNIST()