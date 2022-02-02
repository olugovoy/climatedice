# this script installs all required packages to run DICEinJulia.jl
using Pkg
println("Installation of required packages to solve DICE model")
Pkg.update()
if ! in("RData",keys(Pkg.installed())) Pkg.add("RData") end
if ! in("CSV",keys(Pkg.installed())) Pkg.add("CSV") end
if ! in("JuMP",keys(Pkg.installed())) Pkg.add("JuMP") end
if ! in("Ipopt",keys(Pkg.installed())) Pkg.add("Ipopt") end
if ! in("DataFrames",keys(Pkg.installed())) Pkg.add("DataFrames") end
if ! in("Gadfly",keys(Pkg.installed())) Pkg.add("Gadfly") end
@warn("Packages installed!")
