# this script installs all required packages to run DICEinJulia.jl
println("Installation of required packages to solve DICE model")
if Pkg.installed("RData") == nothing Pkg.add("RData") end
if Pkg.installed("CSV") == nothing Pkg.add("CSV") end
if Pkg.installed("JuMP") == nothing Pkg.add("JuMP") end
if Pkg.installed("Ipopt") == nothing Pkg.add("Ipopt") end
if Pkg.installed("DataFrames") == nothing Pkg.add("DataFrames") end
if Pkg.installed("Gadfly") == nothing Pkg.add("Gadfly") end
Pkg.update()
warn("Now restart Julia!")
