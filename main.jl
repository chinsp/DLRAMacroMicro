###
# Main file for running the macro-micro decomposition
# for the radiative transfer equation. 
#
# developed by Chinmay Patwardhan, 2023
###

include("settings.jl")
include("solver.jl")

s = Settings(); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv

# run solver for various Settings
