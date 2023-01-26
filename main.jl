###
# Main file for running the macro-micro decomposition
# for the radiative transfer equation. 
#
# developed by Chinmay Patwardhan, 2023
###

include("settings.jl")
include("solver.jl")

s = Settings(1001,800); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv

# run solver for various Settings

s.Tend = 1.0;
# s.dt = s.dt/2
Solver = solver(s);
@time t, rho1, g1 = solveFullProblem(Solver);

using PyPlot


fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(Solver.x, rho1)
fig.canvas.draw()