###
# Main file for running the macro-micro decomposition
# for the radiative transfer equation. 
#
# developed by Chinmay Patwardhan, 2023
###

include("settings.jl")
include("solver.jl")

using PyPlot
using BenchmarkTools

<<<<<<< HEAD
close("all")

s = Settings(2002,500,1.0,"hyperbolic"); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv
# The input parameter for setting the matrices for angular discretisation are the same for both Pn and Sn solver
=======
s = Settings(1001,500,1.0,"mixed"); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv
>>>>>>> fc541698d1109c0ba26900cb19ae016a46ef7fcf

# run solver for various Settings

s.Tend = 1.0;
# s.dt = s.dt/10; # For the kinetic regime smaller step size than the one selected is required
Solver = solver(s);
@time t, rho1, g1 = solveFullProblem_Pn(Solver);

<<<<<<< HEAD
s.Tend = 0.2;
=======
s.Tend = 5.0;
>>>>>>> fc541698d1109c0ba26900cb19ae016a46ef7fcf
# s.dt = s.dt/10; # For the kinetic regime smaller step size than the one selected is required
# Solver = solver(s);
# @time t, rho2 = solveLimitingDiff(Solver);


<<<<<<< HEAD
fig1, ax1 = subplots(figsize=(15, 12), dpi=100);
ax1.plot(Solver.x, rho1, label="Macro-Micro");
# ax1.plot(Solver.x, rho2, label="Diffusion limit");
ax1.legend();
fig1.canvas.draw();

# fig2, ax2 = subplots(figsize=(15, 12), dpi=100);
# ax2.semilogy(Solver.x, rho2-rho1, label="Error");
# ax2.legend();
# fig2.canvas.draw();
=======
fig, ax = subplots(figsize=(15, 12), dpi=100);
ax.plot(Solver.x, rho1, label="Macro-Micro");
# ax.plot(Solver.x, rho2, label="Diffusion limit");
ax.legend();
fig.canvas.draw();

# fig, ax = subplots(figsize=(15, 12), dpi=100);
# ax.semilogy(Solver.x, rho2-rho1, label="Error");
# ax.legend();
# fig.canvas.draw();
>>>>>>> fc541698d1109c0ba26900cb19ae016a46ef7fcf
