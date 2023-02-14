###
# Main file for running the macro-micro decomposition
# for the radiative transfer equation. 
#
# developed by Chinmay Patwardhan, 2023
###

include("settings.jl")
include("solver.jl")

using PyPlot
using DelimitedFiles
using BenchmarkTools

close("all")

s = Settings(2001,200,1.0,"hyperbolic"); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv
# The input parameter for setting the matrices for angular discretisation are the same for both Pn and Sn solver

# run solver for various Settings

s.Tend = 1.0;
# s.dt = s.dt/10; # For the kinetic regime smaller step size than the one selected is required
Solver = solver(s);
@time t, rho1, g1 = solveFullProblem_Sn(Solver);

s.Tend = 0.2;
# s.dt = s.dt/10; # For the kinetic regime smaller step size than the one selected is required
# Solver = solver(s);
# @time t, rho2 = solveLimitingDiff(Solver);

# read reference solution
v = readdlm("PlaneSourceRaw", ',')
uEx = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx[i] = 0.0;
    else
        uEx[i] = Float64(v[i])
    end
end
x = collect(range(-1.5,1.5,length=(2*length(v)-1)));
uEx = [uEx[end:-1:2];uEx]

fig1, ax1 = subplots(figsize=(15, 12), dpi=100);
ax1.plot(x,uEx, label="Exact");
ax1.plot(Solver.x, rho1, label="Macro-Micro");
# ax1.plot(Solver.x, rho2, label="Diffusion limit");
ax1.legend();
fig1.canvas.draw();

# fig2, ax2 = subplots(figsize=(15, 12), dpi=100);
# ax2.semilogy(Solver.x, rho2-rho1, label="Error");
# ax2.legend();
# fig2.canvas.draw();