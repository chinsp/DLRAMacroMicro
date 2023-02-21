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

## The solver type is a parameter that is beign introduced to control the stencil matrices that are being used in the solvers. 
# Default is set to 3 i.e. Pn solver with macro-micro decomposition
# 0 - Sn solver for kinetic equation without macro-micro decomposition
# 1 - Sn solver for kinetic equation with macro-micro decomposition
# 2 - Pn solver for kinetic equation without macro-micro decomposition
# 3 - Pn solver for kinetic equation with macro-micro decomposition

SolverType = 1; 

s = Settings(1001,200,1.0,"hyperbolic",SolverType); # Give the number of discretisation points for spatial domain and velocity domain as input i.e., Nx and Nv
# The input parameter for setting the matrices for angular discretisation are the same for both Pn and Sn solver
s.Tend = 1.0;
Solver = solver(s);
@time t, rho1, g1 = solveMMDLRA_SnIMEX(Solver);
# run solver for various Settings



SolverType = 1; 

s1 = Settings(1001,20,1.0,"hyperbolic",SolverType);
s1.Tend = 1.0;
Solver = solver(s1);
# @time t, rho1, g1 = solveMMDLRA_SnIMEX(Solver);
# @time t, rho2 = solveLimitingDiff(Solver);
@time t, rho,g = solveFullProblem_Sn(Solver);


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
# ax1.plot(s1.x,Manufactured1D_rho(s1,s1.Tend,s1.x), label="Exact");
# ax1.plot(Solver.x, 0.5*g_SN*Solver.w, label="SN w/o Macro-Micro");
ax1.plot(Solver.x, rho, label="SN for MM");
ax1.plot(Solver.x, rho1, label="UI with SN for MM");
# ax1.plot(Solver.x, rho2, label="Diffusion limit");
# ax1.set_title("First order upwind scheme, CFL = 0.25")
ax1.legend();
fig1.canvas.draw();

# fig2, ax2 = subplots(figsize=(15, 12), dpi=100);
# ax2.plot(Solver.x, broadcast(abs,rho2-rho1), label="Error");
# ax2.legend();
# fig2.canvas.draw();