__precompile__
include("quadrature.jl")

using ProgressMeter
using LinearAlgebra

struct solver
    # Spatial grid of cell vertices
    x::Array{Float64,1};
    xMid::Array{Float64,1};

    # Solver settings
    settings::Settings;

    # Angular discretisation
    w::Array{Float64,1};
    v::Array{Float64,2};
    vp::Array{Float64,2};
    vm::Array{Float64,2};

    # Stencil matrices for spatial discretisation
    Dp::Array{Float64,2};
    Dm::Array{Float64,2};
    Dc::Array{Float64,2};
    Dcx::Array{Float64,2};

    # Physical parameters
    sigmaA::Float64;
    sigmaS::Float64;

    # Constructor
    function solver(settings)
        x = settings.x;
        xMid = settings.xMid;

        nx = settings.Nx;
        nxC = settings.NxC;

        # Setting up the weights vector
        Nv = settings.Nv;

        quad = Quadrature(Nv,"Gauss");

        w = quad.w[end:-1:1];
        v = Diagonal(quad.v[end:-1:1]);
        
        vp = Diagonal(zeros(Float64,Nv));
        vm = Diagonal(zeros(Float64,Nv));
        for q = 1:Nv
            vp[q,q] = max(v[q,q],0);
            vm[q,q] = min(v[q,q],0);
        end

        dx = settings.dx;
        
        Dp = zeros(Float64,nxC,nxC);
        Dm = zeros(Float64,nxC,nxC);
        Dc = zeros(Float64,nxC,nx);
        Dcx = zeros(Float64,nx,nxC);
        
        # Currently running a second order upwind scheme
        
        for i = 1:nxC
            Dp[i,i] = 3/dx;
            if i-1 > 0
                Dp[i,i-1] = -4/dx;
            end
            if i-2 > 0
                Dp[i,i-2] = 1/2/dx;
            end
        end

        for i = 1:nxC
            Dm[i,i] = -3/dx;
            if i+1 < nxC+1
                Dm[i,i+1] = 4/dx;
            end
            if i+2 < nxC+1
                Dm[i,i+2] = -1/2/dx;
            end
        end

        for i = 1:nxC
            Dc[i,i] = -1/dx;
            Dc[i,i+1] = 1/dx;
        end

        for i = 1:nxC
            Dcx[i,i] = 1/dx;
            Dcx[i+1,i] = -1/dx;
        end

        new(x,xMid,settings,w,v,vp,vm,Dp,Dm,Dc,Dcx,settings.sigmaA,settings.sigmaS);
    end
 end

 function setupIC(obj::solver)
    g0 = zeros(obj.settings.NxC,obj.settings.Nv);
    rho0 = zeros(obj.settings.Nx);
    rho0 = ICrho(obj.settings,obj.x);
    g0 = ICg(obj.settings,obj.xMid);
    return rho0,g0;
 end


 function solveFullProblem(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    Nv = obj.settings.Nv;
    epsilon = obj.settings.epsilon;

    Dp = obj.Dp;
    Dm = obj.Dm;
    Dc = obj.Dc;
    Dcx = obj.Dcx;

    w = obj.w;
    v = obj.v;
    vp = obj.vp;
    vm = obj.vm;

    rho0,g0 = setupIC(obj);
    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    g1 = zeros(size(g0));
    rho1 = zeros(size(rho0));

    Nt = round(Tend/dt);
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    println("Running solver for the full problem")
    
    for k = 1:Nt
        # println(rho0)
        fac = epsilon^2/(epsilon^2 + obj.sigmaS*dt);
        g1 =  g0 + dt*(-(Dp * g0 * vp + Dm * g0 * vm)*(Iden - 1/2 * w * Transpose(unitvec))/epsilon - (Dc * rho0 * Transpose(unitvec) * v)/(epsilon^2));
        
        g1 =  fac * g1;
        
        rho1 = rho0 - dt/2 * Dcx * g1 * v * w;

        g0 = g1;
        rho0 = rho1;
        t = t + dt;
    end
    return t, rho1, g1;
 end