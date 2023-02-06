__precompile__
include("quadrature.jl")

using ProgressMeter
using ProgressBars
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

    rho1::Array{Float64,1};
    g1::Array{Float64,2};
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

        g1 = zeros(nxC,Nv);
        rho1 = zeros(nx);

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
            Dp[i,i] = 3/(2*dx);
            if i-1 > 0
                Dp[i,i-1] = -4/(2*dx);
            end
            if i-2 > 0
                Dp[i,i-2] = 1/(2*dx);
            end
        end

        for i = 1:nxC
            Dm[i,i] = -3/(2*dx);
            if i+1 < nxC+1
                Dm[i,i+1] = 4/(2*dx);
            end
            if i+2 < nxC+1
                Dm[i,i+2] = -1/(2*dx);
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

        new(x,xMid,settings,w,v,vp,vm,Dp,Dm,Dc,Dcx,rho1,g1,settings.sigmaA,settings.sigmaS);
    end
 end

 function setupIC(obj::solver)
    g0 = zeros(obj.settings.NxC,obj.settings.Nv);
    rho0 = zeros(obj.settings.Nx);
    rho0 .= ICrho(obj.settings,obj.x);
    g0 .= ICg(obj.settings,obj.xMid);
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
    g1 = obj.g1;
    rho1 = obj.rho1;

    Nt = round(Tend/dt);
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    println("Running solver for the full problem")
    
    for k = ProgressBar(1:Nt)
        RHS = zeros(NxC,Nv);
        fac = 1 + dt*obj.sigmaS/epsilon^2;
        RHS .= (Dc * rho0 * Transpose(unitvec) * v)./(epsilon^2);
        RHS .=  RHS .+ (Dp * g0 * vp .+ Dm * g0 * vm)*(Iden .- 0.5 * w * Transpose(unitvec))./epsilon;
        RHS .= RHS .+ obj.sigmaA .* g0 ;
        g1 .=  (g0 .- dt*RHS)./fac; 
        
        rho1 = rho0 - dt *(0.5 .* Dcx * g1 * v * w) ;

        rho1[1],rho1[end] = 0,0;

        g0 = g1;
        rho0 = rho1;
        t = t + dt;
    end
    return t, rho1, g1;
 end

 function solveLimitingDiff(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    dx = obj.settings.dx;

    rho0,g0 = setupIC(obj);
    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    g1 = obj.g1;

    Nt = round(Tend/dt);

    println("Running solver for the limiting diffusion equation")

    for k = ProgressBar(1:Nt)
        y = zeros(Nx);
        for i = 1:Nx
            if i == 1
                y[i] = (rho0[i+1] - 2*rho0[i])./ (dx^2 * 3 * obj.sigmaS);
            elseif i == Nx
                y[i] = (-2*rho0[i] + rho0[i-1]) ./ (dx^2 * 3 * obj.sigmaS)   ;
            else
                y[i] = (rho0[i+1] - 2*rho0[i] + rho0[i-1])./ (dx^2 * 3 * obj.sigmaS);
            end
        end

        rho1 = rho0 + dt .* y;

        rho0 = rho1;
        t = t + dt;
    end

    return t, rho0;
 end