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

    ## Angular discretisation

    # Sn discretisation
    w::Array{Float64,1};
    v::Array{Float64,2};
    vp::Array{Float64,2};
    vm::Array{Float64,2};

    # Pn discretisation 
    AFull::Array{Float64,2};
    A::Array{Float64,2};
    absA::Array{Float64,2};
    Abar::Array{Float64,1};

    # Stencil matrices for spatial discretisation
    Dp::Array{Float64,2};
    Dm::Array{Float64,2};
    Dc::Array{Float64,2};
    Dcx::Array{Float64,2};
    Dx::Array{Float64,2};
    Dxx::Array{Float64,2};

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

        # Setting up matrices for Sn solver

        quad = Quadrature(Nv,"Gauss");

        w = quad.w[end:-1:1];
        v = Diagonal(quad.v[end:-1:1]);
        
        vp = Diagonal(zeros(Float64,Nv));
        vm = Diagonal(zeros(Float64,Nv));
        for q = 1:Nv
            vp[q,q] = max(v[q,q],0);
            vm[q,q] = min(v[q,q],0);
        end

        # Setting up the matrices for the Pn solver
        nPN = settings.N + 1; # total number of Legendre polynomials used
        gamma = zeros(nPN); # vector with norms of the Legendre polynomials

        for i = 1:nPN
            gamma[i] = 2/(2*(i-1) + 1);
        end

        A = zeros(Float64,nPN-1,nPN-1); # reduced Flux matrix for the micro equation
        a_norm = zeros(Float64,nPN);

        for i = 1:nPN
            a_norm[i] = i/(sqrt((2i-1)*(2i+1)));
        end

        AFull = Tridiagonal(a_norm[1:end-1],zeros(nPN),a_norm[1:end-1]); # Full flux matrix

        A = AFull[2:end,2:end]; # extractign reduced flux matrix for the micro equations

        M,R = eigen(AFull);

        Mabs = broadcast(abs,M);
        absA1 = R*Diagonal(Mabs)*Transpose(R); # Computing and setting Roe's matrix
        absA = absA1[2:end,2:end];

        Abar = zeros(Float64,nPN-1);
        Abar[1] = gamma[2];

        dx = settings.dx;
        
        # Stencil matrices for the Sn sovler
        Dp = zeros(Float64,nxC,nxC);
        Dm = zeros(Float64,nxC,nxC);

        # Stencil matrices for the Pn solver
        Dx = zeros(Float64,nxC,nxC);
        Dxx = zeros(Float64,nxC,nxC);

        # Stencil matrices for macro-micro interface
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

        Dx = Tridiagonal(-ones(nxC-1)./2/dx, zeros(nxC), ones(nxC-1)./2/dx);
        Dxx = Tridiagonal(ones(nxC-1)./2/dx, -ones(nxC)./dx, ones(nxC-1)./2/dx);

        new(x,xMid,settings,w,v,vp,vm,AFull,A,absA,Abar,Dp,Dm,Dc,Dcx,Dx,Dxx,rho1,g1,settings.sigmaA,settings.sigmaS);
    end
 end

 function setupIC(obj::solver)
    g0 = zeros(obj.settings.NxC,obj.settings.Nv);
    rho0 = zeros(obj.settings.Nx);
    rho0 = ICrho(obj.settings,obj.x);
    g0 = ICg(obj.settings,obj.xMid);
    return rho0,g0;
 end

#IMEX solver
 function solveFullProblem_Sn(obj::solver)
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

    Nt = round(Tend/dt); # Compute the number of steps
    dt = Tend/Nt; # Find the step size 
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    println("Running solver for the Sn solver for the full problem")
    
    for k = ProgressBar(1:Nt)
        fac = epsilon^2/(epsilon^2 + obj.sigmaS*dt);
        RHS = (Dc * rho0 * Transpose(unitvec) * v)/(epsilon^2);
        RHS =  RHS .+ (Dp * g0 * vp + Dm * g0 * vm)*(Iden - 0.5 * w * Transpose(unitvec))/epsilon;
        RHS = RHS .+ obj.sigmaA .* g0;
        g1 =  g0 .- dt*RHS; 
        
        g1 =  fac * g1;
        
        rho1 = rho0 - dt *(0.5 * Dcx * g1 * v * w) ;

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

    return t, rho1;
 end


#IMEX solver
 function solveFullProblem_Pn(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    Nv = obj.settings.Nv;
    epsilon = obj.settings.epsilon;

    A = obj.A;
    absA = obj.absA;
    Abar = obj.Abar;

    Dx = obj.Dx;
    Dxx = obj.Dxx;
    Dc = obj.Dc;
    Dcx = obj.Dcx;

    rho0,g0 = setupIC(obj);
    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    g1 = obj.g1;
    rho1 = obj.rho1;

    Nt = round(Tend/dt); # Computing the number of steps required 
    dt = Tend/Nt; # Adjusting the step size 

    fac = 1 + obj.settings.sigmaS*dt/epsilon^;

    println("Running solver for the Pn solver for the full problem")

    for k =ProgressBar(1:Nt)
        g1 = g0 + dt.*(-Dx*g0*Transpose(A)./epsilon + Dxx*g0*Transpose(A)./epsilon - Dc*rho0*Transpose(Abar) - obj.settings.sigmaA*g0);
        g1 .= g1./fac;

        rho1 = rho0 + dt.*(-0.5*Abar[1]*Dcx*g1 - obj.settings.sigmaA*rho0);

        rho0 .= rho1;
        g0 .= g1;

        t = t+dt;
    end
    return t,rho1,g1;

end
