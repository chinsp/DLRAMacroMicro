__precompile__
include("quadrature.jl")

using ProgressMeter
using ProgressBars
using LinearAlgebra
using FastGaussQuadrature, LegendrePolynomials

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
        nPN = settings.Nv; # total number of Legendre polynomials used
        gamma = zeros(nPN+1); # vector with norms of the Legendre polynomials

        for i = 1:nPN+1
            gamma[i] = 2/(2*(i-1) + 1);
        end

        A = zeros(Float64,nPN,nPN); # reduced Flux matrix for the micro equation
        a_norm = zeros(Float64,nPN+1);

        for i = 1:nPN+1
            a_norm[i] = i/(sqrt((2i-1)*(2i+1)));
        end

        AFull = Tridiagonal(a_norm[1:end-1],zeros(nPN+1),a_norm[1:end-1]); # Full flux matrix

        A = AFull[2:end,2:end]; # extractign reduced flux matrix for the micro equations

        # TFull = zeros(nPN+1,nPN+1) # allocate transformation matrix
        # mu1, w1 = gausslegendre(nPN+1)
        # for k = 1:nPN+1
        #     P = collectPl(mu1[k],lmax = nPN);
        #     for i = 1:nPN+1
        #         TFull[i,k] = P[i-1]*sqrt(w1[k])/sqrt(gamma[i]);
        #     end
        # end
        # T = TFull[2:end,:];
        # #AbsA = T*abs.(diagm([0.0; mu[2:end]]))*T';
        # absA = T*abs.(diagm(mu1))*T';

        S = eigvals(Matrix(AFull));
        V = eigvecs(Matrix(AFull));
        absA_1 = V*abs.(diagm(S))*inv(V);
        absA = absA_1[2:end,2:end];
        

        # Mabs = broadcast(abs,M);
        # absA = R*Diagonal(Mabs)*Transpose(R); # Computing and setting Roe's matrix
        # absA = absA1[2:end,2:end];

        Abar = zeros(Float64,nPN);
        Abar[1] = sqrt(gamma[2]);
        # Abar = AFull[1,2:end];

        dx = settings.dx;
        
        # Stencil matrices for the Sn or Pn sovler based on whether the MM decompositon is used
        if settings.SolverType == 0 || settings.SolverType == 2
            Dp = zeros(Float64,nx,nx);
            Dm = zeros(Float64,nx,nx);
        elseif settings.SolverType == 1 || settings.SolverType == 3
            Dp = zeros(Float64,nxC,nxC);
            Dm = zeros(Float64,nxC,nxC);
        end

        # Stencil matrices for the Pn solver
        Dx = zeros(Float64,nxC,nxC);
        Dxx = zeros(Float64,nxC,nxC);

        # Stencil matrices for macro-micro interface
        Dc = zeros(Float64,nxC,nx);
        Dcx = zeros(Float64,nx,nxC);
        
        # Currently running a second order upwind scheme
        
        m = size(Dp)[1]

        if settings.SpatDisc == "FoUw"
            for i = 1:m
                Dp[i,i] = 1/dx;
                if i-1>0
                    Dp[i,i-1] = -1/dx;
                end
            end

            for i = 1:m
                Dm[i,i] = -1/dx;
                if i+1<m
                    Dm[i,i+1] = 1/dx;
                end
            end
        elseif settings.SpatDisc == "SoUw"
            for i = 1:m
                Dp[i,i] = 3/(2*dx);
                if i-1 > 0
                    Dp[i,i-1] = -4/(2*dx);
                end
                if i-2 > 0
                    Dp[i,i-2] = 1/(2*dx);
                end
            end
            
            for i = 1:m
                Dm[i,i] = -3/(2*dx);
                if i+1 < m+1
                    Dm[i,i+1] = 4/(2*dx);
                end
                if i+2 < m+1
                    Dm[i,i+2] = -1/(2*dx);
                end
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
        Dcx[1,1], Dcx[end,end] = 0,0;

        Dx = Tridiagonal(-ones(nxC-1)./2.0/dx, zeros(nxC), ones(nxC-1)./2.0/dx);
        Dxx = Tridiagonal(ones(nxC-1)./2.0/dx, -ones(nxC)./dx, ones(nxC-1)./2.0/dx);

        new(x,xMid,settings,w,v,vp,vm,AFull,A,absA,Abar,Dp,Dm,Dc,Dcx,Dx,Dxx,rho1,g1,settings.sigmaA,settings.sigmaS);
    end
 end

 function setupIC(obj::solver)
    g0 = zeros(obj.settings.NxC,obj.settings.Nv);
    rho0 = zeros(obj.settings.Nx);
    rho0 = ICrho(obj.settings,obj.x);
    g0 = ICg(obj.settings,obj.xMid,obj.v);
    return rho0,g0;
 end

# Sn solver for the kinetic equation witout macro-micr decomposition
function solveSN_kinetic(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    epsilon = obj.settings.epsilon;
    Nx = obj.settings.Nx
    Nv = obj.settings.Nv;
    epsilon = obj.settings.epsilon;
    epsilon = 1.0; # This solver is only for the kinetic equation thus we override the externally set epsilon

    Dp = obj.Dp;
    Dm = obj.Dm;
    Dc = obj.Dc;
    Dcx = obj.Dcx;

    w = obj.w;
    v = obj.v;
    vp = obj.vp;
    vm = obj.vm;

    rho0,g0 = setupIC(obj);
    g = zeros(Float64,Nx,Nv);
    for i in 1:Nv
        g[:,i] = rho0;
    end
    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    # g1 = obj.g1;
    # rho1 = obj.rho1;

    Nt = round(Tend/dt); # Compute the number of steps
    dt = Tend/Nt; # Find the step size 
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    for k = ProgressBar(1:Nt)
      g .= g .- dt.*Dp*g*vp .- dt.*Dm*g*vm .+ dt*obj.sigmaS.*g*(0.5*w*Transpose(unitvec) - Iden) - dt*obj.sigmaA.*g;

      t = t + dt;
    end
    return t, g;
end

#IMEX solver
 function solveFullProblem_Sn(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
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

    fac = 1 + dt*obj.settings.sigmaS/epsilon^2;

    println("Running solver for the Sn solver for the full problem")
    
    for k = ProgressBar(1:Nt)
        RHS = (Dc * rho0 * Transpose(unitvec) * v)/(epsilon^2);
        RHS =  RHS .+ (Dp * g0 * vp + Dm * g0 * vm)*(Iden - 0.5 * w * Transpose(unitvec))/epsilon;
        RHS = RHS .+ obj.sigmaA .* g0 .- Source_macro(obj.settings,t,obj.settings.xMid,obj.v)./epsilon;
        g1 =  g0 .- dt*RHS; 
        
        g1 =  g1./fac;
        
        rho1 = rho0 - dt *(0.5 * Dcx * g1 * v * w .- Source_micro(obj.settings,t,obj.settings.x)) ;

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
    rho1 = obj.rho1;

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

        rho1 .= rho0 + dt .* y;

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
    # Nx = obj.settings.Nx;
    # NxC = obj.settings.NxC;
    # Nv = obj.settings.Nv;
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

    fac = 1 + obj.settings.sigmaS*dt/epsilon^2;

    println("Running solver for the Pn solver for the full problem")

    for k = ProgressBar(1:Nt)
        g1 .= g0 + dt.*(-Dx*g0*Transpose(A)./epsilon + Dxx*g0*Transpose(absA)./epsilon - Dc*rho0*Transpose(Abar)./epsilon^2 - obj.settings.sigmaA*g0);
        g1 .= g1./fac;

        rho1 .= rho0 + dt.*(-0.5*Dcx*g1*Abar - obj.settings.sigmaA*rho0);

        rho0 .= rho1;
        g0 .= g1;
        
        # rho1[1] = rho1[end-1]; rho1[end] = rho1[2];
        # g1[1,:] = g1[end-1,:]; g1[end,:] = g1[2,:];

        t = t+dt;
    end
    return t,rho1,g1;

end


function solveMMDLRA_Sn(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    Nv = obj.settings.Nv;
    NxC = obj.settings.NxC;
    Nx = obj.settings.Nx;
    epsilon = obj.settings.epsilon;
    r = obj.settings.r;

    Dp = obj.Dp;
    Dm = obj.Dm;
    Dc = obj.Dc;
    Dcx = obj.Dcx;

    w = obj.w;
    v = obj.v;
    vp = obj.vp;
    vm = obj.vm;

    rho0,g0 = setupIC(obj);

    X,s,V = svd(g0);
    X = X[:,1:r];
    V = V[:,1:r];
    S = diagm(s[1:r]);

    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    g1 = obj.g1;
    rho1 = obj.rho1;

    Nt = round(Tend/dt); # Compute the number of steps
    dt = Tend/Nt; # Find the step size 
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    Sigma_S = obj.settings.sigmaS.*I(NxC);
    Sigma_A = obj.settings.sigmaA.*I(NxC);
    Sigma_AF = obj.settings.sigmaA.*I(Nx);

    fac = 1 + dt*obj.settings.sigmaS/epsilon^2;

    println("Running solver for the Sn solver for the full problem")

    for k = ProgressBar(1:Nt)
        # Solving the micro equation in time using DLRA
        K = X*S;
        K .= K .- dt.*(Dp*K*Transpose(V)*vp .+ Dm*K*Transpose(V)*vm)*(Iden - 0.5.*w*Transpose(unitvec))*V./epsilon .- dt.*Dc*rho0*Transpose(unitvec)*v*V./epsilon^2 .- dt.*Sigma_S*K./epsilon^2  .- dt.*Sigma_A*K;

        X1,R1 = qr(K);
        X1 = Matrix(X1)[:,1:r];
        M = Transpose(X1)*X;
        X .= X1;

        Lt = S*Transpose(V);
        Lt .= Lt .- dt.*Transpose(X)*(Dp*X*Lt*vp .+ Dm*X*Lt*vm)*(Iden - 0.5.*w*Transpose(unitvec))./epsilon .- dt.*Transpose(X)*Dc*rho0*Transpose(unitvec)*v./epsilon^2 .- dt.*Transpose(X)*Sigma_S*X*Lt./epsilon^2  .- dt.*Transpose(X)*Sigma_A*X*Lt;

        V1,R2 = qr(Transpose(Lt));
        V1 = Matrix(V1)[:,1:r];
        N = Transpose(V1)*V
        V .= V1;

        S .= M*S*Transpose(N);
        S .= S .- dt.*Transpose(X)*(Dp*X*Lt*vp .+ Dm*X*Lt*vm)*(Iden - 0.5.*w*Transpose(unitvec))*V./epsilon .- dt.*Transpose(X)*Dc*rho0*Transpose(unitvec)*v*V./epsilon^2 .- dt.*Transpose(X)*Sigma_S*X*S./epsilon^2  .- dt.*Transpose(X)*Sigma_A*X*S;

        # Solving the macro equation 

        rho1 .= rho0 .- 0.5*dt.*Dcx*X*S*Transpose(V)*v*w .- Sigma_AF*rho0;

        rho0 = rho1;
        t = t + dt;
    end
    return t, rho1, X*S*Transpose(V);
end

function solveMMDLRA_SnIMEX(obj::solver)
    t = 0.0;
    dt = obj.settings.dt;
    Tend = obj.settings.Tend;
    Nv = obj.settings.Nv;
    NxC = obj.settings.NxC;
    Nx = obj.settings.Nx;
    epsilon = obj.settings.epsilon;
    r = obj.settings.r;

    Dp = obj.Dp;
    Dm = obj.Dm;
    Dc = obj.Dc;
    Dcx = obj.Dcx;

    w = obj.w;
    v = obj.v;
    vp = obj.vp;
    vm = obj.vm;

    rho0,g0 = setupIC(obj);

    X,s,V = svd(g0);
    X = X[:,1:r];
    V = V[:,1:r];
    S = diagm(s[1:r]);

    # println(rho0)
    ## pre=allocating memory for solution of macro and micro equation
    g1 = obj.g1;
    rho1 = obj.rho1;

    Nt = round(Tend/dt); # Compute the number of steps
    dt = Tend/Nt; # Find the step size 
    
    unitvec = ones(Nv);
    Iden = I(Nv);

    Sigma_S = obj.settings.sigmaS.*I(NxC);
    Sigma_A = obj.settings.sigmaA.*I(NxC);
    Sigma_AF = obj.settings.sigmaA.*I(Nx);

    fac = 1 + dt*obj.settings.sigmaS/epsilon^2;

    FacK = inv(I(NxC) + dt.*Sigma_S./epsilon^2);

    println("Running solver for the Sn solver for the full problem")

    for k = ProgressBar(1:Nt)
        # Solving the micro equation in time using DLRA
        K = X*S;
        K .= K .- dt.*(Dp*K*Transpose(V)*vp .+ Dm*K*Transpose(V)*vm)*(Iden - 0.5.*w*Transpose(unitvec))*V./epsilon .- dt.*Dc*rho0*Transpose(unitvec)*v*V./epsilon^2   .- dt.*Sigma_A*K;
        K .= FacK*K;
        X1,R1 = qr(K);
        X1 = Matrix(X1)[:,1:r];
        M = Transpose(X1)*X;
        X .= X1;

        Lt = S*Transpose(V);
        FacL = inv(I(r) + dt.*Transpose(X)*Sigma_S*X./epsilon^2);
        Lt .= Lt .- dt.*Transpose(X)*(Dp*X*Lt*vp .+ Dm*X*Lt*vm)*(Iden - 0.5.*w*Transpose(unitvec))./epsilon .- dt.*Transpose(X)*Dc*rho0*Transpose(unitvec)*v./epsilon^2 .- dt.*Transpose(X)*Sigma_A*X*Lt;
        Lt .= FacL*Lt;
        V1,R2 = qr(Transpose(Lt));
        V1 = Matrix(V1)[:,1:r];
        N = Transpose(V1)*V
        V .= V1;

        S .= M*S*Transpose(N);
        FacS = inv(I(r) + dt.*Transpose(X)*Sigma_S*X./epsilon^2);
        S .= S .- dt.*Transpose(X)*(Dp*X*Lt*vp .+ Dm*X*Lt*vm)*(Iden - 0.5.*w*Transpose(unitvec))*V./epsilon .- dt.*Transpose(X)*Dc*rho0*Transpose(unitvec)*v*V./epsilon^2   .- dt.*Transpose(X)*Sigma_A*X*S;
        S .= FacS*S;
        # Solving the macro equation 

        rho1 .= rho0 .- 0.5*dt.*Dcx*X*S*Transpose(V)*v*w .- Sigma_AF*rho0;

        rho0 = rho1;
        t = t + dt;
    end
    return t, rho1, X*S*Transpose(V);
end