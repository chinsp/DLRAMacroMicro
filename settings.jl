__precompile__

mutable struct Settings
    ## Settings of the staggered grids
    # Number of spatial vertices
    Nx::Int64;
    # Number of cell centres
    NxC::Int64;
    # Start and end point of spatial domain
    a::Float64;
    b::Float64;
    # Grid cell width
    dx::Float64;
    
    ## Settings of temporal domain
    # End time
    Tend::Float64;
    # Time increment width
    dt::Float64;
    # CFL number 
    cfl1::Float64; # For parabolic CFL condition
    cfl2::Float64; # For hyperbolic CFL condition
    # CFL condition type
    cflType::String;

    ## Settings for angular approximation
    # Number of quadrature points
    Nv::Int64;
    
    ## Spatial grid
    x
    xMid

    ## Problem Settings
    problem::String;

    ## Rescaled collision length
    epsilon::Float64    

    ## Initial conditions
    ICType::String;
    BCType::String;

    ## Physical parameters
    sigmaS::Float64; ## Try to change this to get a non-constant value for the scattering coefficients
    sigmaA::Float64;
    
    ## Low-rank approximation parameters
    r::Int; # rank of approximation

    ## Parameters for setting stencil matrices for spatial derivaties
    SolverType::Int64;

    SpatDisc::String;

    function Settings(Nx::Int=1001,Nv::Int=500,epsilon::Float64=1.0,cflType::String="hyperbolic",SolverType::Int=3)
        # Setup spatial grid
        NxC = Nx - 1;
        a = -1; # Starting point for the spatial interval
        b = 1; # End point for the spatial interval

        # Setup temporal discretisation
        Tend = 5;
        cfl1 = 1.0; # CFL condition parabolic
        cfl2 = 0.5; # CFL condition hyperbolic
        # cflType = "parabolic"; # or "parabolic", "mixed"

        # epsilon = 10^-6;

        # Initial conditions
        ICType = "MS" ; # Options "LS" or "MS"

        # Problem 
        problem = "ManufacturedSolution"; #Options are "LineSource" or "ManufacturedSolution"

        x = collect(range(a,stop = b,length = Nx));
        dx = x[2] - x[1];
        xMid = x .+ dx/2;
        xMid = xMid[1:(end-1)];

        println("Number of points for spatial discretisation of macro = ",Nx);
        println("Number of points for spatial discretisation of micro = ",NxC);

        # Setting up the time increment
        if cflType == "parabolic"
            dt = cfl1*(dx^2);
        elseif cflType == "hyperbolic"
            dt = cfl2*dx;
        elseif cflType == "mixed"
            dt = cfl1*(dx^2) + cfl2*epsilon*dx;
        else
            println("Please enter valid type for CFL condition")
        end

        # Physical parameters
        if problem == "LineSource"
            sigmaS = 1.0;
            sigmaA = 0.0;
        elseif problem == "ManufacturedSolution"
            sigmaS = 1.0;
            sigmaA = 0.0;
        end

        BCType = "exact"

        r = 30;

        ## A string variable to describe the type of discretisation used for the spatial derivative
        # FoUw - First-order Upwind scheme
        # SoUw - Second-order Upwind scheme
        SpatDisc = "FoUw";

        new(Nx,NxC,a,b,dx,Tend,dt,cfl1,cfl2,cflType,Nv,x,xMid,problem,epsilon,ICType,BCType,sigmaS,sigmaA,r,SolverType,SpatDisc);
    end
end

function ICrho(obj::Settings,x)
    y = zeros(size(x));
    if obj.ICType == "LS"
        s1 = 0.03;
        s2 = s1^2;
        floor = 1e-4;
        x0 = 0.0;
        for j in eachindex(y)
            y[j] = max(floor,1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2));
        end
    elseif obj.ICType == "MS"
        for j in eachindex(y)
            y[j] = sin(pi*x[j]);
        end
    end
    return y;
end

function ICg(obj::Settings,x,v)
    y = zeros(length(x),obj.Nv);
    m = length(x);
    if obj.ICType == "LS"
        y = y;
    elseif obj.ICType == "MS"
        for i = 1:m
            for j = 1:obj.Nv
                y[i,j] = v[j,j]*sin(pi*x[i]);
            end
        end
    end
    return y;
end

function Source_macro(obj::Settings,t,x,v)
    m = length(x);
    y = zeros(m,obj.Nv);
    if obj.ICType == "MS"    
        for i = 1:m
            for j = 1:obj.Nv
                y[i,j] = exp(-t)*sin(pi*x[i])*v[j,j]*(1/obj.epsilon - obj.epsilon) + pi*exp(-t)*cos(pi*x[i])*(v[j,j]/obj.epsilon + v[j,j]^2 - 1/3);
            end
        end
    end
    return y;
end

function Source_micro(obj::Settings,t,x)
    m = length(x);
    y = zeros(m);
    if obj.ICType == "MS"
        for i = 1:m
            y[i] = exp(-t)*(pi*cos(pi*x[i])/3 - sin(pi*x[i]));
        end
    end
    return y;
end

function Manufactured1D_rho(obj::Settings,t,x)
    y = exp(-t)*sin.(pi*x);
    return y;
end