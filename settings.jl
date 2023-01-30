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

    function Settings(Nx::Int=1001,Nv::Int=500,epsilon::Float64=1.0,cflType::String="hyperbolic")
        # Setup spatial grid
        NxC = Nx - 1;
        a = -1.5; # Starting point for the spatial interval
        b = 1.5; # End point for the spatial interval

        # Setup temporal discretisation
        Tend = 5;
        cfl1 = 0.1; # CFL condition parabolic
        cfl2 = 0.1; # CFL condition hyperbolic
        # cflType = "parabolic"; # or "parabolic", "mixed"

        # epsilon = 10^-6;

        # Initial conditions
        ICType = "LS" ;

        # Problem 
        problem = "LineSource";

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
        end

        BCType = "exact"

        r = 30;

        new(Nx,NxC,a,b,dx,Tend,dt,cfl1,cfl2,cflType,Nv,x,xMid,problem,epsilon,ICType,BCType,sigmaS,sigmaA,r);
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
            y[j] = max(floor,1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2))
        end
    elseif obj.ICType == "ManufacturedSolution"
        println("Not coded yet")
    end
    return y;
end

function ICg(obj::Settings,x)
    y = zeros(length(x),obj.Nv);
    if obj.ICType == "LS"
        y = y;
    elseif obj.ICType == "ManufacturedSolution"
        println("Not coded yet");
    end
    return y;
end