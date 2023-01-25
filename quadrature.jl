__precompile__
# import FastTransforms
using FastGaussQuadrature

struct Quadrature
    Nv::Int;
    w = Array{Float64,1}
    v = Array{Float64,1};

    function Quadrature(Nv,quadtype)
        if quadtype == "Gauss"
            v,w = gausslegendre(Nv);
        else
            println("Entered quadrature not available");
        end
        v = v[end:-1:1];
        w = w[end:-1:1];

        new(Nv,collect(w),collect(v));

    end
 end

