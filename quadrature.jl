__precompile__
# import FastTransforms
using FastGaussQuadrature

struct Quadrature
    Nv::Int;
    w
    v

    function Quadrature(Nv,quadtype)
        if quadtype == "Gauss"
            v,w = gausslegendre(Nv);
        else
            println("Entered quadrature not available");
        end
        v = v[end:-1:1]; # Changes the direction of the domain
        w = w[end:-1:1]; # Adjusts the weights according to the quadrature points

        new(Nv,collect(w),collect(v));

    end
 end

