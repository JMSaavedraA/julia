#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 11


using NPZ, Printf

P = npzread("puntos2D_ej3.npy");

x = P[:,1];
y = P[:,2];

function residual(z::AbstractVector)
    A = z[1];
    ω = z[2];
    ϕ = z[3];
    t = ω*x .+ ϕ;
    R = A*sin.(t) - y;
    return R
end

function jacobianaResidual(z::AbstractVector)
    A = z[1];
    ω = z[2];
    ϕ = z[3];
    t = ω*x .+ ϕ;
    j1 = sin.(t);
    j3 = A*cos.(t);
    j2 = x .* j3;
    J = [j1 j2 j3];
    return J
end


function showFitResults(z, f, k, nP, sol)
    if sol
        println("El algoritmo ha convergido")
    else
        println("El algoritmo no convergió")
    end
    
    println("El valor encontrado es:")
    println(@sprintf("A = %2.2f, ω = %2.2f y ϕ = %2.2f",z[1],z[2],z[3]))
    println(@sprintf("El valor de f(z)=%g",f))
    println(@sprintf("El valor de ||p||=%g",nP))
    println(@sprintf("Se realizaron k=%i iteraciones",k))
end
