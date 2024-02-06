#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization

using LinearAlgebra

function pasoFijo(f::Function,g::Function, x::AbstractVector, maxIter::Int,tol::Number,α::Number)
    #Fixed step size descent directions algorithm with steepest descent
    F = zeros(maxIter+1); #We store f(x_k) to plot it later
    G = zeros(maxIter+1); #We store ∇f(x_k) to plot it later
    k = 0;
    gk = g(x); #Initialize the gradient ∇f(x_k)
    n = norm(gk); #Initialize the norm of the gradient ||∇f(x_k)||
    F[1] = f(x);
    G[1] = n;
    while k<maxIter && n>tol #Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        x -= α*gk; #Update x_k+1 as x_k - α ∇f(x_k)
        k = k+1;
        gk .= g(x); #Update the gradient ∇f(x_k)
        n = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        F[k+1] = f(x);
        G[k+1] = n;
    end
    F = F[1:k+1];
    G = G[1:k+1];
    return x, n, F, G
end

function descensoNewton(f::Function,g::Function, H::Function, x::AbstractVector, maxIter::Int,tol::Number)
    #Descent directions algorithm with Newton's point. α=1
    F=zeros(maxIter+1); #We store f(x_k) to plot it later
    G=zeros(maxIter+1); #We store ∇f(x_k) to plot it later
    k=0;
    gk = g(x); #Initialize the gradient ∇f(x_k)
    B = H(x); #Initialize the Hessian Matrix ∇^2f(x_k)
    gz = norm(gk); #Initialize the norm of the gradient ||∇f(x_k)||
    d = copy(-gk) #Initialize the descent direction
    F[1] = f(x);
    G[1] = gz;
    while k<maxIter && gz>tol #Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        B.=H(x); #Calculate the Hessian Matrix 
        d.= B\gk; #Calculate Newton's point
        x-=d; #Update x_k+1 = x_k + d_k
        k = k+1;
        gk .= g(x); #Update the gradient ∇f(x_k)
        gz = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        F[k+1] = f(x);
        G[k+1] = gz;
    end
    F = F[1:k+1];
    G = G[1:k+1];
    return x, gz, F, G
end
