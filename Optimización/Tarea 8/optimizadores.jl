#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 8

using LinearAlgebra

function descensoNewton(f::Function,x0::AbstractVector,∇::Function, H::Function)
    #Newton's method with backtracking step length search
    maxIter = 1000;
    tol = 1e-10;
    k=0;
    ρ = 0.9;
    x = copy(x0);
    f0 = f(x);
    f1 = f0;
    g = ∇(x); #Initialize the gradient ∇f(x_k)
    n = length(g);
    B = H(x); #Initialize the Hessian Matrix ∇^2f(x_k)
    Ng = norm(g); #Initialize the norm of the gradient ||∇f(x_k)||
    d = copy(-g) #Initialize the descent direction
    maxjter = 20;   # Backtracking iterations
    c1 = 0.1;   # First coeficient for Wolfe conditions
    xt = copy(x)    # Initialize trial step
    while k<maxIter && Ng>tol #Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        B .= H(x); #Calculate the Hessian Matrix 
        # Modify B so it's positive definite:
        r = rank(B);
        while r < n
            # If B is singular
            B .= 0.9B + 0.1I
            r = rank(B);
        end
        S = bunchkaufman(Symmetric(collect(B))); # B=P'UDU'P factorization (L'DL) with permutations
        minD = minimum(diag(S.D));
        if minD < tol
            B = S.P'*S.U*abs.(S.D)*S.U'*S.P; #B is P.D. iff D has positive entries
        end
        d .= -B\g; #Calculate Newton's point
        α  = 1;     # Full step
        # Line search
        xt .= x + α*d;  # Calculate the next point
        f1 = f(xt); # Function value at the updated point
        s = d'*g;   # Directional derivative
        jter = 0;   # Backtracking iterations
        while f1 > f0 + α*c1*s && jter < maxjter # Backtracking steps
            α = α*ρ;    # Backtrack α by ρ
            xt .= x + α*d;  # Update the next point
            f1 = f(xt); # Function value at the updated point
            jter += 1;
        end
        x .= xt; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g .= ∇(x); #Update the gradient ∇f(x_k)
        Ng = norm(g); #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
    end
    return x, Ng, k
end
