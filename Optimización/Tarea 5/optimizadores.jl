#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization

using LinearAlgebra

function backtrackingDescent(f::Function,∇::Function, x0::AbstractVector, maxIter::Int,tol::Number,ρ::Number)
    # Fixed step size descent directions algorithm with steepest descent
    x = copy(x0)
    F = zeros(maxIter+1); # We store f(x_k) to plot it later
    G = zeros(maxIter+1); # We store ∇f(x_k) to plot it later
    k = 0;
    g = ∇(x); # Initialize the gradient ∇f(x_k)
    n = norm(g); # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = f(x);
    f1 = f0;
    F[1] = f0;
    G[1] = n;
    maxjter = 12;   # Backtracking iterations
    c1 = 0.1;   # First coeficient for Wolfe conditions
    xt = copy(x)    # Initialize trial step
    d = copy(g) # Initialize the descent direction
    xChange = n
    while k < maxIter && n > tol && xChange > tol^3 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g     # Steepest descent direction
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
        n = norm(g);    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        F[k+1] = f0;    
        G[k+1] = n;
        xChange = α*norm(d)/maximum([norm(x),1])    # Relative and absolute change for x
    end
    F = F[1:k+1];
    G = G[1:k+1];
    return x, n, F, G
end



function bisectionDescent(f::Function,∇::Function, x0::AbstractVector, maxIter::Int,tol::Number)
    # Fixed step size descent directions algorithm with steepest descent
    x = copy(x0)
    F = zeros(maxIter+1); # We store f(x_k) to plot it later
    G = zeros(maxIter+1); # We store ∇f(x_k) to plot it later
    k = 0;
    g0 = ∇(x); # Initialize the gradient ∇f(x_k)
    n = norm(g0); # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = f(x);
    f1 = f0;
    F[1] = f0;
    G[1] = n;
    maxjter = 30;   # Line search iterations
    c1 = 0.1;   # First coeficient for Wolfe conditions
    c2 = 0.9    # Second coeficient for Wolfe conditions
    xt = copy(x)    # Initialize trial step
    g1 = ∇(xt)  # Gradient at trial step
    d = copy(g0)    # Initialize the descent direction
    xChange = n
    while k < maxIter && n > tol && xChange > tol^2 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g0     # Steepest descent direction
        # Line search
        δ  = 0.0; # Left side of interval
        γ = Inf # Right side of interval
        α = 0.1   # Initial guess
        xt .= x + α*d;  # Calculate the next point
        f1 = f(xt); # Function value at the updated point
        s = d'*g0;   # Directional derivative
        jter = 0;   # Bisection iterations
        goOn = true
        g1 .= ∇(xt)
        while goOn && jter < maxjter
            if f1 > f0 + α*c1*s # Armijo's condition
                γ = α # Update right side
                α = (δ + γ)/2 # Update α
                xt .= x + α*d;  # Calculate the next point
                f1 = f(xt); # Function value at the updated point
                g1 .= ∇(xt) # Gradient at the next point
            elseif d'*g1 < c2*s # Wolfe's condition
                δ = α #Update left side
                # Update α
                if isinf(γ)
                    α = 2*δ
                else
                    α = (δ + γ)/2
                end
                xt .= x + α*d;  # Calculate the next point
                f1 = f(xt); # Function value at the updated point
                g1 .= ∇(xt) # Update the gradient
            else # Stop if both condition are accomplished
                goOn = false
            end
            jter += 1
        end
        x .= xt; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g0 .= g1; #Update the gradient ∇f(x_k)
        n = norm(g0);    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        F[k+1] = f0;    
        G[k+1] = n;
        xChange = α*norm(d)/maximum([norm(x),1])    # Relative and absolute change for x
    end
    F = F[1:k+1];
    G = G[1:k+1];
    return x, n, F, G
end


function zoomDescent(f::Function,∇::Function, x0::AbstractVector, maxIter::Int,tol::Number)
    # Fixed step size descent directions algorithm with steepest descent
    x = copy(x0)
    F = zeros(maxIter+1); # We store f(x_k) to plot it later
    G = zeros(maxIter+1); # We store ∇f(x_k) to plot it later
    k = 0;
    g0 = ∇(x); # Initialize the gradient ∇f(x_k)
    n = norm(g0); # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = f(x);
    f1 = f0;
    F[1] = f0;
    G[1] = n;
    maxjter = 30;   # Line search iterations
    c1 = 1e-4;  # First coeficient for Wolfe conditions
    c2 = 0.9    # Second coeficient for Wolfe conditions
    xt = copy(x)    # Initialize trial step
    aMax = 2    # For zoom algorithm
    g1 = ∇(xt)  # Gradient at trial step
    d = copy(g0)    # Initialize the descent direction
    while k < maxIter && n > tol # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -(g0)     # Steepest descent direction
        # Line search
        s = d'*g0;   # Directional derivative
        jter = 0;   # Line Search iterations
        goOn = true
        a0 = 0.0    # α_0 = 0
        a1 = 1.0    # Initialize α_1
        α = 0.01    # Initialize α
        p0 = f0
        while goOn && jter < maxjter
            xt .= x + a1*d
            f1 = f(xt)
            if f1 > minimum([f0 + a1*c1*s,p0]) # Armijo's condition
                α, f1, g1, s1 = zoom(a0,a1,x,d,f,f0,s,c1,c2,∇)
                goOn = false
            else
                g1 .= ∇(xt)
                s1 = d'*g1
                if abs(s1) <= -c2*s #Second strong Wolfe condition
                    α = a1
                    goOn = false
                elseif s1 >= 0
                    α, f1, g1, s1 = zoom(a1,aMax,x,d,f,f0,s,c1,c2,∇)
                    goOn = false
                else
                    a1 = (a1 + aMax)/2
                end
            end
            jter += 1
        end
        g1 .= g(xt) # If for some reason the gradient was not updated
        x += α*d; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g0 .= g1; #Update the gradient ∇f(x_k)
        n = norm(g0);    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        F[k+1] = f0;    
        G[k+1] = n;
    end
    F = F[1:k+1];
    G = G[1:k+1];
    return x, n, F, G
end

function zoom(a::Number,b::Number,x::AbstractVector,d::AbstractVector,f::Function,f0::Number,s::Number,c1::Number,c2::Number,∇::Function)
    looking = true
    kter = 0
    maxkter = 20    # Maximum number of "Zooms" in the interval
    α = (a+b)/2 # Initialize α
    while looking && kter < maxkter
        c = (a+b)/2 # Midpoint of the interval
        fa = f(x + a*d) # Value at the left side of interval
        fb = f(x + b*d) # Value at the right side of the interval
        fc = f(x + c*d) # Value at the middle of the interval
        if (fa+fb-2*fc) > 1e-6  # Linearity if fa + fb - 2fc = 0
            α = (fa*(a+3*b)+fb*(3*a+b)-4*fc*(a+b))/(4*(fa+fb-2*fc)) # α_min for quadratic interpolation on a, b, (a+b)/2
        else
            α = 0.8*b + 0.2*a   # a + .8(b-a)
        end
        xt = x + α*d    # Trial step
        f1 = f(xt)
        if f1 > f0 + α*c1*s # Armijo's condition
            b = α   # Reduce right side of interval
        else
            g1 = ∇(xt)  # Gradient at the trial step
            s1 = d'*g1  # Φ'(α)
            if abs(s1) <= -c2*s # Second strong Wolfe condition
                looking = false
            elseif (b-a)*s1 >= 0
                b = a
                a = α
            else
                a = α
            end
        end
        kter += 1
    end
    xt = x + α*d
    f1 = f(xt)
    g1 = ∇(xt)
    s1 = d'*g1
    return α, f1, g1, s1
end