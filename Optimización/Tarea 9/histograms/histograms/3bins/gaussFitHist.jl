using DelimitedFiles, LinearAlgebra

function read3dHistogram(fileName::String)
    # Auxiliary function that reads a file containing a segmented 3d histogram.
    N = readdlm(fileName)[1,:];
    H = reshape(readdlm(fileName,'\n')[2:end],N[1],N[2],N[3]);
    return H, N
end



function objectiveGaussianFit(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # Gaussian fit objective function for 3d histogram H
    s = 0.0;
    for i = 1:N[1]
        for j = 1:N[2]
            for k = 1:N[3]
                ss = 0.0;
                for m = 1:n
                    ss += α[m]*exp(- norm([i,j,k] - μ[:,m]) /(2*σ^2));
                end
                s += (H[i,j,k] - ss)^2;
            end
        end
    end
    return s
end

function objectiveGaussianFitFDaGrad(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # ∇ f(α,μ) with respect to α approximation via central Finite Differences
    h = 1e-8 * norm(α)/n;
    ∇ = zeros(n);
    αL = copy(α);
    αR = copy(α);
    h2 = 2*h;
    αL[1] -= h; # Forward step
    αR[1] += h; # Backward step
    ∇[1] = (objectiveGaussianFit(αR,μ,N,n,H,σ) - objectiveGaussianFit(αL,μ,N,n,H,σ))/(h2); # Finite difference
    for i=2:n
        αL[i-1] = α[i-1]; # Restore to original value
        αR[i-1] = α[i-1]; # Restore to original value
        αL[i] -= h; # Forward step
        αR[i] += h; # Backward step
        ∇[i] = (objectiveGaussianFit(αR,μ,N,n,H,σ) - objectiveGaussianFit(αL,μ,N,n,H,σ))/(h2); # Finite difference
    end
    return ∇
end

function objectiveGaussianFitFDMuGrad(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # ∇ f(α,μ) with respect to μ approximation via central Finite Differences
    h = 1e-8 * norm(μ)/(n^2);
    ∇ = zeros(n,n);
    mL = copy(μ);
    mR = copy(μ);
    h2 = 2*h;
    for i=1:n
        for j=1:n
            mL[i,j] = μ[i,j] - h; # Forward step
            mR[i,j] = μ[i,j] + h; # Backward step
            ∇[i,j] = (objectiveGaussianFit(α,mR,N,n,H,σ) - objectiveGaussianFit(α,mL,N,n,H,σ))/(h2); # Finite difference
            mL[i,j] = μ[i,j]; # Restore to original value
            mR[i,j] = μ[i,j]; # Restore to original value
        end
    end
    return ∇
end

function backtrackingDescentAlpha(a0::AbstractVector,μ::AbstractMatrix, maxIter::Int,tol::Number,ρ::Number,σ::Number,H::AbstractArray,N::AbstractVector,u::Int)
    # Backtracking descent directions algorithm with steepest descent for the α subproblem on gaussian fit problem
    a = copy(a0);
    k = 0;
    g = objectiveGaussianFitFDaGrad(a,μ,N,u,H,σ); # Initialize the gradient ∇f(x_k)
    n = norm(g); # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = objectiveGaussianFit(a,μ,N,u,H,σ);
    f1 = f0;
    maxjter = 12;   # Backtracking iterations
    c1 = 1e-3;   # First coeficient for Wolfe conditions
    at = copy(a)    # Initialize trial step
    d = copy(g) # Initialize the descent direction
    xChange = n
    while k < maxIter && n > tol && xChange > tol^3 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g;     # Steepest descent direction
        α  = 1;     # Full step
        # Line search
        at .= a + α*d;  # Calculate the next point
        f1 = objectiveGaussianFit(at,μ,N,u,H,σ); # Function value at the updated point
        s = d'*g;   # Directional derivative
        jter = 0;   # Backtracking iterations
        while f1 > f0 + α*c1*s && jter < maxjter # Backtracking steps
            α = α*ρ;    # Backtrack α by ρ
            at .= a + α*d;  # Update the next point
            f1 = objectiveGaussianFit(at,μ,N,u,H,σ); # Function value at the updated point
            jter += 1;
        end
        a .= at; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g .= objectiveGaussianFitFDaGrad(a,μ,N,u,H,σ); #Update the gradient ∇f(x_k)
        n = norm(g);    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        xChange = α*norm(d)/maximum([norm(a),1])    # Relative and absolute change for x
    end
    return a, n
end


function backtrackingDescentMu(a::AbstractVector,μ0::AbstractMatrix, maxIter::Int,tol::Number,ρ::Number,σ::Number,H::AbstractArray,N::AbstractVector,u::Int)
    # Backtracking descent directions algorithm with steepest descent for the μ subproblem on gaussian fit problem
    μ = copy(μ0);
    k = 0;
    g = objectiveGaussianFitFDMuGrad(a,μ,N,u,H,σ); # Initialize the gradient ∇f(x_k)
    n = norm(g)/u; # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = objectiveGaussianFit(a,μ,N,u,H,σ);
    f1 = f0;
    maxjter = 12;   # Backtracking iterations
    c1 = 0.1;   # First coeficient for Wolfe conditions
    μt = copy(μ);    # Initialize trial step
    d = -copy(g); # Initialize the descent direction
    xChange = n;
    while k < maxIter && n > tol && xChange > tol^3 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g;     # Steepest descent direction
        α = 1;     # Full step
        # Line search
        μt .= μ + α*d;  # Calculate the next point
        f1 = objectiveGaussianFit(a,μt,N,u,H,σ); # Function value at the updated point
        s = d'*g;   # Directional derivative
        jter = 0;   # Backtracking iterations
        while any(f1 .> f0 .+ α*c1*s) && jter < maxjter # Backtracking steps
            α = α*ρ;    # Backtrack α by ρ
            μt .= μ + α*d;  # Update the next point
            f1 = objectiveGaussianFit(a,μt,N,u,H,σ); # Function value at the updated point
            jter += 1;
        end
        μ .= μt; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g .= objectiveGaussianFitFDMuGrad(a,μ,N,u,H,σ); #Update the gradient ∇f(x_k)
        n = norm(g)/u;    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        xChange = α*norm(d)/maximum([norm(a),1]);    # Relative and absolute change for x
    end
    return μ, n
end

function coordinateGaussianFit(αInit::AbstractVector,μInit::AbstractMatrix,σ::Number,H::AbstractArray, N::AbstractVector)
    # Coordinate descent algorithm for the gaussian fit of 3d histogram H.
    α = copy(αInit);
    μ = copy(μInit);
    maxIter = 10; # Maximum iterations for each gradient descent subproblem
    maxKter = 2000; # Maximum iterations of coordinate descent
    k = 0; # Coordinate descent iteration counter
    tol = 1e-6; # Tolerance for both coordinate descent and gradient descent
    ρ1 = 0.5; # Backtracking parameter for α descent
    ρ2 = 0.3; # Backtracking parameter for μ descent
    u = length(α); # Number of elements in the radial base (n).
    goOn = true;
    while k < maxKter && goOn
        αN, _ = backtrackingDescentAlpha(α,μ,maxIter,tol,ρ1,σ,H,N,u); # Gradient descent over α
        α .= αN;
        μN, gMu = backtrackingDescentMu(α,μ,maxIter,tol,ρ2,σ,H,N,u); # Gradient descent over μ
        μ .= μN;
        if gMu < tol                                                # If μ is optimal for this α
            gA = norm(objectiveGaussianFitFDaGrad(α,μ,N,u,H,σ));   # Check if this α is optimal for μ
            goOn = (gA > tol);                                      # In which case, stop the coordinate descent.
        end
        k += 1;
    end
    return α, μ
end