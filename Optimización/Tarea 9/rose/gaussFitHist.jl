# Optimization. Homework 9
# Jose Miguel Saavedra Aguilar

# Required packages: DelimitedFiles, LinearAlgebra, Colors, Images, FileIO


using DelimitedFiles, LinearAlgebra, Colors, Images, FileIO

function read3dHistogram(fileName::String)
    # Auxiliary function that reads a file containing a segmented 3d histogram.
    N = readdlm(fileName)[1,:];
    H = reshape(readdlm(fileName,'\n')[2:end],N[1],N[2],N[3]);
    return H, N
end



function objectiveGaussianFit(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # Gaussian fit objective function for 3d histogram H
    s = 0.0;
    for k = 1:N[3]
        for j = 1:N[2]
            for i = 1:N[1]
                ss = objectiveSum(α,μ,n,σ,[i,j,k]);
                s += (H[i,j,k] - ss)^2;
            end
        end
    end
    s = sqrt(s);
    return s
end

function objectiveSum(α::AbstractVector,μ::AbstractMatrix,n::Int,σ::Number,c::AbstractVector)
    # Inner sum for the Gaussian fit objective function for 3d histogram H
    ss = 0.0;
    for m = 1:n
        ss += α[m]*exp(-((norm(c - μ[:,m])/σ)^2)/2);
    end
    return ss
end

function objectiveGaussianFitFDaGrad(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # ∇ f(α,μ) with respect to α approximation via central Finite Differences
    h = 1e-8 * max(norm(α)/sqrt(n),1);
    ∇ = zeros(n);
    αL = copy(α);
    αR = copy(α);
    h2 = 2*h;
    for i=1:n
        αL[i] -= h; # Backward step
        αR[i] += h; # Forward step
        ∇[i] = (objectiveGaussianFit(αR,μ,N,n,H,σ) - objectiveGaussianFit(αL,μ,N,n,H,σ))/(h2); # Finite difference
        αL[i] = α[i]; # Restore to original value
        αR[i] = α[i]; # Restore to original value
    end
    return ∇
end

function objectiveGaussianFitFDMuGrad(α::AbstractVector,μ::AbstractMatrix,N::AbstractVector,n::Int,H::AbstractArray,σ::Number)
    # ∇ f(α,μ) with respect to μ approximation via central Finite Differences
    h = 1e-8 * max(norm(μ)/(n),1);
    ∇ = zeros(3,n);
    mL = copy(μ);
    mR = copy(μ);
    h2 = 2*h;
    for j=1:n
        for i=1:3
            mL[i,j] = μ[i,j] - h; # Backward step
            mR[i,j] = μ[i,j] + h; # Forward step
            ∇[i,j] = (objectiveGaussianFit(α,mR,N,n,H,σ) - objectiveGaussianFit(α,mL,N,n,H,σ))/(h2); # Finite difference
            mL[i,j] = μ[i,j]; # Restore to original value
            mR[i,j] = μ[i,j]; # Restore to original value
        end
    end
    return ∇
end

function backtrackingDescentAlpha(a0::AbstractVector,μ::AbstractMatrix, maxIter::Int,tol::Number,ρ::Number,σ::Number,H::AbstractArray,N::AbstractVector,n::Int)
    # Backtracking descent directions algorithm with steepest descent for the α subproblem on gaussian fit problem
    a = copy(a0);
    k = 0;
    g = objectiveGaussianFitFDaGrad(a,μ,N,n,H,σ); # Initialize the gradient ∇f(x_k)
    nG = norm(g); # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = objectiveGaussianFit(a,μ,N,n,H,σ);
    f1 = f0;
    maxjter = 15;   # Backtracking iterations
    c1 = 1e-3;   # First coeficient for Wolfe conditions
    at = copy(a)    # Initialize trial step
    d = copy(g) # Initialize the descent direction
    xChange = nG;
    while k < maxIter && nG > tol && xChange > tol^3 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g/max(nG,1e-4);     # Steepest descent direction
        α  = 1;     # Full step
        # Line search
        at .= a + α*d;  # Calculate the next point
        f1 = objectiveGaussianFit(at,μ,N,n,H,σ); # Function value at the updated point
        s = d'*g;   # Directional derivative
        jter = 0;   # Backtracking iterations
        while f1 > f0 + α*c1*s && jter < maxjter # Backtracking steps
            α = α*ρ;    # Backtrack α by ρ
            at .= a + α*d;  # Update the next point
            f1 = objectiveGaussianFit(at,μ,N,n,H,σ); # Function value at the updated point
            jter += 1;
        end
        a .= at; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g .= objectiveGaussianFitFDaGrad(a,μ,N,n,H,σ); #Update the gradient ∇f(x_k)
        nG = norm(g);    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        xChange = α*norm(d)/maximum([norm(a),1])    # Relative and absolute change for x
    end
    return a, nG
end


function backtrackingDescentMu(a::AbstractVector,μ0::AbstractMatrix, maxIter::Int,tol::Number,ρ::Number,σ::Number,H::AbstractArray,N::AbstractVector,n::Int)
    # Backtracking descent directions algorithm with steepest descent for the μ subproblem on gaussian fit problem
    μ = copy(μ0);
    k = 0;
    g = objectiveGaussianFitFDMuGrad(a,μ,N,n,H,σ); # Initialize the gradient ∇f(x_k)
    nG = norm(g)/n; # Initialize the norm of the gradient ||∇f(x_k)||
    f0 = objectiveGaussianFit(a,μ,N,n,H,σ);
    f1 = f0;
    maxjter = 15;   # Backtracking iterations
    c1 = 0.1;   # First coeficient for Wolfe conditions
    μt = copy(μ);    # Initialize trial step
    d = -copy(g); # Initialize the descent direction
    xChange = nG;
    while k < maxIter && nG > tol && xChange > tol^3 # Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        d .= -g/max(nG,1e-4);     # Steepest descent direction
        α = 1;     # Full step
        # Line search
        μt .= μ + α*d;  # Calculate the next point
        f1 = objectiveGaussianFit(a,μt,N,n,H,σ); # Function value at the updated point
        s = d'*g;   # Directional derivative
        jter = 0;   # Backtracking iterations
        while any(f1 .> f0 .+ α*c1*s) && jter < maxjter # Backtracking steps
            α = α*ρ;    # Backtrack α by ρ
            μt .= μ + α*d;  # Update the next point
            f1 = objectiveGaussianFit(a,μt,N,n,H,σ); # Function value at the updated point
            jter += 1;
        end
        μ .= μt; # Update x_k+1 as x_k + α d_k
        k = k+1;
        g .= objectiveGaussianFitFDMuGrad(a,μ,N,n,H,σ); #Update the gradient ∇f(x_k)
        nG = norm(g)/n;    #Update the norm of the gradient ||∇f(x_k)||
        f0 = f1;    # Update the function value at the new point
        xChange = α*norm(d)/maximum([norm(a),1]);    # Relative and absolute change for x
    end
    return μ, nG
end

function coordinateGaussianFit(αInit::AbstractVector,μInit::AbstractMatrix,σ::Number,H::AbstractArray, N::AbstractVector)
    # Coordinate descent algorithm for the gaussian fit of 3d histogram H.
    α = copy(αInit);
    μ = copy(μInit);
    maxIter = 50; # Maximum iterations for each gradient descent subproblem for α
    maxJter = 2; # Maximum iterations for each gradient descent subproblem for μ
    maxKter = 20000; # Maximum iterations of coordinate descent
    k = 0; # Coordinate descent iteration counter
    tol = 1e-8; # Tolerance for both coordinate descent and gradient descent
    ρ1 = 0.5; # Backtracking parameter for α descent
    ρ2 = 0.3; # Backtracking parameter for μ descent
    u = length(α); # Number of elements in the radial base (n).
    goOn = true;
    while k < maxKter && goOn
        αN, _ = backtrackingDescentAlpha(α,μ,maxIter,tol,ρ1,σ,H,N,u); # Gradient descent over α
        α .= αN;
        μN, gMu = backtrackingDescentMu(α,μ,maxJter,tol,ρ2,σ,H,N,u); # Gradient descent over μ
        μ .= μN;
        if gMu < tol                                                # If μ is optimal for this α
            gA = norm(objectiveGaussianFitFDaGrad(α,μ,N,u,H,σ));   # Check if this α is optimal for μ
            goOn = (gA > tol);                                      # In which case, stop the coordinate descent.
        end
        k += 1;
    end
    return α, μ
end

function Hist2Label(H0::AbstractArray,H1::AbstractArray)
    # Check if H_1(c) < H_2(c) for all c in the histogram bins. We shall call this labels
    ϵ = 1e-2; # ϵ
    den = H0 + H1 .+ ϵ;
    hL = (H0 .+ ϵ) ./ den; # H_1(c)
    hR = (H1 .+ ϵ) ./ den; # H_2(c)
    B = hL .< hR
    return B
end

function histSegment(H0::AbstractArray,H1::AbstractArray,mat::AbstractArray,N::AbstractVector)
    # Creates the red/blue image that corresponds to the histogram
    p,q = size(mat[1,:,:]); # Image size
    S = zeros(3,p,q); # New red and blue image
    b = N[1]; # Number of bins
    B = Hist2Label(H0,H1); # Labels assigned by the histograms
    v = zeros(Int,3)
    for j = 1:q
        for i = 1:p
            v .= Int.(round.(mat[:,i,j]*(b-1))) .+ 1 # Assign which bin corresponds to the pixel
            if B[v[1],v[2],v[3]]
                S[3,i,j] = 1; # Blue
            else
                S[1,i,j] = 1; # Red
            end
        end
    end
    return S
end

function bin2Label(α0::AbstractVector,μ0::AbstractMatrix,α1::AbstractVector,μ1::AbstractMatrix,n::Int,σ::Number,N::AbstractVector)
    # Assign the label to the bins via the objective function, so we can compare roughly if the fitment is correct. A bad fit will be very different to the Histogram assigned labels
    ϵ = 1e-2;
    F0 = zeros(N[1],N[2],N[3]);
    F1 = copy(F0);
    for k = 1:N[3]
        for j = 1:N[2]
            for i = 1:N[1]
                c = [i,j,k]; # Bin [i,j,k]
                F0[i,j,k] = objectiveSum(α0,μ0,n,σ,c) + ϵ; # f(bin,α1,μ1)
                F1[i,j,k] = objectiveSum(α1,μ1,n,σ,c) + ϵ; # f(bin,α2,μ2)
            end
        end
    end
    B = F0 .< F1 # Equivalent to F(bin,α1,μ1) < F(bin,α2,μ2)
    return B
end

function binSegment(α0::AbstractVector,μ0::AbstractMatrix,α1::AbstractVector,μ1::AbstractMatrix,mat::AbstractArray,n::Int,σ::Number,N::AbstractVector)
    # Creates the red/blue image that corresponds to the bins with the objective function
    p,q = size(mat[1,:,:]); # Size of the image
    S = zeros(3,p,q); # New red/blue image
    b = N[1]; # Number of bins
    B = bin2Label(α0,μ0,α1,μ1,n,σ,N); # Labels obtained from the objective function
    v = zeros(Int,3);
    for j = 1:q
        for i = 1:p
            v .= Int.(round.(mat[:,i,j]*(b-1))) .+ 1 # Assign which bin corresponds to the pixel
            if B[v[1],v[2],v[3]]
                S[3,i,j] = 1; # Blue
            else
                S[1,i,j] = 1; # Red
            end
        end
    end
    return S
end

function fitSegment(α0::AbstractVector,μ0::AbstractMatrix,α1::AbstractVector,μ1::AbstractMatrix,mat::AbstractArray,n::Int,σ::Number)
    # Creates the red/blue image that corresponds for each pixel with the objective function
    ϵ = 1e-2; # ϵ
    p,q = size(mat[1,:,:]); # Size of the image
    S = zeros(3,p,q); # New red/blue image
    b = N[1]; # Number of bins
    for j = 1:q
        for i = 1:p
            c = mat[:,i,j]*b # RGB of the pixel
            F0 = objectiveSum(α0,μ0,n,σ,c) + ϵ; # f(c,α1,μ1)
            F1 = objectiveSum(α1,μ1,n,σ,c) + ϵ; # f(c,α2,μ2)
            if F0 < F1
                S[3,i,j] = 1; # Blue
            else
                S[1,i,j] = 1; # Red
            end
        end
    end
    return S
end