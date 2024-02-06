# Jose Miguel Saavedra Aguilar
# Master in Applied Mathematics
# CIMAT. Optimization. Homework 7

using MLDatasets, LinearAlgebra, StatsBase

# load training set
train_x, train_y = MNIST(split=:train)[:]
train_x = reshape(train_x,(784,60000))
train_ones = falses(60000)
train_ones[train_y.==1] .= 1

# load test set
test_x,  test_y  = MNIST(split=:test)[:]
test_x = reshape(test_x,(784,10000))
test_ones = falses(10000)
test_ones[test_y.==1] .= 1

function πi(β::AbstractVector,X::AbstractMatrix,n::Int)
    #π_i(β) function for all n entries of x
    Pi = 1 ./ (1 .+ exp.(-[ones(n) X']*β));
    return Pi
end

function objectiveFunction(β::AbstractVector,n::Int,X::AbstractMatrix,y::AbstractVector)
    # The objective function of the logistic regression
    Pi = πi(β,X,n); #π_i(β)
    s = (sum(y.*log.(Pi) + (1 .- y) .* log.(1 .- Pi)))/n; # 1/n ∑ y_i log (π_i) + (1-y_i)log(1-π_i)
    return s
end

function objectiveGradient(β::AbstractVector,m::Int,n::Int,X::AbstractMatrix,y::AbstractVector)
    # The gradient of the objective function of the logistic regression
    Pi = πi(β,X,n) - y; 
    g = ([ones(1,n); X]*Pi)/n; # 1/n ∑ π_i z_i
    return g
end

function objectiveHessian(β::AbstractVector,m::Int,n::Int,X::AbstractMatrix,y::AbstractVector)
    Pi = πi(β,X,n)
    H = Symmetric((([ones(1,n); X])*(((1 .- Pi).*Pi).*[ones(1,n); X]'))/n); # 1/n ∑ π_i(1-π_i) z_i z_i'
    return H
end


function pasoFijo(β0::AbstractVector, maxIter::Int,tol::Number,α::Number,n::Int)
    #Fixed step size descent directions algorithm with steepest descent
    β = copy(β0);
    F = zeros(maxIter+1); # We store f(x_k) to plot it later
    G = zeros(maxIter+1); # We store ∇f(x_k) to plot it later
    k = 0;
    m,N = size(train_x); # The size of our dataset
    r = sample(1:N,n,replace=false); # Subsample the indices
    X = train_x[:,r]; # Dataset subsample
    y = train_ones[r]; # Response variables subsample
    gk = objectiveGradient(β,m+1,n,X,y); #Initialize the gradient ∇f(x_k)
    nk = norm(gk); #Initialize the norm of the gradient ||∇f(x_k)||
    F[1] = objectiveFunction(β,n,X,y);
    G[1] = nk;
    while k<maxIter && nk>tol #Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        β -= α*gk; #Update x_k+1 as x_k - α ∇f(x_k)
        k = k+1;
        r .= sample(1:N,n,replace=false); # Subsample the indices
        X .= train_x[:,r]; # Dataset subsample
        y .= train_ones[r]; # Response variables subsample
        gk .= objectiveGradient(β,m+1,n,X,y); #Update the gradient ∇f(x_k)
        nk = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        if nk < tol
            # If we suspect we have a minimizer, resample and test again
            r .= sample(1:N,n,replace=false); # Subsample the indices
            X .= train_x[:,r]; # Dataset subsample
            y .= train_ones[r]; # Response variables subsample
            gk .= objectiveGradient(β,m+1,n,X,y); #Update the gradient ∇f(x_k)
            nk = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        end
        F[k+1] = objectiveFunction(β,n,X,y);
        G[k+1] = nk;
    end
    return β, nk, k, F, G
end

function descensoNewton(β0::AbstractVector, maxIter::Int,tol::Number,α::Number,n::Int)
    #Descent directions algorithm with Newton's point. α=1
    β = copy(β0);
    F = zeros(maxIter+1); # We store f(x_k) to plot it later
    G = zeros(maxIter+1); # We store ∇f(x_k) to plot it later
    k=0;
    m,N = size(train_x);
    r = sample(1:N,n,replace=false); # Subsample the indices
    X = train_x[:,r]; # Dataset subsample
    y = train_ones[r]; # Response variables subsample
    gk = objectiveGradient(β,m+1,n,X,y); #Initialize the gradient ∇f(x_k)
    B = objectiveHessian(β,m+1,n,X,y); #Initialize the Hessian Matrix ∇^2f(x_k)
    gz = norm(gk); #Initialize the norm of the gradient ||∇f(x_k)||
    d = copy(-gk) #Initialize the descent direction
    F[k+1] = objectiveFunction(β,n,X,y);
    G[k+1] = gz;
    while k<maxIter && gz>tol #Check if the maximum number of iterations has been reached or if ||∇f(x_k)||≈0
        B .= 0.9*objectiveHessian(β,m+1,n,X,y)+0.1I/n; #Calculate the Hessian Matrix
        d = - B\gk # Solve Bd_k = -∇h(β,X)
        β += α*d; # Update x_k+1 = x_k + α d_k
        k = k+1;
        r .= sample(1:N,n,replace=false); # Subsample the indices
        X .= train_x[:,r]; # Dataset subsample
        y .= train_ones[r]; # Response variables subsample
        gk .= objectiveGradient(β,m+1,n,X,y); #Update the gradient ∇f(x_k)
        gz = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        if gz < tol
            r .= sample(1:N,n,replace=false); # Subsample the indices
            X .= train_x[:,r]; # Dataset subsample
            y .= train_ones[r]; # Response variables subsample
            gk .= objectiveGradient(β,m+1,n,X,y); #Update the gradient ∇f(x_k)
            gz = norm(gk); #Update the norm of the gradient ||∇f(x_k)||
        end
        F[k+1] = objectiveFunction(β,n,X,y);
        G[k+1] = gz;
    end
    return β, gz, k, F, G
end

function guess(β::AbstractVector,X::AbstractMatrix)
    # Predict with weights β for dataset X
    m,n = size(X);
    Pi = πi(β,X,n); #π_i(β)
    g = BitArray(round.(Pi)); # 1(π_i(β)>0.5)
    return g
end

function guessError(g::AbstractVector,y::AbstractVector)
    # Error function, correct guesses over n
    n = length(g);
    gE = sum(abs.(g-y))/n; 
    return gE
end