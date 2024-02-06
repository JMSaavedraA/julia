#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 11

using LinearAlgebra, Plots

function LevenbergMarquart(R::Function,J::Function,z0::AbstractVector,μ0::Number)
    # Levenberg-Marquart non-linear least squares method
    maxIter = 5000; # Iteration limits
    τ = 1e-8;   # Tolerance
    # Trust region parameters
    γ0 = 0.25;
    γ1 = 0.75;
    η = 0.05;
    # Maximum and minimum penalty
    μMax = 1e+7;
    μMin = 1e-4;
    # Initialization
    z = copy(z0);
    r = R(z);
    j = J(z);
    f0 = r'*r/2;
    f1 = f0
    A = j'*j;
    g = j'*r;
    μ = max(μMin,min(μ0,maximum(diag(A))));
    k = 0;
    ρ = 1.0;
    # Calculate the firts step
    p = -(A+μ*I)\g;
    nP = norm(p);
    sol = nP < τ;
    while k<maxIter && !sol
        z += p; # Update the solution
        r .= R(z);  # Calculate the new residual
        f1 = r'*r/2;    # Update the value of f(x)
        ρ = (f0 - f1) / (p'*(0.5*μ*p - g)); # Calculate ρ
        if ρ < η    # We do the usual trust region analysis
            # In case the solution is very bad, we reject this solution
            z -= p;
            if μ < μMax
                # We increase the penalty parameter
                μ = min(μMax,2*μ);
                r .= R(z);
            else
                # In the case where we are at the maximum penalty, we update z to z-g
                z -= g;
                μ = max(μMin,min(μ0,minimum(diag(A)))); # We reset the penalty parameter
                # Update everything
                j .= J(z);
                A .= j'*j;
                g .= j'*r;
                f0 = f1;
            end
        else    # We accept the new point
            # Update everything
            j .= J(z);
            A .= j'*j;
            g .= j'*r;
            if ρ < γ0   # Check if we should increase the penalty
                μ = min(μMax,2*μ)
            elseif ρ >γ1    # Check if we should decrease the penalty
                μ = max(μMin,μ/3)
            end
            f0 = f1;
        end
        k += 1;
        # Calculate the next direction
        p .= -(A+μ*I)\g;
        nP = norm(p);
        sol = nP < τ    # Check stopping criteria
    end
    return z, f0, k, nP, sol
end

