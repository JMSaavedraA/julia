#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 8

using LinearAlgebra

function dogleg(Δ::Number,g::AbstractVector,B::AbstractMatrix)
    # Dogleg method for the trust region problem
    α = (g'*g)/(g'*B*g); # α for the unconstrained optimum
    pU = -α*g; # The unconstrained optimum α·g
    NpU = norm(pU); # ||α·g||
    if NpU > Δ # pU not in the trust region
        pS = Δ*pU/NpU; # Δ·g/||g|| Approximate solution is Cauchy's point
        c = "C"; # Cauchy's point case
    else
        pN = -B\g; # Newton's point
        if norm(pN) <= Δ # Is Newton's point in the trust region
            pS = pN; # The exact solution is Newton's point
            c = "N"; # Newton's point case
        else
            w = pN-pU;
            # We get the polinomial w'w t^2 + pU'w t + pU'pU = Δ^2
            p1 = w'*w;
            p2 = pU'*w;
            p3 = (Δ + NpU) * (Δ - NpU);
            t = (-p2 + sqrt(p2^2 + p1*p3)) / p1; # Positive root for p1 t^2 + p2 t - p3 = 0
            pS = pU + t*w; # The dogleg point
            c = "D"; # Dogleg case
        end
    end
    return pS, c
end



function checkAcceptance(F::AbstractMatrix,g::AbstractVector,γ::Number,NF::AbstractVector)
    # Check if g is accepted into filter F
    accept = false;
    _,n = size(F);
    if n > 0
        i = 0;
        fi = F[:,1]
        while !accept && i < n
            # Check until rejected or checked all the filter
            i += 1;
            fi .= F[:,1];
            accept = any(abs.(g) .<= (abs.(fi) - γ*NF)) # |g[j]| < |v_i[j]| - γ ||v_i|| for some j
        end
    end
    return accept
end

function updateFilter(F::AbstractMatrix,g::AbstractVector,NF::AbstractVector)
    # Add g to filter F and remove all dominated vectors from filter
    _,n = size(F);
    if n > 0
        f = collect(1:n);
        fi = F[:,1];
        for i = n:-1:1
            fi .= F[:,i];
            if all(abs.(g) .< abs.(fi))
                # Remove from filter if is dominated
                deleteat!(f,i);
            end
        end
        F = F[:,f];
        NF = NF[f];
    end
    # Add g to filter F
    F = hcat(F,g);
    append!(NF,norm(g));
    return F, NF
end



function trustRegionRFTR(f::Function,x0::AbstractVector,∇::Function, H::Function)
    # Retrospective Filter Trust-Region with dogleg approximation of the model
    # Algorithm's parameters
    tol=1e-10;  # Stopping tolerance
    maxIter=50000;  # Maximum iterations
    η1 = 1/4;
    η2 = 3/4;
    γ1 = 1/4;
    γ3 = 3.5;
    ΔMax=1.0;   # Maximum trust radius
    ΔMin=1e-4;  # Minimum trust radius
    n = length(x0);  # the length of x_k
    γ = min(1e-3,1/(2*sqrt(n))); # Parameter for the filter acceptance (Similar to c1 for the Armijo condition)
    # Initializations:
    x = copy(x0);   # x_k
    xs = copy(x);   # x_k + d_k
    Δ = 1;  # Initial Δ = 1
    k = 0;  # Iteration counter
    g = ∇(x);   # ∇f (x_k)
    gs = copy(g);   # ∇f (x_k + d_k)
    B = H(x);   # ∇^2 f(x_k)
    Ng = norm(g);   # ||∇f (x_k)||
    f0 = f(x);  # f(x_k)
    fOld = f0;  # f(x_{k-1})
    f1 = f0;    # f(x_{k+1})
    F = zeros(n,0); # The filter, initilly empty
    NF = zeros(0);  # Norm of the filter's vectors
    changed = true; # Boolean indicating we updated x
    while (Ng>tol && k < maxIter)
        if changed
            # if x_k ̸= x_{k-1} we will modify B so it's positive definite:
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
        end
        d,_ = dogleg(Δ,g,B); # Dogleg approximation of the quadratic model
        xs .= x + d;    # x_k + d_k
        f1 = f(xs);     # f(x_k + d_k)
        ρ = (f1-f0)/(.5*d'*B*d + g'*d); # The trust ratio ρ
        changed = false;    # Set to not changed unless we change
        if ρ >= η1
            # Case 1. Update to x_k + d_k
            x .= xs;    # Update x_k
            g .= ∇(x);  # Update ∇(x_k)
            fOld = f0;  # Update f(x_{k-1})
            f0 = f1;    # Update f(x_k)
            Ng = norm(g);   # Update ||∇(x_k)||
            B .= H(x);  # Update ∇^2(x_k)
            changed = true; # Set Boolean to true
        else
            gs .= ∇(xs); # Calculate ∇(x_k + d_k)
            accept = checkAcceptance(F,gs,γ,NF);    # Determine if gs is acceptable for F
            if accept
                # Case 2. Update to x_k + d_k and add gs to filter
                F,NF = updateFilter(F,g,NF);
                x .= xs;    # Update x_k
                g .= gs;    # Update ∇(x_k)
                fOld = f0;  # Update f(x_{k-1})
                f0 = f1;    # Update f(x_k)
                Ng = norm(g);   # Update ||∇(x_k)||
                B .= H(x);  # Update ∇^2(x_k)
                changed = true; # Set Boolean to true
            elseif (Δ == ΔMin)
                # If we are at the minimum trust radius, we will accept even when gs was not acceptable
                x .= xs;    # Update x_k
                g .= gs;    # Update ∇(x_k)
                fOld = f0;  # Update f(x_{k-1})
                f0 = f1;    # Update f(x_k)
                Ng = norm(g);   # Update ||∇(x_k)||
                B .= H(x);  # Update ∇^2(x_k)
                changed = true; # Set Boolean to true
            end
        end
        k += 1;#Actualizamos el contador
        if changed
            # If x_k ̸= x_{k-1} we calculate the retrospective trust ratio ̂ρ
            ρ = (f0-fOld)/(-.5*d'*B*d + g'*d); # The retrospective trust ratio ̂ρ
            if ρ < η1
                # Case 1. Reduce Δ
                Δ = max(γ1*norm(d),ΔMin)
            elseif ρ > η2
                # Case 3. Increase Δ
                Δ = min(max(γ3*norm(d),Δ),ΔMax)
            end
        else
            # x_k + d_k was not good enough, so we decrease Δ
            Δ = max(γ1*Δ,ΔMin)
        end
    end
    return x, Ng, k
end