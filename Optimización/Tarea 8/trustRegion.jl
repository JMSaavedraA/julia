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

function naive(Δ::Number,g::AbstractVector,B::AbstractMatrix)
  # Not-so-naive method for the trust region problem
  α = (g'*g)/(g'*B*g); # α for the unconstrained optimum
  pU = -α*g; # The unconstrained optimum α·g
  NpU = norm(pU); # ||α·g||
  if NpU > Δ # pU not in the trust region
    pS = Δ*pU/NpU; # Δ·g/||g|| Approximate solution is Cauchy's point
    c = "C"; # Cauchy's point case
  else
    pN = -B\g; # Newton's point
    NpN = norm(pN)
    if NpN <= Δ # Is Newton's point in the trust region
      pS = pN; # The exact solution is Newton's point
      c = "N"; # Newton's point case
    else
      pS = Δ*pN/NpN; # Δ·pN/||pN|| The best solution on Newton's point direction
      c = "D"; # Descent direction case
    end
  end
  return pS, c
end



function trustRegionDogleg(f::Function,x0::AbstractVector,∇::Function, H::Function)
  # Trust-Region method with dogleg approximation of the model
  # Algorithm's parameters
  tol=1e-10;  # Stopping tolerance
  maxIter=50000;  # Maximum iterations
  η1 = 1/8;
  η2 = 7/8;
  λ1 = 1/4;
  λ2 = 2;
  ΔMax=1.0;   # Maximum trust radius
  ΔMin=1e-4;  # Minimum trust radius
  n = length(x0);  # the length of x_k
  x = copy(x0);   # x_k
  xs = copy(x);   # x_k + d_k
  Δ = 1;  # Initial Δ = 1
  k = 0;  # Iteration counter
  g = ∇(x);   # ∇f (x_k)
  B = H(x);   # ∇^2 f(x_k)
  Ng = norm(g);   # ||∇f (x_k)||
  f0 = f(x);  # f(x_k)
  f1 = f0;    # f(x_{k+1})
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
    d,c = dogleg(Δ,g,B); # Dogleg approximation of the quadratic model
    xs .= x + d;    # x_k + d_k
    f1 = f(xs);     # f(x_k + d_k)
    ρ = (f1-f0)/(.5*d'*B*d + g'*d); # The trust ratio ρ
    changed = true; # Set to changed unless we don't change
    if ρ < η1 # Case 1. Reject due to bad trust ratio
      if Δ == ΔMin
        # If we are at the minimum trust radius, we will accept even when ρ was bad
        x .= xs;    # Update x_k
        g .= gs;    # Update ∇(x_k)
        f0 = f1;    # Update f(x_k)
        Ng = norm(g);   # Update ||∇(x_k)||
        B .= H(x);  # Update ∇^2(x_k)
      else
        # Reduce the trust radius
        Δ = max(λ1 * Δ , ΔMin);
        changed = false;  # Boolean to not changed
      end
    elseif (ρ > η2 && c != "N")
      # Case 3. good ρ and ||d_k||=Δ. Update to x_k + d_k and increase Δ
      x .= xs;    # Update x_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
      Δ = min(λ2 * Δ , ΔMax); # Increase Δ
    else
      # Case 2. good ρ and ||d_k||<Δ. Update to x_k + d_k
      x .= xs;    # Update x_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
    end
    k += 1;
  end
  return x, Ng, k
end

function trustRegionDogleg(f::Function,x0::AbstractVector,∇::Function, H::Function)
  # Trust-Region method with not-so-naive approximation of the model
  # Algorithm's parameters
  tol=1e-10;  # Stopping tolerance
  maxIter=50000;  # Maximum iterations
  η1 = 1/8;
  η2 = 7/8;
  λ1 = 1/4;
  λ2 = 2;
  ΔMax=1.0;   # Maximum trust radius
  ΔMin=1e-4;  # Minimum trust radius
  n = length(x0);  # the length of x_k
  x = copy(x0);   # x_k
  xs = copy(x);   # x_k + d_k
  Δ = 1;  # Initial Δ = 1
  k = 0;  # Iteration counter
  g = ∇(x);   # ∇f (x_k)
  B = H(x);   # ∇^2 f(x_k)
  Ng = norm(g);   # ||∇f (x_k)||
  f0 = f(x);  # f(x_k)
  f1 = f0;    # f(x_{k+1})
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
    d,c = naive(Δ,g,B); # Not-so-naive approximation of the quadratic model
    xs .= x + d;    # x_k + d_k
    f1 = f(xs);     # f(x_k + d_k)
    ρ = (f1-f0)/(.5*d'*B*d + g'*d); # The trust ratio ρ
    changed = true; # Set to changed unless we don't change
    if ρ < η1 # Case 1. Reject due to bad trust ratio
      if Δ == ΔMin
        # If we are at the minimum trust radius, we will accept even when ρ was bad
        x .= xs;    # Update x_k
        g .= gs;    # Update ∇(x_k)
        f0 = f1;    # Update f(x_k)
        Ng = norm(g);   # Update ||∇(x_k)||
        B .= H(x);  # Update ∇^2(x_k)
      else
        # Reduce the trust radius
        Δ = max(λ1 * Δ , ΔMin);
        changed = false;  # Boolean to not changed
      end
    elseif (ρ > η2 && c != "N")
      # Case 3. good ρ and ||d_k||=Δ. Update to x_k + d_k and increase Δ
      x .= xs;    # Update x_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
      Δ = min(λ2 * Δ , ΔMax); # Increase Δ
    else
      # Case 2. good ρ and ||d_k||<Δ. Update to x_k + d_k
      x .= xs;    # Update x_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
    end
    k += 1;
  end
  return x, Ng, k
end