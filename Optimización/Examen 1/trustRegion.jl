#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Partial Exam 1

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

function exactTrustRegion(Δ::Number,g::AbstractVector,B::AbstractMatrix,tol::Number)
  # Exact solution to the trust region subproblem via Newton's method
  pS = - B \ g; # Start with λ=0, ie Newton's point
  nPs = norm(pS);
  q = B \ pS;
  # Check if Newton's point is feasible
  if nPs < Δ
    c = "N";
  else
    # Update to (p⊤p / p⊤q) (||p||-Δ/Δ)
    u = (nPs^2)*(nPs - Δ)/(pS'*q*Δ);
    λ = u;
    j = 0;
    while abs(u) > tol && j < 3
      pS .= -(B + λ*I) \ (g / Δ); # p_k = -1/Δ (B+λI)^-1 g
      nPs = norm(pS);
      q .= (B+λ*I) \ pS; # q_k = (B+λI)^-1 p_k
      u = (nPs^2)*(nPs - 1)/(pS'*q); # The update u_k=(p_k⊤p_k / p_k⊤q_k) (||p_k||-1)
      λ += u; # Update λ_k+1 = λ_k + u_k
      j +=1;
    end
    c = "C";
    pS = Δ*pS/norm(pS); # Sometimes pS won't be in the trust region and that is very important for 
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
    ρ = (f0-f1)/abs(.5*d'*B*d + g'*d); # The trust ratio ρ
    changed = true; # Set to changed unless we don't change
    if ρ < η1 # Case 1. Reject due to bad trust ratio
      if Δ == ΔMin
        # If we are at the minimum trust radius, we will accept even when ρ was bad
        x .= xs;    # Update x_k
        g .= ∇(xs);    # Update ∇(x_k)
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
      g .= ∇(xs);    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
      Δ = min(λ2 * Δ , ΔMax); # Increase Δ
    else
      # Case 2. good ρ and ||d_k||<Δ. Update to x_k + d_k
      x .= xs;    # Update x_k
      g .= ∇(xs);    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
    end
    k += 1;
  end
  convergence = (Ng < tol)
  return x, Ng, k, convergence
end


function trustRegionExacta(f::Function,x0::AbstractVector,∇::Function, H::Function)
  # Trust-Region method with not-so-naive approximation of the model
  # Algorithm's parameters
  tol=1e-4;  # Stopping tolerance
  maxIter=100000;  # Maximum iterations
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
    d,c = exactTrustRegion(Δ,g,B,tol); # Get the "exact" solution to the subproblem
    xs .= x + d;    # x_k + d_k
    f1 = f(xs);     # f(x_k + d_k)
    ρ = (f0-f1)/abs(.5*d'*B*d + g'*d); # The trust ratio ρ
    changed = true; # Set to changed unless we don't change
    if ρ < η1 # Case 1. Reject due to bad trust ratio
      if Δ == ΔMin
        # If we are at the minimum trust radius, we will accept even when ρ was bad
        x .= xs;    # Update x_k
        g .= ∇(xs);    # Update ∇(x_k)
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
      g .= ∇(xs);    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
      Δ = min(λ2 * Δ , ΔMax); # Increase Δ
    else
      # Case 2. good ρ and ||d_k||<Δ. Update to x_k + d_k
      x .= xs;    # Update x_k
      g .= ∇(xs);    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      B .= H(x);  # Update ∇^2(x_k)
    end
    k += 1;
  end
  convergence = (Ng < tol)
  return x, Ng, k, convergence
end