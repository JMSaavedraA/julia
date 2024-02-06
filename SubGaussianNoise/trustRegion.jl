#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics


include("ECNoise.jl")
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

function ∇(x0::AbstractVector, f::Function, h::AbstractFloat)
  # ∇ f(α,μ) with respect to α approximation via central Finite Differences
  g = copy(x0);
  xL = copy(x0);
  xR = copy(x0);
  h2 = 2*h;
  for i in eachindex(x0)
      xL[i] -= h; # Backward step
      xR[i] += h; # Forward step
      g[i] = (f(xR) - f(xL))/(h2); # Finite difference
      xL[i] = x0[i]; # Restore to original value
      xR[i] = x0[i]; # Restore to original value
  end
  return g
end

function trustRegionDogleg(f::Function,x0::AbstractVector,σ::Number)
  # Trust-Region method for noisy functions with dogleg approximation of the model and BFGS Update
  # Algorithm's parameters
  maxIter = 10000;  # Maximum iterations
  h = 10^(floor(log10(sqrt(σ))));
  η1 = 1/8;
  η2 = 7/8;
  λ1 = 1/4;
  λ2 = 2;
  s = 5;
  r = 2*s/(1 - η2);
  ΔMax = 1.0;   # Maximum trust radius
  ΔMin = 1e-4;  # Minimum trust radius
  n = length(x0);  # the length of x_k
  x = copy(x0);   # x_k
  xs = copy(x);   # x_k + d_k
  Δ = 1;  # Initial Δ = 1
  k = 0;  # Iteration counter
  g = ∇(x,f,h);   # ∇f (x_k)
  Ng = norm(g);   # ||∇f (x_k)||
  gs = copy(g)/Ng; # ∇f(x_k + d_k)
  funNoise(t) = f(x0 + t*gs); # Auxiliary function to calculate the noise of f
  y = copy(g); # ∇f(x_k+1) - # ∇f(x_k)
  B = Matrix(I, n,n);   # Initial guess is the Identity Matrix, a SPD Matrix
  f0 = f(x);  # f(x_k)
  f1 = f0;    # f(x_{k+1})
  tol = 1e-6 + sqrt(n)*σ/(h) # Tolerance update to consider the noise of ∇f(x_k) approximation
  while (Ng>(tol)  && k < maxIter)
    d,c = dogleg(Δ,g,B); # Dogleg approximation of the quadratic model
    xs .= x + d;    # x_k + d_k
    f1 = f(xs);     # f(x_k + d_k)
    gs .= ∇(xs,f,h); # ∇f(x_k + d_k)
    ρ = ((f1-f0 - r*σ*s)/(.5*d'*B*d + g'*d - r*σ*s)); # The trust ratio ρ
    if ρ < η1 # Case 1. Reject due to bad trust ratio
      if Δ == ΔMin
        # If we are at the minimum trust radius, we will accept even when ρ was bad
        x .= xs;    # Update x_k
        y .= gs-g;  # Update y_k
        g .= gs;    # Update ∇(x_k)
        f0 = f1;    # Update f(x_k)
        Ng = norm(g);   # Update ||∇(x_k)||
        u = B*d;    # Auxiliary B*d_k
        B += (y*y') / (dot(y,d)) - (u*u') / dot(d,u);  # BFGS Update of the Hessian Matrix
      else
        # Reduce the trust radius
        Δ = max(λ1 * Δ , ΔMin);
      end
    elseif (ρ > η2 && c != "N")
      # Case 3. good ρ and ||d_k||=Δ. Update to x_k + d_k and increase Δ
      x .= xs;    # Update x_k
      y .= gs-g;  # Update y_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      u = B*d;    # Auxiliary B*d_k
      B += (y*y') / (dot(y,d)) - (u*u') / dot(d,u);  # BFGS Update of the Hessian Matrix
      Δ = min(λ2 * Δ , ΔMax); # Increase Δ
    else
      # Case 2. good ρ and ||d_k||<Δ. Update to x_k + d_k
      x .= xs;    # Update x_k
      y .= gs-g;  # Update y_k
      g .= gs;    # Update ∇(x_k)
      f0 = f1;    # Update f(x_k)
      Ng = norm(g);   # Update ||∇(x_k)||
      u = B*d;    # Auxiliary B*d_k
      B += (y*y') / (dot(y,d)) - (u*u') / dot(d,u);  # BFGS Update of the Hessian Matrix
    end
    k += 1;
  end
  return x, Ng, k, tol
end
