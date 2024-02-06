include("gaussFitHist.jl")

H0, N0 = read3dHistogram("H_0.txt");
H1, N1 = read3dHistogram("H_1.txt");

n = 3;
α = ones(n);
μ = ones(n,n);
σ = 1;

α0, μ0 = coordinateGaussianFit(α,μ,σ,H0,N0);

α1, μ1 = coordinateGaussianFit(α,μ,σ,H1,N1);

