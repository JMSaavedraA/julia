# Optimization. Homework 9
# Jose Miguel Saavedra Aguilar

# Required packages: DelimitedFiles, LinearAlgebra, Colors, Images, FileIO

include("gaussFitHist.jl")

# Load the Histograms and the Image
H0, N = read3dHistogram("H_0.txt");
H1, _ = read3dHistogram("H_1.txt");
img = load("rose.bmp");
mat = channelview(img);

# Parameters
n = 3; # Radial base dimension
α = zeros(n); # Initial guess for α
μ = 2*randn(3,n); # Initial guess for μ
σ = 1; # σ parameter

# Fit α1, μ1
α0, μ0 = coordinateGaussianFit(α,μ,σ,H0,N);

# Fit α2, μ2
α1, μ1 = coordinateGaussianFit(α,μ,σ,H1,N);

# Produce the red/blue images
SH = histSegment(H0,H1,mat,N);
SF = fitSegment(α0,μ0,α1,μ1,mat,n,σ);
SB = binSegment(α0,μ0,α1,μ1,mat,n,σ,N);
sF = colorview(RGB,SF);
sH = colorview(RGB,SH);
sB = colorview(RGB,SB);

# Save the red/blue images
save("fitted.png",sF)
save("binned.png",sB)
save("histogram.png",sH)