using ClusterManagers
ntasks = parse(Int, ENV["SLURM_NTASKS"])
addprocs_slurm(ntasks)

using Distributed, SharedArrays, DelimitedFiles, Printf, Random

k=parse(Int64,ARGS[1]);

instance=ARGS[2];

semilla=parse(Int64,ARGS[3]);

Random.seed!(semilla);

nombreArchivo=ARGS[4]

function leerVector(nombre::String)
    b = convert(Array{Int64,1},readdlm(nombre,'\n')[2:end]);
    return b
end

function guardarMatriz(A::AbstractMatrix,m::Int,n::Int,nombre::String)
    open(nombre, "w") do f # "w" for writing JULIA
    write(f, @sprintf("%i %i\n",m,n)) # \n for newline
    for i=1:m-1
        for j=1:n-1
            write(f, @sprintf("%i ",A[i,j]))
        end
        write(f, @sprintf("%i\n",A[i,n]))
    end
    for j=1:n-1
        write(f, @sprintf("%i ",A[n,j]))
    end
    write(f, @sprintf("%i",A[n,n]))
    end
end

function loadInstance(nombre::String)
    r = readdlm(nombre,' ');
    D = convert(Array{Int64,2},r[3:end,3:end-1]);
    c = convert(Array{Int64,1},r[3:end,2]);
    M = convert(Array{Int64,1},r[2,1:2]);
    m=M[1];
    n=M[2];
    return D, c, m, n
end
    
D, c, m, n = loadInstance(instance)

@everywhere begin
    using Distributed, SharedArrays, Random, StatsBase, DelimitedFiles
    using Printf
    
    function guardarVector(x0::AbstractVector,m::Integer,nombre::String)
    open(nombre, "w") do f # "w" for writing JULIA
    write(f, string(m)*" 1\n"); # \n for newline
    for i=1:m-1
        write(f, string(x0[i])*"\n");
    end
    write(f, string(x0[m]));
    end
end

function calculateCost(x::AbstractVector,D::AbstractMatrix,c::AbstractVector)
    tC = sum(minimum(D[:,x],dims=2));
    tC += sum(c[x]);
    return tC
end

function generateRandomSolution(m::Integer,k::Integer)
    xRand=sample(Random.GLOBAL_RNG,1:m,k);
    return xRand
end
end

for i in procs()
    @spawnat i global Dtemp=D
    @spawnat i global ctemp=c
    @spawnat i global mtemp=m
    @spawnat i global ntemp=n
    @spawnat i global ktemp=k
end

z=100;

solCost=zeros(Int64,z);
sols=zeros(Int64,z,k);

pmap(1:z) do i
    try
        bestSol=collect(1:ktemp);
        bestCost=calculateCost(bestSol,Dtemp,ctemp);
        for j ∈ 1:100000
            xRand=generateRandomSolution(mtemp,ktemp);
            thisCost=calculateCost(xRand,Dtemp,ctemp);
            if thisCost<bestCost
                bestCost=thisCost;
                bestSol=xRand;
            end
        end
        guardarVector(bestSol,ktemp,"sol"*string(i)*".txt");
    catch e
        println("error en instancia", i)
        rethrow(e)
    end
end

function calculateCost(x::AbstractVector,D::AbstractMatrix,c::AbstractVector)
    tC = sum(minimum(D[:,x],dims=2));
    tC += sum(c[x]);
    return tC
end

for i=1:z
b=leerVector("sol"*string(i)*".txt")
sols[i,:]=b
solCost[i]=calculateCost(b,D,c)
rm("sol"*string(i)*".txt")
end


guardarVector(solCost,z,"costos"*nombreArchivo);
guardarMatriz(sols,z,k,"soluciones"*nombreArchivo);
