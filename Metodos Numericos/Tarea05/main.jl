include("Solvers.jl")
include("metodosTridiagonales.jl")
include("Potencia.jl")
include("problemaEliptico.jl")
using Printf

nodos=1000
A=-construyeMatrizElipticaTridiagonal(nodos+2)
n=5

println(@sprintf("Se considera el caso de %i nodos",nodos))

Pteor=zeros(Float64,nodos,n)
Lteor=zeros(Float64,n)
N=(nodos+1)
h=1/N
A=A*N*N

for j=nodos+1-n:nodos
    Lteor[nodos+1-j]=2*N*N*(cospi(j*h)-1)
    for i=1:nodos
        Pteor[i,nodos+1-j]=sinpi(i*j*h)
    end    
    nAux=norma2(Pteor[:,nodos+1-j],nodos)
    Pteor[:,nodos+1-j]=Pteor[:,nodos+1-j]/nAux
end

Aaux=reconstruyeTridiagonal(A,nodos)
x0=zeros(Float64,nodos)
for i=1:500
    x0[2*i-1]=1
    x0[2*i]=-1
end

println(@sprintf("Se obtienen los %i eigenpares más grandes con el Método de la Potencia",n))

@time P,D=metodoPotenciaNtridiagonal(A,x0,nodos,1e-13,1000000,n)

guardarMatriz(P,nodos,n,"PPotencia.txt")
guardarVector(D,n,"lambdaPotencia.txt")

EL=zeros(Float64,n)
EP=zeros(Float64,n)
for i=1:n
    EP[n+1-i]=norma2(P[:,i]-Pteor[:,i],nodos)
    println(@sprintf("El error del %iº eigenvector más grande es %0.5e",i,EP[n+1-i]))
    EL[n+1-i]=abs(Lteor[i]-D[i])
    println(@sprintf("El error del %iº eigenvalor más grande es %0.5e",i,EL[n+1-i]))
end


println(@sprintf("Se obtienen los %i eigenpares más grandes con el Método de la Potencia Inversa",n))

for j=1:n
    Lteor[j]=N*N*2*(cospi(j*h)-1)
    for i=1:nodos
        Pteor[i,j]=sinpi(i*j*h)
    end
    nAux=norma2(Pteor[:,j],nodos)
    Pteor[:,j]=Pteor[:,j]/nAux
end


@time P,D=metodoPotenciaInversaNtridiagonal(A,x0,nodos,1e-13,1000000,n)

guardarMatriz(P,nodos,n,"PPotenciaInversa.txt")
guardarVector(D,n,"lambdaPotenciaInversa.txt")

EL=zeros(Float64,n)
EP=zeros(Float64,n)
for i=1:n
    EP[i]=norma2(P[:,i]-Pteor[:,i],nodos)
    println(@sprintf("El error del %iº eigenvector más pequeño es %0.5e",i,EP[i]))
    EL[i]=abs(Lteor[i]-D[i])
    println(@sprintf("El error del %iº eigenvalor más pequeño es %0.5e",i,EL[i]))
end