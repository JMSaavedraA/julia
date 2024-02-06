include("Solvers.jl")
include("problemaEliptico.jl")
include("metodoJacobi.jl")
using Printf


#Este progrma resuelve el problema de la tarea, solo cambiamos el número de nodos en la siguiente linea si deseamos. Para n=100 toma aproximadamente 2 segundos, para n=1000 toma muchas horas
nodos=100
n=nodos
A=-construyeMatrizEliptica(nodos+2)
Pteor=zeros(Float64,nodos,nodos)
Lteor=zeros(Float64,nodos)
N=(nodos+1)
h=1/N
A=A*N*N

for j=1:nodos
    Lteor[j]=N*N*2*(cospi(j*h)-1)
    for i=1:nodos
        Pteor[i,j]=sinpi(i*j*h)
    end
    nAux=norma2(Pteor[:,j],nodos)
    Pteor[:,j]=Pteor[:,j]/nAux
end


@time P,D,k,Aaux=eigenJacobi(A,nodos,1e-10,8000000)

guardarMatriz(P,nodos,nodos,@sprintf("PJacobi%i.txt",nodos))
guardarVector(D,nodos,@sprintf("lambdaJacobi%i.txt",nodos))

EL=zeros(Float64,nodos)
EP=zeros(Float64,nodos)
for i=1:nodos
    e1=norma2(P[:,i]-Pteor[:,i],nodos)
    e2=norma2(P[:,i]+Pteor[:,i],nodos)
    EP[i]=minimum([e1,e2])
    #Las lineas siguientes que están comentadas imprimen el valor del error de cada eigenpar, pero es mucha información incluso para n=100, pero si desean ver los errores se puede descomentar
#    println(@sprintf("El error del %iº eigenvector más pequeño es %0.5e",i,EP[i]))
    EL[i]=abs(Lteor[i]-D[i])
#    println(@sprintf("El error del %iº eigenvalor más pequeño es %0.5e",i,EL[i]))
end