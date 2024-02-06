include("Solvers.jl")
include("metodosTridiagonales.jl")

function construyeMatrizEliptica(m::Int)
    #Construye la matriz del problema Elíptico en 1D
    M=zeros(Float64, m-2, m-2)
    M[1,1]=2
    M[1,2]=-1
    M[m-2,m-3]=-1
    M[m-2,m-2]=2
    for i=2:m-3
        M[i,i-1]=-1
        M[i,i]=2
        M[i,i+1]=-1
    end
    return M
end



function construyeMatrizElipticaTridiagonal(m::Int)
    #Construye la matriz del problema Elíptico en 1D pero en forma compacta tridiagonal
    M=zeros(Float64, m-2, 3)
    M[1,2]=2
    M[1,3]=-1
    M[m-2,1]=-1
    M[m-2,2]=2
    for i=2:m-3
        M[i,1]=-1
        M[i,2]=2
        M[i,3]=-1
    end
    return M
end


function construyeVectorEliptico(m::Int)
    #Construye el vector del problema Elíptico en 1D
    v=zeros(Float64, m-2)
    v[m-2]=2*m*((m-2)/(m-1))
    for i=1:m-3
        v[i]=-2/(m-1)
    end
    return v
end


function problemaElipticoDiagonal(nodos::Int)
    #Función que construye el problema Eliptico y lo resuelve con Cholesky para matrices tridiagonales para el número de nodos indicado
    m=nodos-2
    A=construyeMatrizElipticaTridiagonal(nodos)
    b=construyeVectorEliptico(nodos)
    L=CholeskyLDLTridiagonal(A,m)
    d=zeros(Float64,m)
    #Obtenemos L, D
    for i=1:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    y=resuelveTridiagonalInferior(L,b,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    xAux=resuelveTridiagonalSuperior(L,y,m)
    #Construimos la solución completa con los extremos que ya son conocidos
    x=zeros(Float64,nodos)
    for i=1:m
        x[i+1]=xAux[i]/(nodos-1)
    end
    x[nodos]=2.0
    return x
end

function problemaElipticoJacobi(nodos::Int)
    #Función que construye el problema Eliptico y lo resuelve con el método de Jacobi para matrices tridiagonales para el número de nodos indicado
    m=nodos-2
    A=construyeMatrizElipticaTridiagonal(nodos)
    b=construyeVectorEliptico(nodos)
    h=1.6
    xInit=zeros(Float64,m)
    for i=1:m
        xInit[i]=h*i
    end
    #Resolvemos con el metodo de Jacobi para matrices tridiagonales
    xAux, Cercania, k=metodoJacobiTridiagonal(A,b,m,20000000,1e-8,xInit)
    #Construimos la solución completa con los extremos que ya son conocidos
    x=zeros(Float64,nodos)
    for i=1:m
        x[i+1]=xAux[i]/(nodos-1)
    end
    x[nodos]=2.0
    return x,Cercania,k
end

function problemaElipticoGaussSeidel(nodos::Int)
    #Función que construye el problema Eliptico y lo resuelve con el método de Gauss-Seidel para matrices tridiagonales para el número de nodos indicado
    m=nodos-2
    A=construyeMatrizElipticaTridiagonal(nodos)
    b=construyeVectorEliptico(nodos)
    h=1.6
    xInit=zeros(Float64,m)
    for i=1:m
        xInit[i]=h*i
    end
    #Resolvemos con el metodo de Gauss-Seidel para matrices tridiagonales
    xAux, Cercania, k=metodoGaussSeidelTridiagonal(A,b,m,20000000,1e-8,xInit)
    #Construimos la solución completa con los extremos que ya son conocidos
    x=zeros(Float64,nodos)
    for i=1:m
        x[i+1]=xAux[i]/(nodos-1)
    end
    x[nodos]=2.0
    return x,Cercania,k
end