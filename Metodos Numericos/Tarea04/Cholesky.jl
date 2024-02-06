include("Solvers.jl")
function CholeskyLDL(Aaux::AbstractMatrix,m::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A simetrica de tamaño mxm
    L=similar(Aaux)
    for i=1:m
        sumaExterna=0.0
        for k=1:i-1
            vAux=0.0
            suma=0.0
            for j=1:k-1
                suma+=L[i,j]*L[j,j]*L[k,j]
            end
            vAux=(Aaux[i,k]-suma)
            L[i,k]=vAux/L[k,k]
            sumaExterna+=vAux*vAux/L[k,k]
        end
        L[i,i]=Aaux[i,i]-sumaExterna
    end
    return L
end

function CholeskyLDLNdiagonal(Aaux::AbstractMatrix,m::Int,n::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A bandada simetrica de tamaño de banda n y A de tamaño mxm
    L=zeros(Float64,m,n) #Guardamos en una matriz mxn solo la parte inferior de L, igualmente así se espera que sea A
    for i=1:m
        iShift=minimum([n-i,0])#Calculamos que entrada de L estamos calculando
        sumaExterna=0.0
        for k=1-iShift:i-1
            kShift=minimum([n-k,0])
            vAux=0.0
            suma=0.0
            for j=1-iShift:k-1
                jShift=minimum([n-j,0])
                suma+=L[i,j+iShift]*L[j,j+jShift]*L[k,j+kShift]
            end
            vAux=(Aaux[i,k+iShift]-suma)
            L[i,k+iShift]=vAux/L[k,k+kShift]
            sumaExterna+=vAux*vAux/L[k,k+kShift]
        end
        L[i,i+iShift]=Aaux[i,i+iShift]-sumaExterna
    end
    return L
end

function formaNdiagonal(A::AbstractMatrix,m::Int,n::Int)
    #Función que nos regresa la diagonal de la matriz A bandada de tamaño de banda n y tamaño de A mxm
    Aaux=zeros(Float64, m, n)
    for i=1:m
        iShift=minimum([n-i,0])        
        for k=1-iShift:i
            Aaux[i,k+iShift]=A[i,k]
        end
    end
    return Aaux
end

function reconstruyeNdiagonalInferior(A::AbstractMatrix,m::Int,n::Int)
    #Nos permite pasar una matriz bandada Triangular inferior de su forma compacta a su forma completa (con ceros) 
    Aaux=zeros(Float64, m, m)
    for i=1:m
        iShift=minimum([n-i,0])        
        for k=1-iShift:i
            Aaux[i,k]=A[i,k+iShift]
        end
    end
    return Aaux
end


function formaNdiagonalSuperior(B::AbstractMatrix,m::Int,n::Int)
    #Nos permite pasar una matriz bandada Triangular superior de su forma completa a su forma compacta (sin ceros) 
    Baux=zeros(Float64, m, n)
    for i=1:m-n
        Baux[i,1:n]=B[i,i:i+n-1]
    end
    Baux[m-n+1:m,1:n]=B[m-n+1:m,m-n+1:m]
    return Baux
end

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

function transponerNdiagonalInferior(A::AbstractMatrix,m::Int,n::Int)
    #Nos permite pasar una matriz bandada Triangular inferior a una Triangular superior ambas de forma compacta.
    J=formaNdiagonalSuperior(transpose(reconstruyeNdiagonalInferior(A,m,n)),m,n)
    return J
end

function construyeMatrizElipticaDiagonal(m::Int)
    #Construye la matriz del problema Elíptico en 1D pero en forma compacta tridiagonal
    M=zeros(Float64, m-2, 2)
    M[1,1]=2
    M[1,2]=0.0
    M[m-2,1]=-1
    M[m-2,2]=2
    for i=2:m-3
        M[i,1]=-1
        M[i,2]=2
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

function resuelveNdiagonalInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,n::Int)
    #Eliminación gaussiana para Aaux triangular inferior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    for i=1:n
        x[i]=(baux[i]-sum(x[1:i-1].*Aaux[i,1:i-1]))/Aaux[i,i]
    end
    for i=n+1:m
        x[i]=(baux[i]-sum(x[i-n+1:i-1].*Aaux[i,1:n-1]))/Aaux[i,n]
    end
    return x
end

function resuelveNdiagonalSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,n::Int)
    #Eliminación gaussiana para Aaux triangular superior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    for i in 1:n
        j=m-i+1
        k=n-i+1
        x[j]=(baux[j]-sum(x[j+1:end].*Aaux[j,k+1:end]))/Aaux[j,k]
    end
    for i=n+1:m
        j=m-i+1
        x[j]=(baux[j]-sum(x[j+1:j+n-1].*Aaux[j,2:n]))/Aaux[j,1]
    end
    return x
end

function norma2(x::AbstractVector,m::Int)
    #La norma del vector x de tamaño m
    suma=0
    for i=1:m
        suma+=x[i]*x[i]
    end
    n=sqrt(suma)
    return n
end

function metodoJacobi(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Jacobi para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        xNew=similar(xOld)
        for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xOld[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xOld[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoGaussSeidel(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Gauss-Seidel para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xNew[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xNew[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoJacobiTridiagonal(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Jacobi para un problema en general Ax=y, A tridiagonal simétrica de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        xNew=similar(xOld)
        xNew[1]=(y[1]-A[2,1]*xOld[2])/A[1,1]
        for j=2:m-1
            xNew[j]=(y[j]-A[j,1]*xOld[j-1]-A[j+1,1]*xOld[j+1])/A[j,2]
        end
        xNew[m]=(y[m] - A[m,1]*xOld[m-1])/A[m,2]
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i#Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoGaussSeidelTridiagonal(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Gauss-Seidel para un problema en general Ax=y, A tridiagonal simétrica de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        xNew[1]=(y[1]-A[2,1]*xNew[2])/A[1,1]
        for j=2:m-1
            xNew[j]=(y[j]-A[j,1]*xNew[j-1]-A[j+1,1]*xNew[j+1])/A[j,2]
        end
        xNew[m]=(y[m] - A[m,1]*xNew[m-1])/A[m,2]
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i#Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function resuelveDiagonalCompacta(A::AbstractMatrix,b::AbstractVector,m::Int,n::Int)
    #Solucion de Dx=b, donde D es la diagonal de A, una matriz bandada de banda n y tamaño mxm
    x=similar(b)
    for i=1:n
        x[i]=b[i]/A[i,i]
    end
    for i=n+1:m
        x[i]=b[i]/A[i,n]
    end
    return x
end

function problemaElipticoDiagonal(nodos::Int)
    #Función que construye el problema Eliptico y lo resuelve con Cholesky para matrices tridiagonales para el número de nodos indicado
    m=nodos-2
    A=construyeMatrizElipticaDiagonal(nodos)
    b=construyeVectorEliptico(nodos)
    L=CholeskyLDLNdiagonal(A,m,2)
    d=zeros(Float64,m)
    #Obtenemos L, D
    d[1]=L[1,1]
    L[1,1]=1.0
    for i=2:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    #Obtenemos L^t
    Ltrasposed=transponerNdiagonalInferior(L,m,2)
    y=resuelveNdiagonalInferior(L,b,m,2)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    xAux=resuelveNdiagonalSuperior(Ltrasposed,y,m,2)
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
    A=construyeMatrizElipticaDiagonal(nodos)
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
    A=construyeMatrizElipticaDiagonal(nodos)
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