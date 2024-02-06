include("Solvers.jl")
include("metodosTridiagonales.jl")

function metodoPotencia(A::AbstractMatrix,xOld::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int)
    #El método de la Potencia para una matriz en general
    k=1;
    (nInf,i)=findmax(abs.(xOld))
    xOld=xOld/nInf                                          #Normalizamos el vector original (Norma Infinito)
    xNew=copy(xOld)                                         #Inicializamos el nuevo vector
    error=10.0*tol                                          #Inicializamos el error
    lambda=0.0                                              #Inicializamos lambda
    while error>tol && k<maxIter
        xNew=A*xOld                                         #Iteración del Método de la Potencia
        lambda=xNew[i]/xOld[i]                              #Actualizar lambda
        (nInf,i)=findmax(abs.(xNew))    
        error=maximum(abs.(xNew - lambda*xOld))/nInf        #Calculamos el error
        xNew=xNew/nInf                                      #Normalizamos el nuevo vector (Norma Infinito)
        xOld.=xNew                                          #Actualizamos el vector anterior
        k=k+1   
    end
    n2=norma2(xNew,m)                                       #Regresamos el vector normalizado (Norma 2)
    xNew=xNew/n2
    return xNew,lambda
end

function metodoPotenciaTridiagonal(A::AbstractMatrix,xOld::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int)
    #El método de la Potencia para una matriz tridiagonal
    k=0;
    (nInf,i)=findmax(abs.(xOld))
    xOld=xOld/nInf                                          #Normalizamos el vector original (Norma Infinito)
    xNew=copy(xOld)                                         #Inicializamos el nuevo vector
    error=10.0*tol                                          #Inicializamos el error
    lambda=0.0                                              #Inicializamos lambda
    while error>tol && k<maxIter
        xNew[1]=A[1,2]*xOld[1] + A[1,3]*xOld[2]             #Iteración del Método de la Potencia
        for k=2:m-1
            xNew[k]=A[k,1]*xOld[k-1] + A[k,2]*xOld[k] + A[k,3]*xOld[k+1]
        end
        xNew[m]=A[m,1]*xOld[m-1] + A[m,2]*xOld[m]
        lambda=xNew[i]/xOld[i]                              #Actualizar lambda
        (nInf,i)=findmax(abs.(xNew))
        error=maximum(abs.(xNew - lambda*xOld))/nInf        #Calculamos el error
        xNew=xNew/nInf                                      #Normalizamos el nuevo vector (Norma Infinito)
        xOld.=xNew                                          #Actualizamos el vector anterior
        k=k+1
    end
    n2=norma2(xNew,m)                                       #Regresamos el vector normalizado (Norma 2)
    xNew=xNew/n2
    return xNew,lambda
end

function metodoPotenciaInversa(A::AbstractMatrix,xOld::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int)
    #El método de la Potencia Inversa para una matriz en general
    k=1;
    L=CholeskyLDL(A,m)                                  #Hacemos factorización de Cholesky
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,i]
        L[i,i]=1.0
    end
    (nInf,i)=findmax(abs.(xOld))
    xOld=xOld/nInf                                      #Normalizamos el vector original (Norma Infinito)
    xNew=copy(xOld)                                     #Inicializamos el nuevo vector
    error=10.0*tol                                      #Inicializamos el error
    lambda=0.0                                          #Inicializamos lambda
    while error>tol && k<maxIter
        y=resuelveTriangularInferior(L,xOld,m)          #Iteración del Método de la Potencia Inversa
        for i=1:m
            y[i]=y[i]/d[i]
        end
        xNew=resuelveTriangularSuperior(L,y,m)          #Resolviendo para LDL'x1=x0
        (nInf,i)=findmax(abs.(xNew))
        lambda=xOld[i]/xNew[i]                          #Actualizar lambda
        error=maximum(abs.(xNew - xOld/lambda))/nInf    #Calculamos el error
        xNew=xNew/nInf                                  #Normalizamos el nuevo vector (Norma Infinito)
        xOld.=xNew                                      #Actualizamos el vector anterior
        k=k+1
    end
    n2=norma2(xNew,m)                                   #Regresamos el vector normalizado (Norma 2)
    xNew=xNew/n2
    return xNew,lambda
end

function metodoPotenciaInversaTridiagonal(A::AbstractMatrix,xOld::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int)
    #El método de la Potencia Inversa para una matriz tridiagonal
    k=1;
    L=CholeskyLDLTridiagonal(A,m)                       #Hacemos factorización de Cholesky
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    (nInf,i)=findmax(abs.(xOld))
    xOld=xOld/nInf                                      #Normalizamos el vector original (Norma Infinito)
    xNew=copy(xOld)                                     #Inicializamos el nuevo vector
    error=10.0*tol                                      #Inicializamos el error
    lambda=0.0                                          #Inicializamos lambda
    while error>tol && k<maxIter
        y=resuelveTridiagonalInferior(L,xOld,m)         #Iteración del Método de la Potencia Inversa
        for i=1:m
            y[i]=y[i]/d[i]
        end
        xNew=resuelveTridiagonalSuperior(L,y,m)         #Resolviendo para LDL'x1=x0
        (nInf,i)=findmax(abs.(xNew))
        lambda=xOld[i]/xNew[i]                          #Actualizar lambda
        error=maximum(abs.(xNew - xOld/lambda))/nInf    #Calculamos el error
        xNew=xNew/nInf                                  #Normalizamos el nuevo vector (Norma Infinito)
        xOld.=xNew                                      #Actualizamos el vector anterior
        k=k+1
    end
    n2=norma2(xNew,m)
    xNew=xNew/n2                                        #Regresamos el vector normalizado (Norma 2)
    return xNew,lambda
end

function metodoPotenciaN(A::AbstractMatrix,x0::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int,n::Int)
    #El método de la Potencia generalizado para los n eigenpares más grandes de una matriz en general
    P=zeros(Float64,m,n)                                        #Inicializamos la matriz de eigenvectores
    L=zeros(Float64,n)                                          #Inicializamos el vector de eigenvalores
    xOld=similar(x0)                                            #Inicializamos el vector anterior
    xNew=similar(x0)                                            #Inicializamos el vector nuevo
    h=0                                                         #h nos ayudará si llegamos a algunos casos degenerados
    for s=1:n
        xOld.=x0                                                #Iniciamos en x0
        k=1;
        (nInf,i)=findmax(abs.(xOld))
        xOld=xOld/nInf                                          #Normalizamos el vector original (Norma Infinito)
        xNew.=xOld                                              #Inicializamos el nuevo vector
        error=10.0*tol                                          #Inicializamos el error
        lambda=0.0                                              #Inicializamos lambda
        while error>tol && k<maxIter
            xNew=A*xOld                                         #Iteración del Método de la Potencia
            lambda=xNew[i]/xOld[i]                              #Actualizar lambda
            (nInf,i)=findmax(abs.(xNew))
            error=maximum(abs.(xNew - lambda*xOld))/nInf        #Calculamos el error
            xNew=xNew/nInf                                      #Normalizamos el nuevo vector (Norma Infinito)
            xOld.=xNew                                          #Actualizamos el vector anterior
            j=1
            k=k+1
            if (k==maxIter && error>tol && h<m)                 #Caso degenerado en que no convergemos (posiblemente lambda_s=-lambda_{s+1})
                h+=1
                x0=zeros(Float64,m)                             #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                x0[h]=1
                xOld.=x0
                i=h
                j=1
                k=1                                             #Volvemos a empezar con este nuevo vector
            end
            while j<s
                a=0.0                                           #Vamos a quitar el error numérico de los eigenvectores que ya tenemos
                for l=1:m
                    a+=P[l,j]*xOld[l]                           #Cálculo de la proyección
                end
                for l=1:m
                    xOld[l]=xOld[l]-a*P[l,j]                    #Quitamos la proyección
                end
                nOld=norma2(xOld,m)                             #Cálculo de la norma 2 del vector sin la proyección
                if nOld<tol                                     #Caso degenerado en que el vector era combinación lineal de los eigenvectores anteriores
                    h+=1
                    x0=zeros(Float64,m)                         #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                    x0[h]=1
                    xOld.=x0
                    i=h
                    j=1
                    k=1                                         #Volvemos a empezar con este nuevo vector
                    error=10.0*tol
                else                                            #No estamos en casos degenerados
                    xOld=xOld/nOld                              #Normalizamos el vector sin la proyección
                    j=j+1
                end
            end
        end
        if xOld[1]<-tol
            xOld=-xOld                                          #Por fines comparativos, ponemos el vector de forma que su primera entrada sea no negativa
        end
        z=0                                                     #Vamos a ordenar el nuevo eigenpar (puede no ser el más pequeño si x0 no tiene proyección en todos)
        continuar=true
        while z<s && continuar
            z+=1
            if L[z]>lambda
                continuar=false
            end
        end
        u=maximum([z,2])
        for w=s:-1:u
            L[w]=L[w-1]
            for l=1:m
                P[l,w]=P[l,w-1]
            end
        end
        L[z]=lambda
        for l=1:m
            P[l,z]=xOld[l]
        end
    end
    return P,L
end

function metodoPotenciaNtridiagonal(A::AbstractMatrix,x0::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int,n::Int)
    #El método de la Potencia generalizado para los n eigenpares más grandes de una matriz tridiagonal
    P=zeros(Float64,m,n)                                        #Inicializamos la matriz de eigenvectores
    L=zeros(Float64,n)                                          #Inicializamos el vector de eigenvalores
    h=0                                                         #h nos ayudará si llegamos a algunos casos degenerados
    xOld=similar(x0)                                            #Inicializamos el vector anterior
    xNew=similar(x0)                                            #Inicializamos el vector nuevo
    for s=1:n
        xOld.=x0                                                #Iniciamos en x0
        k=1;
        (nInf,i)=findmax(abs.(xOld))
        xOld=xOld/nInf                                          #Normalizamos el vector original (Norma Infinito)
        xNew.=xOld                                              #Inicializamos el nuevo vector
        error=10.0*tol                                          #Inicializamos el error
        lambda=0.0                                              #Inicializamos lambda
        while error>tol && k<maxIter
            xNew[1]=A[1,2]*xOld[1] + A[1,3]*xOld[2]             #Iteración del Método de la Potencia
            for l=2:m-1
                xNew[l]=A[l,1]*xOld[l-1] + A[l,2]*xOld[l] + A[l,3]*xOld[l+1]
            end
            xNew[m]=A[m,1]*xOld[m-1] + A[m,2]*xOld[m]
            lambda=xNew[i]/xOld[i]                              #Actualizar lambda
            (nInf,i)=findmax(abs.(xNew))
            error=maximum(abs.(xNew - lambda*xOld))/nInf        #Calculamos el error
            xNew=xNew/nInf                                      #Normalizamos el nuevo vector (Norma Infinito)
            xOld.=xNew                                          #Actualizamos el vector anterior
            j=1
            k=k+1
            if (k==maxIter && error>tol && h<m)                 #Caso degenerado en que no convergemos (posiblemente lambda_s=-lambda_{s+1})
                h+=1
                x0=zeros(Float64,m)                             #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                x0[h]=1
                xOld.=x0
                i=h
                j=1
                k=1                                             #Volvemos a empezar con este nuevo vector
            end
            while j<s
                a=0.0                                           #Vamos a quitar el error numérico de los eigenvectores que ya tenemos
                for l=1:m
                    a+=P[l,j]*xOld[l]                           #Cálculo de la proyección
                end
                for l=1:m
                    xOld[l]=xOld[l]-a*P[l,j]                    #Quitamos la proyección
                end
                nOld=norma2(xOld,m)                             #Cálculo de la norma 2 del vector sin la proyección
                if nOld<tol                                     #Caso degenerado en que el vector era combinación lineal de los eigenvectores anteriores
                    h+=1
                    x0=zeros(Float64,m)                         #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                    x0[h]=1
                    xOld.=x0
                    i=h
                    j=1
                    k=1                                         #Volvemos a empezar con este nuevo vector
                    error=10.0*tol
                else                                            #No estamos en casos degenerados
                    xOld=xOld/nOld                              #Normalizamos el vector sin la proyección
                    j=j+1
                end
            end
        end
        if xOld[1]<-tol                                         #Por fines comparativos, ponemos el vector de forma que su primera entrada sea no negativa
            xOld=-xOld
        end
        z=0
        continuar=true
        while z<s && continuar                                  #Vamos a ordenar el nuevo eigenpar (puede no ser el más pequeño si x0 no tiene proyección en todos)
            z+=1
            if L[z]>lambda
                continuar=false
            end
        end
        u=maximum([z,2])
        for w=s:-1:u
            L[w]=L[w-1]
            for l=1:m
                P[l,w]=P[l,w-1]
            end
        end
        L[z]=lambda
        for l=1:m
            P[l,z]=xOld[l]
        end
    end
    return P,L
end

function metodoPotenciaInversaN(A::AbstractMatrix,x0::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int,n::Int)
    #El método de la Potencia Inversa generalizado para los n eigenpares más pequeños de una matriz en general
    P=zeros(Float64,m,n)                                    #Inicializamos la matriz de eigenvectores
    D=zeros(Float64,n)                                      #Inicializamos el vector de eigenvalores
    L=CholeskyLDL(A,m)                                      #Hacemos factorización de Cholesky
    h=0                                                     #h nos ayudará si llegamos a algunos casos degenerados
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,i]
        L[i,i]=1.0
    end
    Lt=transpose(L)
    xOld=similar(x0)                                        #Inicializamos el vector anterior
    xNew=similar(x0)                                        #Inicializamos el vector nuevo
    for s=1:n
        xOld.=x0                                            #Iniciamos en x0
        k=1;
        (nInf,i)=findmax(abs.(xOld))
        xOld=xOld/nInf                                      #Normalizamos el vector original (Norma Infinito)
        xNew.=xOld                                          #Inicializamos el nuevo vector
        error=10.0*tol                                      #Inicializamos el error
        lambda=0.0                                          #Inicializamos lambda
        while error>tol && k<maxIter
            y=resuelveTriangularInferior(L,xOld,m)          #Iteración del Método de la Potencia Inversa
            for i=1:m
                y[i]=y[i]/d[i]
            end
            xNew=resuelveTriangularSuperior(Lt,y,m)         #Resolviendo para LDL'x1=x0
            (nInf,i)=findmax(abs.(xNew))
            lambda=xOld[i]/xNew[i]                          #Actualizar lambda
            error=maximum(abs.(xNew - xOld/lambda))/nInf    #Calculamos el error
            xNew=xNew/nInf                                  #Normalizamos el nuevo vector (Norma Infinito)
            xOld.=xNew                                      #Actualizamos el vector anterior
            j=1
            k=k+1
            if (k==maxIter && error>tol && h<m)             #Caso degenerado en que no convergemos (posiblemente lambda_s=-lambda_{s+1})
                h+=1
                x0=zeros(Float64,m)                         #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                x0[h]=1
                xOld.=x0
                i=h
                j=1
                k=1                                         #Volvemos a empezar con este nuevo vector
            end
            while j<s
                a=0.0                                       #Vamos a quitar el error numérico de los eigenvectores que ya tenemos
                for l=1:m
                    a+=P[l,j]*xOld[l]                       #Cálculo de la proyección
                end
                for l=1:m
                    xOld[l]=xOld[l]-a*P[l,j]                #Quitamos la proyección
                end
                nOld=norma2(xOld,m)                         #Cálculo de la norma 2 del vector sin la proyección
                if nOld<tol                                 #Caso degenerado en que el vector era combinación lineal de los eigenvectores anteriores
                    h+=1
                    x0=zeros(Float64,m)                     #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                    x0[h]=1
                    xOld.=x0
                    i=h
                    j=1
                    k=1                                     #Volvemos a empezar con este nuevo vector
                    error=10.0*tol
                else                                        #No estamos en casos degenerados
                    xOld=xOld/nOld                          #Normalizamos el vector sin la proyección
                    j=j+1
                end
            end
        end
        if xOld[1]<-tol
            xOld=-xOld                                      #Por fines comparativos, ponemos el vector de forma que su primera entrada sea no negativa
        end
        z=0
        continuar=true
        while z<s && continuar                              #Vamos a ordenar el nuevo eigenpar (puede no ser el más grande si x0 no tiene proyección en todos)
            z+=1
            if D[z]<lambda
                continuar=false
            end
        end
        u=maximum([z,2])
        for w=s:-1:u
            D[w]=D[w-1]
            for l=1:m
                P[l,w]=P[l,w-1]
            end
        end
        D[z]=lambda
        for l=1:m
            P[l,z]=xOld[l]
        end
    end
    return P,D
end

function metodoPotenciaInversaNtridiagonal(A::AbstractMatrix,x0::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int,n::Int)
    #El método de la Potencia Inversa generalizado para los n eigenpares más pequeños de una matriz tridiagonal
    P=zeros(Float64,m,n)                                        #Inicializamos la matriz de eigenvectores
    D=zeros(Float64,n)                                          #Inicializamos el vector de eigenvalores
    L=CholeskyLDLTridiagonal(A,m)                               #Hacemos factorización de Cholesky
    h=0                                                         #h nos ayudará si llegamos a algunos casos degenerados
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    xOld=similar(x0)                                            #Inicializamos el vector anterior
    xNew=similar(x0)                                            #Inicializamos el vector nuevo
    for s=1:n
        xOld.=x0                                                #Iniciamos en x0
        k=1;
        (nInf,i)=findmax(abs.(xOld))
        xOld=xOld/nInf                                          #Normalizamos el vector original (Norma Infinito)
        xNew.=xOld                                              #Inicializamos el nuevo vector
        error=10.0*tol                                          #Inicializamos el error
        lambda=0.0                                              #Inicializamos lambda
        while error>tol && k<maxIter
            y=resuelveTridiagonalInferior(L,xOld,m)             #Iteración del Método de la Potencia Inversa
            for i=1:m
                y[i]=y[i]/d[i]
            end
            xNew=resuelveTridiagonalSuperior(L,y,m)             #Resolviendo para LDL'x1=x0
            (nInf,i)=findmax(abs.(xNew))
            lambda=xOld[i]/xNew[i]                              #Actualizar lambda
            error=maximum(abs.(xNew - xOld/lambda))/nInf        #Calculamos el error
            xNew=xNew/nInf                                      #Normalizamos el nuevo vector (Norma Infinito)
            xOld.=xNew                                          #Actualizamos el vector anterior
            j=1
            k=k+1
            if (k==maxIter && error>tol && h<m)                 #Caso degenerado en que no convergemos (posiblemente lambda_s=-lambda_{s+1})
                h+=1
                x0=zeros(Float64,m)                             #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                x0[h]=1
                xOld.=x0
                i=h
                j=1
                k=1                                             #Volvemos a empezar con este nuevo vector
            end
            while j<s
                a=0.0                                           #Vamos a quitar el error numérico de los eigenvectores que ya tenemos
                for l=1:m
                    a+=P[l,j]*xOld[l]                           #Cálculo de la proyección
                end
                for l=1:m
                    xOld[l]=xOld[l]-a*P[l,j]                    #Quitamos la proyección
                end
                nOld=norma2(xOld,m)                             #Cálculo de la norma 2 del vector sin la proyección
                if nOld<tol                                     #Caso degenerado en que el vector era combinación lineal de los eigenvectores anteriores
                    h+=1
                    x0=zeros(Float64,m)                         #Si estamos en casos degenerados, ponemos un vector que no hayamos usado
                    x0[h]=1
                    xOld.=x0
                    i=h
                    j=1
                    k=1                                         #Volvemos a empezar con este nuevo vector
                    error=10.0*tol
                else                                            #No estamos en casos degenerados
                    xOld=xOld/nOld                              #Normalizamos el vector sin la proyección
                    j=j+1
                end
            end
        end
        if xOld[1]<-tol                                         #Por fines comparativos, ponemos el vector de forma que su primera entrada sea no negativa
            xOld=-xOld
        end
        z=0
        continuar=true
        while z<s && continuar                                  #Vamos a ordenar el nuevo eigenpar (puede no ser el más grande si x0 no tiene proyección en todos)
            z+=1
            if D[z]<lambda
                continuar=false
            end
        end
        u=maximum([z,2])
        for w=s:-1:u
            D[w]=D[w-1]
            for l=1:m
                P[l,w]=P[l,w-1]
            end
        end
        D[z]=lambda
        for l=1:m
            P[l,z]=xOld[l]
        end
    end
    return P,D
end