include("Solvers.jl")

function coeficientesMinimosCuadrados(F::Function,x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #F es una funcion R->R^n con las que se interpola, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de x, n es el numero de funciones con que se interpola
    Y=zeros(m,n)
    for i=1:m
        Y[i,:]=F(x[i])
    end
    A=zeros(Float64,n,n)
    b=zeros(Float64,n)
    multiplicaMatrizVector(Y',y,b,n,m)
    multiplicaMatrices(Y',Y,A,n,m,n)
    #gradienteConjugado(A,b,n,1e-10,2*n)
    a=resuelveCholesky(A,b,n)
    return a
end

function coeficientesMinimosCuadradosPolinomios(x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de x, n es el orden del polinomio
    Y=zeros(m,n+1)
    for i=1:m
        Y[i,1]=1
        for j=1:n
            Y[i,j+1]=x[i]^j
        end
    end
    A=zeros(Float64,n+1,n+1)
    b=zeros(Float64,n+1)
    multiplicaMatrizVector(Y',y,b,n+1,m)
    multiplicaMatrices(Y',Y,A,n+1,m,n+1)
    a=resuelveCholesky(A,b,n+1)
    return a
end

function minimosCuadrados(F::Function,x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #F es una funcion R->R^n con las que se interpola, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de x, n es el numero de funciones con que se interpola
    a=coeficientesMinimosCuadrados(F,x,y,m,n)
    p=zeros(Float64,m)
    for i=1:m
        p[i]=a'*F(x[i])
    end
    return p
end

function minimosCuadradosPolinomios(x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de x, n es el orden del polinomio
    a=coeficientesMinimosCuadradosPolinomios(x,y,m,n)
    p=zeros(Float64,m)
    for i=1:m
        pi=a[1]
        for j=1:n
            pi+=a[j+1]*(x[i]^j)
        end
        p[i]=pi
    end
    return p
end

function interpolarPolinomios(z::AbstractVector,x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de x, n es el orden del polinomio
    a=coeficientesMinimosCuadradosPolinomios(x,y,n,n-1)
    p=zeros(Float64,m)
    for i=1:m
        pi=a[1]
        for j=1:n-1
            pi+=a[j+1]*(z[i]^j)
        end
        p[i]=pi
    end
    return p
end

function polinomioLagrange(z::AbstractFloat,x::AbstractVector,y::AbstractVector,n::Int)
    #z es el punto donde se evalua el polinomio, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de z, n es el tamaño de x
    p=0.0
    for i=1:n
        li=1.0
        for j=1:i-1
            li=li*(z-x[j])/(x[i]-x[j])
        end
        for j=i+1:n
            li=li*(z-x[j])/(x[i]-x[j])
        end
        p+=y[i]*li
    end
    return p
end

function polinomioLagrange(z::AbstractVector,x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #z es el vector donde se evalua el polinomio, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de z, n es el tamaño de x
    p=zeros(Float64,m)
    for i=1:n
        li=zeros(Float64,m).+1.0
        for j=1:i-1
            li=li.*((z.-x[j])./(x[i]-x[j]))
        end
        for j=i+1:n
            li=li.*((z.-x[j])./(x[i]-x[j]))
        end
        p+=y[i]*li
    end
    return p
end

function diferenciaDividida(d::AbstractVector,x::AbstractVector,n::Int,k::Int)
    #Se evalua una diferencia dividida para el polinomio de Newton.
    #d es la diferencia dividida anterior, x es el vector de puntos interpolados, n es el tamaño de x, k es el orden de la diferencia
    D=zeros(Float64,n-k)
    for i=1:n-k
        D[i]=(d[i+1]-d[i])/(x[i+k]-x[i])
    end
    return D
end

function polinomioNewton(z::AbstractVector,x::AbstractVector,y::AbstractVector,m::Int,n::Int)
    #z es el vector donde se evalua el polinomio, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de z, n es el tamaño de x
    D=zeros(Float64,n,n)
    D[:,1]=y
    for j=1:n-1
        d=D[1:n+1-j,j]
        dr=diferenciaDividida(d,x,n,j)
        D[1:n-j,j+1]=dr
    end
    p=zeros(Float64,m).+D[1,1]
    for k=1:n-1
        pk=zeros(Float64,m).+D[1,k+1]
        for i=1:k
            pk=pk.*(z.-x[i])
        end
        p.=p+pk
    end
    return p
end

function polinomioNewton(z::AbstractFloat,x::AbstractVector,y::AbstractVector,n::Int)
    #z es el punto donde se evalua el polinomio, x es el vector de puntos donde se interpola, y es f(x), m es el tamaño de z, n es el tamaño de x
    D=zeros(Float64,n,n)
    D[:,1]=y
    for j=1:n-1
        d=D[1:n+1-j,j]
        dr=diferenciaDividida(d,x,n,j)
        D[1:n-j,j+1]=dr
    end
    p=D[1,1]
    for k=1:n-1
        pk=D[1,k+1]
        for i=1:k
            pk=pk*(z-x[i])
        end
        p=p+pk
    end
    return p
end