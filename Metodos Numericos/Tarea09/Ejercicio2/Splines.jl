include("metodosTridiagonales.jl")

function splineLineal(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int)
    #Spline de continuidad c=0
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        h[i]=xd-x[i]
        t[i]=y[i+1]-y[i]
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        yi=y[j]
        xi=x[j]
        di=t[j]/h[j]
        f[i]=yi+di*(z[i]-xi)
    end
    return f    
end
        
function splineCuadraticoDerecha(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int,d::AbstractFloat)
    #Spline de continuidad c=1 con derivada dada en x_n
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    H=zeros(Float64,n,3) #Matriz tridiagonal del sistema para S'
    b=zeros(Float64,n) #Lado derecho del sistema para S'
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        hi=xd-x[i]
        h[i]=hi
        H[i,2]=hi #Construimos la matriz del sistema
        H[i,3]=hi
        ti=y[i+1]-y[i]
        t[i]=ti
        b[i]=ti
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    H[n,2]=1 #En este caso añadimos el renglón final de la matriz del sistema, de forma que sea tridiagonal y triangular Superior
    b.=2*b
    b[n]=d #Añadimos el valor de la derivada en el último punto
    sPrime=resuelveTridiagonalSuperior(H,b,n)
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        si=sPrime[j]
        sd=sPrime[j+1]
        yi=y[j]
        xi=x[j]
        hi=h[j]
        f[i]=yi+si*(z[i]-xi)+(sd-si)/(2*hi)*(z[i]-xi)^2
    end
    return f
end
function splineCuadraticoDerecha(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int)
    #Spline de continuidad c=1 con derivada aproximada en x_n
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    H=zeros(Float64,n,3) #Matriz tridiagonal del sistema para S'
    b=zeros(Float64,n) #Lado derecho del sistema para S'
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        hi=xd-x[i]
        h[i]=hi
        H[i,2]=hi #Construimos la matriz del sistema
        H[i,3]=hi
        ti=y[i+1]-y[i]
        t[i]=ti
        b[i]=ti
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    H[n,2]=1 #En este caso añadimos el renglón final de la matriz del sistema, de forma que sea tridiagonal y triangular Superior
    b.=2*b
    b[n]=t[n-1]/h[n-1]  #Añadimos el valor aproximado de la derivada en el último punto
    sPrime=resuelveTridiagonalSuperior(H,b,n)
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        si=sPrime[j]
        sd=sPrime[j+1]
        yi=y[j]
        xi=x[j]
        hi=h[j]
        f[i]=yi+si*(z[i]-xi)+(sd-si)/(2*hi)*(z[i]-xi)^2
    end
    return f
end

function splineCuadraticoIzquierda(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int,d::AbstractFloat)
    #Spline de continuidad c=1 con derivada dada en x_1
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    H=zeros(Float64,n,3) #Matriz tridiagonal del sistema para S'
    b=zeros(Float64,n) #Lado derecho del sistema para S'
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        hi=xd-x[i]
        h[i]=hi
        H[i+1,1]=hi #Construimos la matriz del sistema
        H[i+1,2]=hi
        ti=y[i+1]-y[i]
        t[i]=ti
        b[i+1]=ti
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    H[1,2]=1 #En este caso añadimos el primer renglón de la matriz del sistema, de forma que sea tridiagonal y triangular inferior
    b.=2*b
    b[1]=d #Añadimos el valor de la derivada en el primer punto
    sPrime=resuelveTridiagonalInferior(H,b,n)
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        si=sPrime[j]
        sd=sPrime[j+1]
        yi=y[j]
        xi=x[j]
        hi=h[j]
        f[i]=yi+si*(z[i]-xi)+(sd-si)/(2*hi)*(z[i]-xi)^2
    end
    return f
end
function splineCuadraticoIzquierda(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int)
    #Spline de continuidad c=1 con derivada dada en x_1
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    H=zeros(Float64,n,3) #Matriz tridiagonal del sistema para S'
    b=zeros(Float64,n) #Lado derecho del sistema para S'
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        hi=xd-x[i]
        h[i]=hi
        H[i+1,1]=hi #Construimos la matriz del sistema
        H[i+1,2]=hi
        ti=y[i+1]-y[i]
        t[i]=ti
        b[i+1]=ti
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    H[1,2]=1 #En este caso añadimos el primer renglón de la matriz del sistema, de forma que sea tridiagonal y triangular inferior
    b.=2*b
    b[1]=t[1]/h[1] #Añadimos el valor aproximado de la derivada en el primer punto
    sPrime=resuelveTridiagonalInferior(H,b,n)
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        si=sPrime[j]
        sd=sPrime[j+1]
        yi=y[j]
        xi=x[j]
        hi=h[j]
        f[i]=yi+si*(z[i]-xi)+(sd-si)/(2*hi)*(z[i]-xi)^2
    end
    return f
end

function splineCubico(z::AbstractVector, x::AbstractVector, y::AbstractVector,m::Int,n::Int,ddI::AbstractFloat,ddD::AbstractFloat)
    #Spline de continuidad c=2 con la segunda derivada dada en x_1 y x_n
    I=zeros(typeof(m),m).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    h=zeros(Float64,n-1) #h_i
    H=zeros(Float64,n-2,3) #Matriz tridiagonal del sistema para S''
    b=zeros(Float64,n-2) #Lado derecho del sistema para S'
    t=copy(h) #t_i
    z=sort(z) #Ordenamos el vector z donde se evalúa el Spline
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    j=1
    zj=z[j]
    for i=1:n-1
        xd=x[i+1]
        h[i]=xd-x[i]
         t[i]=y[i+1]-y[i]
        while zj<=xd && j<=m
            I[j]=i
            j+=1
            if j<=m
                zj=z[j]
            end
        end
    end
    #Construimos la matriz tridiagonal y simétrica del sistema y el vector del lado derecho
    for i=1:n-2
        h0=h[i]
        h1=h[i+1]
        t0=t[i]/h0
        t1=t[i+1]/h1
        H[i,1]=h0
        H[i,2]=2*(h0+h1)
        H[i,3]=h1
        b[i]=6*(t1-t0)
    end
    H[1,1]=0
    H[n-2,3]=0
    b[1]=b[1]-h[1]*ddI
    b[n-2]=b[n-2]-h[n-1]*ddD

    sPrimePrime=zeros(Float64,n)
    sPrimePrime[2:n-1]=resuelveCholeskyTridiagonal(H,b,n-2) #Resolvemos el sistema por Cholesky
    sPrimePrime[1]=ddI #La segunda derivada en x_1 está dada
    sPrimePrime[n]=ddD #La segunda derivada en x_n está dada
    #Cálculo del Spline
    f=zeros(Float64,m)
    for i=1:m
        j=I[i]
        s0=sPrimePrime[j]
        s1=sPrimePrime[j+1]
        hi=h[j]
        y0=y[j]
        y1=y[j+1]
        z0=z[i]-x[j]
        z1=x[j+1]-z[i]
        f[i]=(s0*z1^3 + s1*z0^3)/(6*hi) + (y1/hi - s1*hi/6)*z0 + (y0/hi - s0*hi/6)*z1
    end
    return f
end