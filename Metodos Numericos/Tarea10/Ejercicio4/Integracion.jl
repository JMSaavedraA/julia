function integralMonteCarlo(f::Function,a::AbstractVector,b::AbstractVector, m::Int, k::Int,tol::AbstractFloat,n::Int,fMax::AbstractFloat,fMin::AbstractFloat)
    h=b-a #Ancho del intervalo
    L=fMax-fMin #Alto de la función
    A=L*prod(h) #Volumen del paralelepípedo D×L
    I=zeros(Float64,k)
    eAprox=10*tol
    i=0
    iOld=1.0
    R=0
    while i<k && eAprox>tol
        for j=1:m
            x=h.*rand(n)+a
            y=L*rand()
            fx=f(x)-fMin
            R+=(fx>y)
        end
        i+=1
        I[i]=R/(m*i) #Aproximación actual
        eAprox=abs((iOld*i/R)-1)
        iOld=R/i
    end
    res=I[i]*A + fMin*prod(h)
    return res
end

function integralMonteCarlo(f::Function,a::AbstractFloat,b::AbstractFloat, m::Int, k::Int,tol::AbstractFloat,fMax::AbstractFloat,fMin::AbstractFloat)
    h=b-a #Ancho del intervalo
    L=fMax-fMin #Alto de la función
    A=L*h #Área del rectángulo D×L
    I=zeros(Float64,k)
    eAprox=10*tol
    i=0
    iOld=1.0
    R=0
    while i<k && eAprox>tol
        for j=1:m
            x=h*rand()+a
            y=L*rand()
            fx=f(x)-fMin
            R+=(fx>y)
        end
        i+=1
        I[i]=R/(m*i) #Aproximación actual
        eAprox=abs((iOld*i/R)-1)
        iOld=R/i
    end
    res=I[i]*A + fMin*h
    return res
end