function integralTrapecio(f::Function,a::Number,b::Number,n::Int)
    x=LinRange(a,b,n+1)
    h=x[2]-x[1]
    fx=h*f.(x)
    I=sum(fx)
    I-=fx[1]/2
    I-=fx[n+1]/2
    return I
end

function integralSimpson(f::Function,a::Number,b::Number,n::Int)
    x=LinRange(a,b,2*n+1)
    h=2*(x[2]-x[1])/3
    fx=h*f.(x)
    I=fx[1]/2
    for i=1:n
        I+=2*fx[2*i]
        I+=fx[2*i+1]
    end
    I-=fx[2*n+1]/2
    return I
end

function integralTresOctavos(f::Function,a::Number,b::Number,n::Int)
    x=LinRange(a,b,3*n+1)
    h=3*(x[2]-x[1])/8
    fx=h*f.(x)
    I=fx[1]
    for i=1:n
        I+=3*fx[3*i-1]
        I+=3*fx[3*i]
        I+=2*fx[3*i+1]
    end
    I-=fx[3*n+1]
    return I
end

function integralMilne(f::Function,a::Number,b::Number,n::Int)
    x=LinRange(a,b,4*n+1)
    h=4*(x[2]-x[1])/45
    fx=h*f.(x)
    I=7*fx[1]/2
    for i=1:n
        I+=16*fx[4*i-2]
        I+=6*fx[4*i-1]
        I+=16*fx[4*i]
        I+=7*fx[4*i+1]
    end
    I-=7*fx[4*n+1]/2
    return I
end

function integralTrapecio(f::Function,a::Number,b::Number,n::Int,nMax::Int,tol::AbstractFloat)
    oldI=integralTrapecio(f,a,b,n)
    n+=1
    newI=integralTrapecio(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    while eAbs>tol && n<nMax
    newI=integralTrapecio(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    end
    return newI
end

function integralSimpson(f::Function,a::Number,b::Number,n::Int,nMax::Int,tol::AbstractFloat)
    oldI=integralSimpson(f,a,b,n)
    n+=1
    newI=integralSimpson(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    while eAbs>tol && n<nMax
    newI=integralSimpson(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    end
    return newI
end

function integralTresOctavos(f::Function,a::Number,b::Number,n::Int,nMax::Int,tol::AbstractFloat)
    oldI=integralTresOctavos(f,a,b,n)
    n+=1
    newI=integralTresOctavos(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    while eAbs>tol && n<nMax
    newI=integralTresOctavos(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    end
    return newI
end

function integralMilne(f::Function,a::Number,b::Number,n::Int,nMax::Int,tol::AbstractFloat)
    oldI=integralMilne(f,a,b,n)
    n+=1
    newI=integralMilne(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    while eAbs>tol && n<nMax
    newI=integralMilne(f,a,b,n)
    eAbs=abs(newI-oldI)
    oldI=newI
    n+=1
    end
    return newI
end