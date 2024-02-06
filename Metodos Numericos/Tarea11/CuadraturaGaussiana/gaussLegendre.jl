function cuadraturaGaussiana(f::Function,a::Number,b::Number,n::Int)
    W=[2.0 0.0 0.0;1.0 1.0 0.0;5/9 8/9 5/9]
    U=[0.0 0.0 0.0;-1/sqrt(3) 1/sqrt(3) 0;-sqrt(6/10) 0.0 sqrt(6/10)]
    d=(b-a)/2
    I=0.0
    for i=1:n
        I+=d*W[n,i]*f(d*(U[n,i]+1)+a)
    end
    return I
end

function cuadraturaGaussiana(f::Function,a::Number,b::Number,n::Int,m::Int)
    W=[2.0 0.0 0.0;1.0 1.0 0.0;5/9 8/9 5/9]
    U=[0.0 0.0 0.0;-1/sqrt(3) 1/sqrt(3) 0;-sqrt(6/10) 0.0 sqrt(6/10)]
    h=(b-a)/(2*m)
    I=0.0
    for i=1:m
        for j=1:n
            I+=h*W[n,j]*f(h*(U[n,j]+2*i-1)+a)
        end
    end
    return I
end