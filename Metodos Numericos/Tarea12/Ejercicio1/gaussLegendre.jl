function cuadraturaGaussiana(f::Function,a::Number,b::Number,n::Int)
    #cuadratura de Gauss-Legendre para integrar f en [a,b] por un n puntos
    W=[2.0 0.0 0.0;1.0 1.0 0.0;5/9 8/9 5/9]
    U=[0.0 0.0 0.0;-1/sqrt(3) 1/sqrt(3) 0;-sqrt(6/10) 0.0 sqrt(6/10)]
    d=(b-a)/2 #tamaño del intervalo
    I=0.0 #Inicialización de la integral
    for i=1:n #ciclo sobre los puntos
        I+=d*W[n,i]*f(d*(U[n,i]+1)+a)
    end
    return I
end

function cuadraturaGaussiana(f::Function,a::Number,b::Number,n::Int,m::Int)
    #cuadratura de Gauss-Legendre para integrar f en [a,b] por un n puntos en m+1 intervalos equiespaciados
    W=[2.0 0.0 0.0;1.0 1.0 0.0;5/9 8/9 5/9]
    U=[0.0 0.0 0.0;-1/sqrt(3) 1/sqrt(3) 0;-sqrt(6/10) 0.0 sqrt(6/10)]
    h=(b-a)/(2*m)#tamaño del intervalo
    I=0.0 #Inicialización de la integral
    for i=1:m #ciclo sobre los intervalos
        for j=1:n #ciclo sobre los puntos
            I+=h*W[n,j]*f(h*(U[n,j]+2*i-1)+a)
        end
    end
    return I
end