include("newtonCotes.jl")

function integralRomberg(f::Function,a::Number,b::Number,n::Int)
    T=zeros(Float64,n+1,n+1)
    
    for i=1:n+1
        T[i,1]=integralTrapecio(f,a,b,2^(i-1))
        for j=2:i
            c=4^(j-1)
            T[i,j]=(c*T[i,j-1]-T[i-1,j-1])/(c-1)
        end
    end
    return T[n+1,n+1]
end

function integralRichardson(f::Function,metodoIntegral::Function,a::Number,b::Number,n::Int)
    A=zeros(Float64,n+1,n+1)
    
    for i=1:n+1
        A[i,1]=metodoIntegral(f,a,b,2^(i-1))
        for j=2:i
            c=2^(j-1)
            A[i,j]=(c*A[i,j-1]-A[i-1,j-1])/(c-1)
        end
    end
    return A[n+1,n+1]
end