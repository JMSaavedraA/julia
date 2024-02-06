include("newtonCotes.jl")

function integralRomberg(f::Function,a::Number,b::Number,n::Int)
    T=zeros(Float64,n,n)
    
    for i=1:n
        T[i,1]=integralTrapecio(f,a,b,2^(i-1))
        for j=2:i
            c=4^(j-1)
            T[i,j]=(c*T[i,j-1]-T[i-1,j-1])/(c-1)
        end
    end
    return T[n,n]
end

function integralRichardson(f::Function,metodoIntegral::Function,a::Number,b::Number,n::Int)
    A=zeros(Float64,n,n)
    
    for i=1:n
        A[i,1]=metodoIntegral(f,a,b,2^(i-1))
        for j=2:i
            c=2^(j-1)
            A[i,j]=(c*A[i,j-1]-A[i-1,j-1])/(c-1)
        end
    end
    return A[n,n]
end