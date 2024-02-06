include("Solvers.jl")

function gradienteConjugado(A::AbstractMatrix,b::AbstractVector,m::Int,tol::AbstractFloat,maxIter::Int)
    r=-copy(b)
    p=copy(b)
    w=zeros(Float64,m)
    x=zeros(Float64,m)
    k=0
    nNew=norma2(r,m)
    sNew=nNew^2
    while k<maxIter && nNew>tol
        sOld=sNew
        for i=1:m
            w[i]=0.0
            for j=1:m
                w[i]=w[i]+A[i,j]p[j]
            end
        end
        denom=0.0
        for i=1:m
            denom+=p[i]*w[i]
        end
        alfa=sOld/denom
        for i=1:m
            x[i]=x[i]+alfa*p[i]
        end
        for i=1:m
            r[i]=r[i]+alfa*w[i]
        end
        sNew=0.0
        for i=1:m
            sNew+=r[i]*r[i]
        end
        beta=sNew/sOld
        for i=1:m
            p[i]=-r[i]+beta*p[i]
        end
        nNew=sqrt(sNew)
        k+=1
    end
    return x
end