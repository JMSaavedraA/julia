include("Solvers.jl")
function factorizaQR(A::AbstractMatrix,m::Int)
    Q=zeros(Float64,m,m)
    R=zeros(Float64,m,m)
    for i=1:m
        for j=1:i-1
            aAux=0.0
            for l=1:m
                aAux+=Q[l,j]*A[l,i]
            end
            R[j,i]=aAux
        end
        nAux=0.0
        for l=1:m
            aAux=A[l,i]
            for j=1:i-1
                aAux-=R[j,i]*Q[l,j]
            end
            Q[l,i]=aAux
            nAux+=aAux*aAux
        end
        nAux=sqrt(nAux)
        R[i,i]=nAux
        for l=1:m
            Q[l,i]=Q[l,i]/nAux
        end
    end
    return Q,R
end

function resuelveQR(A::AbstractMatrix,b::AbstractVector,m::Int)
    Q,R=factorizaQR(A,m)
    y=zeros(Float64,m)
    for i=1:m
        for j=1:m
            y[i]+=Q[j,i]*b[j]
        end
    end
    x=resuelveTriangularSuperior(R,y,m)
    return x
end