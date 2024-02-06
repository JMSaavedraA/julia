include("metodoJacobi.jl")
include("Potencia.jl")
include("factQR.jl")
include("gradienteConjugado.jl")


function iteracionSubespaciosPotencia(A::AbstractMatrix,m::Int,n::Int, maxIter::Int, tol::AbstractFloat, iterJacobi::Int, iterPotencia::Int,x0::AbstractMatrix)
    I=copy(x0)
    iter=0
    Phi,_=metodoPotenciaN(A,I,m,tol,iterPotencia,n)
    B=Phi'*A*Phi
    zOld=Phi[:,1]
    zNew=copy(zOld)
    y=10*tol
    while iter<maxIter && y>tol
        Z,_,_,_=eigenJacobi(B,n,tol,iterJacobi)
        I.=A*Phi*Z
        factorizaQR(I,m,n,Phi)
        B.=Phi'*A*Phi
        zNew.=Phi[:,1]
        nNew=norma2(zNew,m)
        y=norma2(zNew-zOld,m)/nNew
        zOld.=zNew
        iter+=1
    end
    d=vectorDiagonal(B,n)
    return Phi,d
end

function iteracionSubespaciosPotenciaInversa(A::AbstractMatrix,m::Int,n::Int, maxIter::Int, tol::AbstractFloat, iterJacobi::Int, iterPotencia::Int,x0::AbstractMatrix)
    L=CholeskyLDL(A,m)
    D=zeros(Float64,m)
    for i=1:m
        D[i]=L[i,i]
        L[i,i]=1.0
    end
    I=copy(x0)
    iter=0
    Phi,_=metodoPotenciaInversaNcholesky(L,D,I,m,tol,iterJacobi,n)
    B=Phi'*A*Phi
    zOld=Phi[:,1]
    zNew=copy(zOld)
    b=copy(zOld)
    y=10*tol
    while iter<maxIter && y>tol
        Z,_,_,_=eigenJacobi(B,n,tol,iterJacobi)
        I.=Phi*Z
        for i=1:n
            b.=I[:,i]
            I[:,i]=gradienteConjugado(A,b,m,tol,3*m)
        end
        factorizaQR(I,m,n,Phi)
        B.=Phi'*A*Phi
        zNew.=Phi[:,1]
        nNew=norma2(zNew,m)
        y=norma2(zNew-zOld,m)/nNew
        zOld.=zNew
        iter+=1
    end
    d=vectorDiagonal(B,n)
    return Phi,d
end