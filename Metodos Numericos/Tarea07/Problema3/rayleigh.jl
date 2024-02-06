include("Solvers.jl")

function eigenRayleigh(A::AbstractMatrix,vAprox::AbstractVector,lambdaAprox::AbstractFloat,m::Int,tol::AbstractFloat,maxIter::Int)
    vOld=copy(vAprox)
    nNew=norma2(vOld,m)
    vOld=vOld/nNew
    vNew=copy(vOld)
    lambdaOld=lambdaAprox
    lambdaNew=lambdaOld
    exact=10*tol
    B=copy(A)
    for j=1:m
        B[j,j]-=lambdaOld
    end
    i=0
    while i<maxIter && exact>tol
        vNew.=resuelveGaussJordan(B,vOld,m)
        aAux=0.0
        nNew=0.0
        for j=1:m
            aAux+=vNew[j]*vOld[j]
            nNew+=vNew[j]*vNew[j]
        end
        lambdaNew+=1/aAux
        exact=abs(lambdaOld-lambdaNew)
        for j=1:m
            B[j,j]=B[j,j] +lambdaOld-lambdaNew
        end
        lambdaOld=lambdaNew
        vOld.=vNew/sqrt(nNew)
        i+=1
    end
    return vOld,lambdaOld
end