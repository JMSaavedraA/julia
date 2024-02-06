include("metodosTridiagonales.jl")
include("problemaExamen.jl")
include("Solvers.jl")
using Printf

m=5000
A=crearMatrizExamenPentadiagonal(m)
b=crearVectorExamen(m)
L=CholeskyLDLPentadiagonal(A,m)
d=similar(b)
for i=1:m
    d[i]=L[i,3]
    L[i,3]=1.0
end
@time x=resuelveCholeskyPentadiagonalDada(L,d,b,m)
bEstimado=multiplicaNdiagonalVector(A,x,m,2)
guardarVector(x,m,"xExamen.txt")
eAbs=norma2(b-bEstimado,m)
eRel=eAbs/norma2(b,m)
println(@sprintf("El error absoluto de aproximación es de %g",eAbs))
println(@sprintf("El error relativo de aproximación es de %g",eRel))