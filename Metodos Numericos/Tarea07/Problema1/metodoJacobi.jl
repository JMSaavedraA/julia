include("Solvers.jl")

function eigenJacobi(A::AbstractMatrix, m::Int, tol::AbstractFloat,maxIter::Int)
    #Metodo de Jacobi para obtener los eigenpares de una matriz simétrica por medio de rotaciones
    Aaux=copy(A) #Vamos a rotar sobre la matriz original
    B=Identidad(m) #Iniciamos la matriz de los eigenvectores de A como la identidad
    v,i,j=maxAbsolutoFueraDiagonal(Aaux,m)
    k=0
    while v>tol && k<maxIter
        t=atan(2*v,Aaux[i,i]-Aaux[j,j])/2 #Obtenemos el ángulo por medio de atan2
        ct=cos(t) #Solo evaluamos una vez el seno y coseno del ángulo
        st=sin(t)
        for l=1:m
            ali=Aaux[l,i]
            alj=Aaux[l,j]
            Aaux[l,i]=ali*ct+alj*st #Actualizamos AR'
            Aaux[l,j]=alj*ct-ali*st
            bli=B[l,i]
            blj=B[l,j]
            B[l,i]=bli*ct+blj*st #Actualizamos BR'
            B[l,j]=blj*ct-bli*st
        end
        for l=1:m
            ali=Aaux[i,l]
            alj=Aaux[j,l]
            Aaux[i,l]=ali*ct+alj*st #Actualizamos RAR'
            Aaux[j,l]=alj*ct-ali*st
        end
        v,i,j=maxAbsolutoFueraDiagonal(Aaux,m) #Buscamos el valor más grande fuera de la diagonal
        k=k+1
    end
    a=vectorDiagonal(Aaux,m)
    p=sortperm(abs.(a)) #Ordenamos de menor a mayor
    Aaux.=Aaux[:,p]
    B.=B[:,p]
    a.=a[p]
    return B,a,k,Aaux
end
