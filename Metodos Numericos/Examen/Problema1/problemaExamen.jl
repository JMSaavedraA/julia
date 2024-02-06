function crearMatrizExamenPentadiagonal(m::Int)
    A=zeros(Float64,m,5)
    for i=1:m
        A[i,1]=-4
        A[i,2]=-8
        A[i,3]=40
        A[i,4]=-8
        A[i,5]=-4
    end
    A[1,1]=0.0
    A[1,2]=0.0
    A[2,1]=0.0
    A[m-1,5]=0.0
    A[m,4]=0.0
    A[m,5]=0.0
    return A
end

function crearVectorExamen(m::Int)
    b=zeros(Float64,m).+100
    b[1]=20
    b[2]=50
    b[m-1]=50
    b[m]=20
    return b
end