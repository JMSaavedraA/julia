using Printf

function optimizadorNewton(f::Function,x::AbstractFloat,k::Int,tol::AbstractFloat,h::AbstractFloat)
    i=1
    fX=f(x)/h
    fXP=f(x+h)/h
    fXM=f(x-h)/h
    op1=(fXP - fXM)
    op2=(fXP + fXM - 2.0*fX)
    fPX=op1/(2.0*h)
    fDX=op2/(h^2)
    afPX=abs(fPX)
    foundCritical=afPX<tol
    hasFailed=abs(fDX)<tol
    if foundCritical
        @printf("El valor introducido x=%0.20f es un punto crítico de f(x)", x)
        println()
    end
    while i<k && !foundCritical && !hasFailed
        x-=(h*op1)/(2*op2)
        fX=f(x)
        fXP=f(x+h)
        fXM=f(x-h)
        op1=(fXP - fXM)
        op2=(fXP + fXM - 2.0*fX)
        fPX=op1/(2.0*h)
        fDX=op2/(h^2)
        afPX=abs(fPX)
        foundCritical=afPX<tol
        hasFailed=abs(fDX)<tol
        i+=1
    end
    if foundCritical
        if fDX>tol
            @printf("El mínimo encontrado es x=%0.20f", x)
            println()
        elseif fDX<-tol
            @printf("El máximo encontrado es x=%0.20f", x)
            println()
        else
            @printf("Punto crítico encontrado, x=%0.10f, f'(x)=%0.10f, f''(x)=%0.10f", x,fPX,fDX)
            println()
        end
        return x
    else
        println("Error")
        @printf("x=%0.10f, , f'(x)=%0.10f, f''(x)=%0.10f, i=%i", x,fPX,fDX,i)
        println()
    end
end