using Printf

function Bisección(f::Function,a::AbstractFloat ,b::AbstractFloat,k::Int,tol::AbstractFloat)
    i=1;
    fA=f(a);
    fB=f(b);
    hasFailed= (fA*fB)>tol;
    
    foundRoot= (abs(fA)<tol || abs(fB)<tol);
    
    while i<k && !foundRoot && !hasFailed
        c=(a+b)/2.0;
        fC=f(c);
        foundRoot=abs(fC)<tol;
        if (fA*fC)>0
            a=c;
            fA=fC;
        else
            b=c;
            fB=fC;
        end
        i+=1;
    end
    if i>1
        if foundRoot
            @printf("La raíz encontrada es x=%0.20f",c)
            println()
            @printf("Se realizaron %i iteraciones", i)
            println()
            return c
        elseif hasFailed
            print("Error")
            @printf("a=%0.20f, f(a)=%0.20f", a,fA)
            println()
            @printf("b=%0.20f, f(b)=%0.20f", b,fB)
            println()
            @printf("c=%0.20f, f(c)=%0.20f", c,fC)
            println()
            @printf("Se realizaron %i iteraciones", i)
            println()
        else
            println("Raíz no encontrada")
            @printf("a=%0.20f, f(a)=%0.20f", a,fA)
            println()
            @printf("b=%0.20f, f(b)=%0.20f", b,fB)
            println()
            @printf("c=%0.20f, f(c)=%0.20f", c,fC)
            println()
            @printf("Se realizaron %i iteraciones", i)
            println()
        end
    else
        if foundRoot
            if abs(fA)<tol
                c=a;
                fC=fA;
                @printf("El valor introducido a=%0.20f es una raíz de f", a)
                println()
            else
                c=b;
                fC=fB;
                @printf("El valor introducido b=%0.20f es una raíz de f", b)
                println()
            end
        else
            println("Los valores de f(a) y f(b) tienen el mismo signo")
        end

    end
end
    
function newtonRaphson(f::Function,fPrime::Function,x::AbstractFloat,k::Int,tol::AbstractFloat)
    i=1;
    bX=x;
    fX=f(x);
    afX=abs(fX);
    bFX=afX;
    foundRoot=afX<tol;
    fPX=fPrime(x);
    hasFailed=abs(fPX)<tol;
    if foundRoot
        @printf("El valor introducido x=%0.20f es una raíz de f", x)
    end
    while i<k && !foundRoot && !hasFailed
        x-=(fX/fPX);
        fX=f(x);
        afX=abs(fX);
        foundRoot=afX<tol;
        if afX<bFX
            bFX=afX;
            bX=x;
        end
        fPX=fPrime(x);
        hasFailed=abs(fPX)<tol;
        i+=1;
    end
    if foundRoot
        @printf("La raíz encontrada es x=%0.20f", x)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
        return x
    elseif hasFailed
        println("Error")
        @printf("x=%0.20f, f(x)=%0.20f, f'(x)=%0.20f", x,fX,fPX)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
    else
        println("Raíz no encontrada")
        @printf("Mejor x=%0.20f, f(x)=%0.20f", bX,bFX)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
    end
end

function newtonSinDerivada(f::Function,x::AbstractFloat,k::Int,tol::AbstractFloat,h::AbstractFloat)
    i=1;
    bX=x;
    fX=f(x);
    fXP=f(x+h);
    fXM=f(x-h);
    afX=abs(fX);
    bFX=afX;
    foundRoot=afX<tol;
    fPX=(fXP - fXM)/(2.0*h);
    hasFailed=abs(fPX)<tol;
    if foundRoot
        @printf("El valor introducido x=%0.20f es una raíz de f", x)
        println()
    end
    while i<k && !foundRoot && !hasFailed
        x-=(fX/fPX);
        fX=f(x);
        fXP=f(x+h);
        fXM=f(x-h);
        afX=abs(fX);
        foundRoot=afX<tol;
        if afX<bFX
            bFX=afX;
            bX=x;
        end
        fPX=(fXP - fXM)/(2.0*h);
        hasFailed=abs(fPX)<tol;
        i+=1;
    end
    if foundRoot
        @printf("La raíz encontrada es x=%0.20f", x)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
        return x
    elseif hasFailed
        println("Error")
        @printf("x=%0.20f, f(x)=%0.20f, f'(x)=%0.20f", x,fX,fPX)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
    else
        println("Raíz no encontrada")
        @printf("Mejor x=%0.20f, f(x)=%0.20f", bX,bFX)
        println()
        @printf("Se realizaron %i iteraciones", i)
        println()
    end
end

