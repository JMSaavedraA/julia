using Plots, LaTeXStrings

n=1001 #Tamaño del vector de x
x=LinRange(-11,11,n); y1=sqrt.(12 .+ x.^2); y2=8 ./ x; #Vector de x y de la curva de nivel f_2(x,y)=16
X=[x;x[n:-1:1]]; Y1=[y1;-y1[n:-1:1]]; #Vector auxiliar de x y de la curva de nivel f_1(x,y)=12
x0=[-2;2];y0=2*x0 #Puntos de intersección

plt1=plot(X,Y1, xlims=(-10,10),ylims=(-7,7),size=(800,600), aspect_ratio = 1,color_palette=:seaborn_bright,label=L"f_1(x_1,x_2)=12"); #Gráfica de la curva de nivel f_1(x,y)=12
plot!(x,y2,label=L"f_2(x_1,x_2)=16");#Gráfica de la curva de nivel f_2(x,y)=16
scatter!(x0,y0,label=L"f(x_1,x_2)=[12,16]^T")#Puntos donde f(x,y)=[12,16]'
savefig(plt1,"ej1.png")