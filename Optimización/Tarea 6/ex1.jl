using Plots

x = LinRange(-1.2,1.2,200)
y = LinRange(-1.2,1.2,200)


f(x, y) = (x^2 + y^2 -1)^2 + (y^2 -1)^2

plt1=plot(x,y,z,xlims=[-1.2,1.2],ylims=[-1.2,1.2],zlims=[0,2],st=:surface,color=:plasma,size=[800,600])
savefig(plt1,"ex1.png")
