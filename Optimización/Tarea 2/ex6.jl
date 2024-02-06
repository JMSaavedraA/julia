using Plots;

f(x, y) = 100(y-x^2)^2 + (1-x)^2

x = range(0.5, 1.5, length=160)
y = range(0.5, 1.5, length=120)
z = @. f(x', y)
plt1=contour(x, y, z,c=:berlin,size=(800,600))
savefig(plt1,"ex6.png")