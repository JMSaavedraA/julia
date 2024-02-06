using Plots;

f(x, y) = (x^2 + y^2 -1)^2 + (y^2 -1)^2

x = range(-1.5, 1.5, length=160)
y = range(-1.5, 1.5, length=160)
z = @. f(x', y)
plt1=contour(x, y, z,c=:acton,size=(800,600))
savefig(plt1,"ex1.png")