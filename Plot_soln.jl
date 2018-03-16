using Plots, Gadfly
plot(soln.f)

x = [30, 60, 90, 135, 200]
y = [0.846, 3.84, 13.46, 20.48, 51.48]
plot(x,y, Geom.point, Geom.smooth,
     Guide.xlabel("Stimulus"), Guide.ylabel("Response"), Guide.title("Dog Training"))
