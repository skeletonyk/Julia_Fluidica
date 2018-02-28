using Plots
using Util
using PyPlot
pyplot()
m,n = 181,181
a = reshape(full(job.ec[1,:]),181,181)
#x,y =meshgrid( 1.0:m, 0.5:n+0.5 )
x = 1:181.0
y = 0.5:181.5

p1 =Plots.contour(x,x,a,fill=true)
Plots.plot(p1)
#Plots.plot(x,y,a)
