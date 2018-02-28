mutable struct geo
    xyorig
    xy
    Geo() = new()
end

function create_circle(x)
    # get points

    orig = x[:,1]
    print(x)
    rad = norm(x[:,1] - x[:,2])
    println(rad)
    ss = linspace(0,2*π,101)

    obj = Geo()
    obj.xyorig = x
    obj.xy = [orig[1,1] + rad*cos(ss); orig[2,1] + rad*sin(ss)];
end
