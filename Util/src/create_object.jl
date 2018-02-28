function create_circle(x, ds)
    # get points
    orig = x[:,1]
    rad = norm(x[:,1] - x[:,2])

    len = 2*pi*rad;
    np = floor(len/ds);
    err = abs(ds - len/np);

    ss = linspace(0,2*pi,np+1)';

    obj = Geo()
    obj.xyorig = x
    obj.xy = [orig[1,1] + rad*cos.(ss); orig[2,1] + rad*sin.(ss)];
    obj.xy = obj.xy[:,1:end-1]
    return obj

end
