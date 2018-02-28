
function et_mat(mn,bdy,useddf)

    m = mn[1]
    n = mn[2]

    nb = size(bdy,1)

    # first determine the compactness of the discrete delta function()
    rr = 0:0.5:10
    r = rr[findfirst(x->abs.(x)<1e-14, useddf(rr))]
    r = r+ 0.5
    # x-part: interpolate to x-cell()-edges

    xx,yx =meshgrid( 1.0:m, 0.5:n+0.5 )
    xx = reshape(xx',m*(n+1),1)
    yx = reshape(yx',m*(n+1),1)

    ex = spzeros(m*(n+1),nb)

    for k=1:nb
        newpts = find((abs.(bdy[k,1]-xx) .<= r) .& (abs.(bdy[k,2]-yx) .<= r))
        ex[newpts,k] = useddf(bdy[k,1]-xx[newpts]).*useddf( bdy[k,2]-yx[newpts])
    end

    # y-part: interpolate to y-cell()-edges

    xy,yy =meshgrid( 0.5:m+0.5,1.0:n )
    xy = reshape(xy',(m+1)*n,1)
    yy = reshape(yy',(m+1)*n,1)

    ey = spzeros((m+1)*n,nb)

    for k=1:nb
        newpts = find((abs.(bdy[k,1]-xy) .<= r) .& (abs.(bdy[k,2]-yy) .<= r))
        ey[newpts,k] = useddf(bdy[k,1]-xy[newpts]).*useddf(bdy[k,2]-yy[newpts])
    end

    return ex, ey

end
