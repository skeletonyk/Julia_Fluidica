function velocity(w :: Array{Float64,2}, time :: Float64, lgfhat :: Array{Complex64,2},
    frame_rot :: Function, frame_xvel :: Function, frame_yvel :: Function,
    wcrossx_x :: Array{Float64,2}, wcrossx_y :: Array{Float64,2},
    FFT_big :: FFTW_Type, IFFT_big :: IFFTW_Type)

       # velocity in various reference frames
       ufr = frame_rot(time).* wcrossx_x + frame_xvel(time).*ones(size(wcrossx_x))
       vfr = frame_rot(time).* wcrossx_y + frame_yvel(time).*ones(size(wcrossx_y))

       M = size(w,1)+2;
       N = size(w,2)+2;
       #w2 = [ zeros(1,N); [ zeros(M-2,1) w zeros(M-2,1)]; zeros(1,N)]
       #println("apply_lgf_biiig")
       s = -apply_lgf_big(w,lgfhat, FFT_big, IFFT_big)
       u, v = curl(s)
       u = u - ufr
       v = v - vfr

       return u,v
end

@inline function ctnonlin(w :: Array{Float64,2},time::Float64,absvel::Function)
    #nx,ny = nonlin(w,time,absvel);
    #curltx(nx)+curlty(ny);
    u,v = absvel(w, time);  # accelerating frame velocity computed with a buffer layer of cells added to w [make sure to set up vel_a with appropriate lgf]

    #wx = avg_wx(w);  # averages in y to x-cell edges
    #wy = avg_wy(w);  # averages in x to y-cell edges
    # average velocities
    # vv = avg(v);     # averages v 4 nearest neighbors to x-cell edges
    # uu = avg(u);     # averages u 4 nearest neighbors to y-cell edges

    curltx(-avg_wx(w).*avg(v))+ curlty(avg_wy(w).*avg(u))

end

function nonlin(w,time,vel_a)

    # absolute velocities
    u, v = vel_a(w, time);  # accelerating frame velocity computed with a buffer layer of cells added to w [make sure to set up vel_a with appropriate lgf]

    # average vorticity
    wx = avg_wx(w);  # averages in y to x-cell edges
    wy = avg_wy(w);  # averages in x to y-cell edges
    # average velocities
    # vv = avg(v);     # averages v 4 nearest neighbors to x-cell edges
    # uu = avg(u);     # averages u 4 nearest neighbors to y-cell edges

    return -wx.*avg(v), wy.*avg(u)
end

function avg_wx(w :: Array{Float64,2})
    # average vorticity [w] from vertices to x-faces [average in y dir()]
    m = size(w,1)
    n = size(w,2)

    # assumes vorticity is zero on a layer of cells on top and bottom of domain
    wx = zeros(m,n+1)
    wx[1:m,1] = 0.5*w[1:m,1]
    wx[1:m,2:n] = 0.5*(w[1:m,1:n-1]+w[1:m,2:n])
    wx[1:m,n+1] = 0.5*w[1:m,n]

    wx
end

function avg_wy(w :: Array{Float64,2})
    # average vorticity [w] from vertices to y-faces [average in x dir()]
    m = size(w,1)
    n = size(w,2)

    wy = zeros(m+1,n)
    # assumes vorticity is zero on a layer of cells on left and right of domain
    wy[1,1:n] = 0.5*w[1,1:n]
    wy[2:m,1:n] = 0.5*(w[1:m-1,1:n]+w[2:m,1:n])
    wy[m+1,1:n] = 0.5*w[m,1:n]

    wy
end

function avg(f :: Array{Float64,2})
    # average 4 nearest neighbors
    m = size(f,1)
    n = size(f,2)
    0.25*(f[1:m-1,1:n-1]+f[2:m,1:n-1]+f[1:m-1,2:n]+f[2:m,2:n])
end


function vel_out(soln:: Soln,comp,lgfhat,frame_rot, frame_xvel, frame_yvel, wcrossx_x,wcrossx_y,FFT_big,IFFT_big)
# velocity [comp=1,2,3 for u, v, or speed]  (frame="l','a','r" for lab, accelerating, or relative vel.).
    #println("test")
    #println(soln.v)
    u,v = velocity(soln.w,soln.t,lgfhat,
        frame_rot, frame_xvel, frame_yvel,
        wcrossx_x,wcrossx_y,FFT_big,IFFT_big)
    if comp==1
        a = avg_x(u)
    elseif comp==2
        a = avg_y(v)
    else
        a = sqrt.( avg_x(u).^2 + avg_y(v).^2 )
    end

    return a
end

function avg_x(u :: Array{Float64,2})
    # average quantity from [extended] xfaces to vertices [averages in y dirn]
    m = size(u,1)
    n = size(u,2)-1

    0.5*(u[2:m-1,1:n]+u[2:m-1,2:n+1])
end

function avg_y(v :: Array{Float64,2})
    # average quantity from [extended] yfaces to vertices [averages in x dirn]
    m = size(v,1)-1
    n = size(v,2)

    0.5*(v[1:m,2:n-1]+v[2:m+1,2:n-1])
end
