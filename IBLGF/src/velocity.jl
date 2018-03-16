function velocity(w :: Array{Float64,2},
    depot :: Depot,#s_buffer :: Array{Float64,2}, s_buffer_2 ::Array{Complex{Float64},2}, s_buffer_3 ::Array{Complex{Float64},2},
    time :: Float64, lgfhat :: Array{Complex{Float64},2},
    frame_rot :: Function, frame_xvel :: Function, frame_yvel :: Function,
    wcrossx_x :: Array{Float64,2}, wcrossx_y :: Array{Float64,2},
    FFT_big :: FFTW_Type, IFFT_big :: IFFTW_Type, L_buffer :: Array{Float64,2})

       # velocity in various reference frames
       ω = frame_rot(time)
       u_x = frame_xvel(time)
       u_y = frame_yvel(time)
       depot.ufr .= ω .* wcrossx_x .+ u_x
       depot.vfr .= ω .* wcrossx_y .+ u_y

       apply_lgf_big(w, depot, lgfhat, FFT_big, IFFT_big, L_buffer)
       curl( L_buffer, depot.u_buffer, depot.v_buffer) # should be -L_buffer !! Carry the minus sign to next 2 lines
       depot.u .= .-depot.u_buffer .- depot.ufr
       depot.v .= .-depot.v_buffer .- depot.vfr
       nothing
end

function ctnonlin(w :: Array{Float64,2}, dp :: Depot, time::Float64, absvel::Function, L_buffer :: Array{Float64,2}, r1_out :: Array{Float64,2})
    #nx,ny = nonlin(w,time,absvel);
    #curltx(nx)+curlty(ny);
    #println("Nonlin")
    absvel(w, time, L_buffer);  # accelerating frame velocity computed with a buffer layer of cells added to w [make sure to set up vel_a with appropriate lgf]

    avg_wx(w, dp.wx_out) # shoulw be -w, we carry the minus with wx_out
    avg_wy(w, dp.wy_out)

    avg(dp.v, dp.avg_v)
    avg(dp.u, dp.avg_u)

    curltx(dp.wx_out.*dp.avg_v, dp.t1)
    curlty(dp.wy_out.*dp.avg_u, dp.t2)

    # Ke : -C X N(w) with a minus sign
    r1_out .= dp.t1 .- dp.t2
    nothing
end

function nonlin(w,time,vel_a)

    # absolute velocities
    #println("Nonlin")
    u, v = vel_a(w, time);  # accelerating frame velocity computed with a buffer layer of cells added to w [make sure to set up vel_a with appropriate lgf]

    # average vorticity
    wx = avg_wx(w);  # averages in y to x-cell edges
    wy = avg_wy(w);  # averages in x to y-cell edges
    # average velocities
    # vv = avg(v);     # averages v 4 nearest neighbors to x-cell edges
    # uu = avg(u);     # averages u 4 nearest neighbors to y-cell edges

    return -wx.*avg(v), wy.*avg(u)
end

function avg_wx(w :: Array{Float64,2}, wx_out::Array{Float64,2})
    # average vorticity [w] from vertices to x-faces [average in y dir()]
    m = size(w,1)
    n = size(w,2)

    # assumes vorticity is zero on a layer of cells on top and bottom of domain
    #wx = zeros(m,n+1)
    #wx[1:m,1] = 0.5*w[1:m,1]
    #wx[1:m,2:n] = 0.5*(w[1:m,1:n-1]+w[1:m,2:n])
    #wx[1:m,n+1] = 0.5*w[1:m,n]

    #wx
    #[ 0.25 * ( f[i,j] + f[i+1,j] + f[i,j+1] + f[i+1,j+1]) for i=1:size(f,1)-1, j=1:size(f,2)-1 ]

    @inbounds for i = 1:size(w,1), j=1:size(w,2)-1
        wx_out[i,j+1] = 0.5 .* (w[i,j] .+ w[i,j+1])
    end
    wx_out[:,1] .= 0.5 .* w[:,1]#[wx[i,1] = 0.5 * (w[i,1]) for i = 1:size(w,1) ]
    wx_out[:, n+1] .= 0.5 .* (w[:,n])#[wx[i, n+1] = 0.5 * (w[i,n]) for i = 1:size(w,1) ]

    nothing
end

function avg_wy(w :: Array{Float64,2}, wy_out::Array{Float64,2})
    # average vorticity [w] from vertices to y-faces [average in x dir()]
    m = size(w,1)
    n = size(w,2)
#
    #wy = zeros(m+1,n)
    # assumes vorticity is zero on a layer of cells on left and right of domain
    #wy[1,1:n] = 0.5*w[1,1:n]
    #wy[2:m,1:n] = 0.5*(w[1:m-1,1:n]+w[2:m,1:n])
    #wy[m+1,1:n] = 0.5*w[m,1:n]

    #wy
    @inbounds for i = 1:size(w,1)-1, j=1:size(w,2)
        wy_out[i+1, j] = 0.5 .* (w[i,j] .+ w[i+1,j])
    end
    wy_out[1,:]   .= 0.5 * w[1,:]#[wy[1  , j] = 0.5 * (w[1,j]) for j = 1:n ]
    wy_out[m+1,:] .= 0.5 * w[m,:]#[wy[m+1, j] = 0.5 * (w[m,n]) for j = 1:n ]

    nothing
end

function avg(f :: Array{Float64,2}, avg_out :: Array{Float64,2})
    # average 4 nearest neighbors
    #m = size(f,1)
    #n = size(f,2)
    #0.25*(f[1:m-1,1:n-1]+f[2:m,1:n-1]+f[1:m-1,2:n]+f[2:m,2:n])

    @inbounds for i=1:size(f,1)-1, j=1:size(f,2)-1
        avg_out[i,j] = 0.25 * ( f[i,j] + f[i+1,j] + f[i,j+1] + f[i+1,j+1])
    end
    nothing
end


function vel_out(soln:: Soln,comp,lgfhat,frame_rot, frame_xvel, frame_yvel,
     wcrossx_x,wcrossx_y,FFT_big,IFFT_big, w_buffer :: AbstractArray{Float64,2},
    s_buffer :: Array{Float64,2}, s_buffer_2 ::Array{Complex{Float64},2}, s_buffer_3 ::Array{Complex{Float64},2},
    L_buffer :: Array{Float64,2} )
# velocity [comp=1,2,3 for u, v, or speed]  (frame="l','a','r" for lab, accelerating, or relative vel.).
    #println("test")
    #println(soln.v)
    u,v = velocity(soln.w,w_buffer,
        s_buffer, s_buffer_2, s_buffer_3,
        soln.t,lgfhat,
        frame_rot, frame_xvel, frame_yvel,
        wcrossx_x,wcrossx_y,FFT_big,IFFT_big, L_buffer)

    if comp==1
        a = avg_x(u)
    elseif comp==2
        a = avg_y(v)
    else
        a = sqrt.( avg_x(u).^2 + avg_y(v).^2 )
    end
    return a
end

function avg_x(u :: Array{Float64,2}, avg_x_out :: Array{Float64,2})
    # average quantity from [extended] xfaces to vertices [averages in y dirn]
    m = size(u,1)
    n = size(u,2)-1

    #tmp = 0.5*(u[2:m-1,1:n]+u[2:m-1,2:n+1])
    @inbounds for i = 1:m-2, j=1:n
        avg_x_out[i,j] = 0.5 .* (u[i+1,j] + u[i+1, j+1] )
    end
    nothing
end

function avg_y(v :: Array{Float64,2}, avg_y_out :: Array{Float64,2})
    # average quantity from [extended] yfaces to vertic es [averages in x dirn]
    m = size(v,1)-1
    n = size(v,2)

    #0.5*(v[1:m,2:n-1]+v[2:m+1,2:n-1])
    @inbounds for i = 1:m, j=1:n-2
        avg_x_out[i,j] = 0.5 .* (u[i,j+1] + u[i+1, j+1] )
    end
end
