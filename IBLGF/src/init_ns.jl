function init_ns(job:: Job_static)
    # set up function handles
    m = job.dims[1]
    n = job.dims[2]
    println()
    println("-------------- FFT_plans -----------------------------------------")
    tic()
    tmp = rfft(rand(Float64, 2*m+3,2*n+3))
    const FFT_big  = plan_rfft( rand(Float64, 2*m+3,2*n+3), flags = FFTW.MEASURE)#plan_fft(zeros(2*m+3,2*n+3), flags=FFTW.MEASURE)
    const IFFT_big = plan_irfft(tmp, 2*m+3 , flags = FFTW.MEASURE)#plan_ifft(zeros(2*m+3,2*n+3), flags=FFTW.MEASURE)

    tmp = rfft(rand(Float64, 2*m-1,2*n-1))
    const FFT      = plan_rfft( rand(Float64, 2*m-1,2*n-1), flags = FFTW.MEASURE) #plan_fft(zeros(2*m-1,2*n-1), flags=FFTW.MEASURE)
    const IFFT     = plan_irfft(tmp, 2*m-1 , flags = FFTW.MEASURE)#plan_ifft(zeros(2*m-1,2*n-1), flags=FFTW.MEASURE)

    #FFT_dst2d  = plan_fft(zeros(2*(m+1),2*(n+1)), flags=FFTW.MEASURE)

    const R2R = plan_r2r(rand(Float64, m,n), FFTW.RODFT00, flags=FFTW.MEASURE)
    const IR2R = plan_r2r(rand(Float64, m,n), FFTW.RODFT00, flags=FFTW.MEASURE)

    depot = depot_ini(m,n)
    #w_buffer = zeros(Float64, 2*m+3, 2*n+3)
    #s_buffer = zeros(Float64, 2*m+3, 2*n+3)
    #s_buffer_2 = rfft(s_buffer)
    #s_buffer_3 = similar(s_buffer_2)
    #w_buffer_small = zeros(Float64, 2*m-1, 2*n-1)
    #s_buffer_small = zeros(Float64, 2*m-1, 2*n-1)
    #s_buffer_small_2 = rfft(s_buffer_small)
    #s_buffer_small_3 = similar(s_buffer_small_2)

    toc()
    println("-------------- FFT_plans completed -------------------------------")
    println()
    println("-------------- building g_times ----------------------------------")
    tic()
    solver = Solver()
    solver.dt = job.dt

    const solver.g_times = (w ::Array{Float64,2}, L_buffer ::Array{Float64,2})->
    apply_lgf(w, depot, job.ghat, FFT, IFFT, L_buffer)

    #(IFFT * ( (FFT*([ w spzeros(size(w,1),size(w,2)-1); spzeros(size(w,1)-1,2 * size(w,2)-1)]) .*job.ghat )))[size(w,1):end,size(w,2):end] #apply_lgf(w, job.ghat, FFT, IFFT)
    #job.t = schur_comp(job.dims, job.ec, solver.g_times)
    #const job.t = Symmetric(job.t)
    #const solver.tinv_times = z :: AbstractVector ->  job.t\z

    toc()
    println("-------------- building g_times completed ------------------------")
    println()
    println("-------------- extra operators for ab2/cn scheme------------------")
    tic()
    # extra operators for ab2/cn scheme
        lam = zeros(m,n)
        for i=1:m
            for j=1:n
                lam[i,j] = 1./(1 - (.5*job.dt/job.r)*( 2*cos(pi*i/(m+1)) + 2*cos(pi*j/(n+1)) - 4))
            end
        end

        #DST =1/sqrt(2*(n+1) * 2 * (m+1)) * plan_r2r(lam, FFTW.RODFT00)
        solver.ainv_times   = w :: Array{Float64,2} -> 1/(4 * (n+1) * (m+1)) * IR2R * (lam .* (R2R*w))  #4/( (n+1) * (m+1)) * dst2d(lam.*dst2d(w, FFT_dst2d),FFT_dst2d) #r2r(lam.*r2r(w, FFTW.RODFT00), FFTW.RODFT00) #4/( (n+1) * (m+1)) * dst( (dst( lam.*(dst( (dst(w))')') ))')'; #r2r(lam.*r2r(w, FFTW.RODFT01), FFTW.RODFT10) #DST * ( lam .* (DST * w) );

        L_buffer = zeros(m,n)
        solver.la_inv_times = (w :: Array{Float64,2}, L_buffer :: Array{Float64,2}) -> solver.g_times(solver.ainv_times(w), L_buffer)
        z = schur_comp(job.dims, job.ec, solver.la_inv_times)
        z = factorize(Symmetric(z))
        solver.zinv_times = w :: Array{Float64,2} -> z\(-w)
        #job.cno = cond(job.z)
        solver.visc = (w :: Array{Float64,2}, w_out :: Array{Float64,2}) -> lap(w, w_out,0.5*job.dt/job.r)
    toc()
    println("-------------- extra operators for ab2/cn scheme completed--------")
    # extra operators FOR ab2/cn scheme / end

    const solver.b_times = z :: Array{Float64,2}  -> job.ec * reshape(z,m*n,1);
    const solver.bt_times = z :: Array{Float64,2} -> reshape(job.ec' * z, m,n);

    wcrossx_xb = - (job.bdy[:,2]-job.center[2]);
    wcrossx_yb = (job.bdy[:,1]-job.center[1]);

    # ----------------------------------------
    ~,yx = meshgrid(0.0:m+1,0.5:n+0.5)
    xy,~ = meshgrid(0.5:m+0.5,0.0:n+1)
    wcrossx_x :: Array{Float64,2} = - (yx' - job.center[2])
    wcrossx_y :: Array{Float64,2} = (xy' - job.center[1])

    # nonlinear terms

    vel_a =
    (w :: Array{Float64,2}, t :: Float64, L_buffer)->
    velocity(w, depot,#w_buffer, s_buffer, s_buffer_2, s_buffer_3
             t, job.ghatbig,
             job.frame_rot, job.frame_xvel, job.frame_yvel,
             wcrossx_x, wcrossx_y,FFT_big,IFFT_big, L_buffer)

    solver.vel_a = vel_a
    solver.r1 = (w :: Array{Float64,2}, depot :: Depot, t :: Float64, r1_vel_buffer :: Array{Float64,2}, r1 :: Array{Float64,2}) ->
    ctnonlin(w, depot, t, vel_a, r1_vel_buffer, r1)

    # Last bit
    solver.r2 = (t :: Float64) -> reshape([ -(job.frame_rot(t)* wcrossx_xb + job.frame_xvel(t)*wcrossx_xb);
                         -(job.frame_rot(t)* wcrossx_yb + job.frame_yvel(t)*wcrossx_yb) ], size(job.bdy,1)*2,1)

    solver.step = (soln:: Soln, h :: Solver_static, first::Bool,
     r1 :: Array{Float64,2}, rhs :: Array{Float64,2}, wstar :: Array{Float64,2},
     tmp_Linv_wstar:: Array{Float64,2}, r1_vel_buffer:: Array{Float64,2}) ->
     cn_ab2(soln, h, first, r1, rhs, wstar, tmp_Linv_wstar, r1_vel_buffer, depot)

    solver.speed = (soln :: Soln,frame) -> vel_out(soln,3,job.ghatbig,
            job.frame_rot, job.frame_xvel, job.frame_yvel,
            wcrossx_x, wcrossx_y,FFT_big,IFFT_big,w_buffer,s_buffer, s_buffer_2, s_buffer_3, L_buffer)
    solver.cfl = soln::Soln -> solver.dt * maximum(maximum(vel_out(soln,3,job.ghatbig,
            job.frame_rot, job.frame_xvel, job.frame_yvel,
            wcrossx_x,wcrossx_y, FFT_big,IFFT_big,w_buffer,s_buffer, s_buffer_2, s_buffer_3, L_buffer) ))

    solver_static = Solver_static(solver.dt,
            solver.g_times,
            solver.ainv_times,
            solver.la_inv_times,
            solver.zinv_times,
            solver.visc,
            solver.b_times,
            solver.bt_times,
            solver.r1,
            solver.r2,
            solver.step,
            solver.speed,
            solver.cfl,
            solver.vel_a)

    return solver_static
end

function apply_lgf(w :: Array{Float64,2}, dp::Depot,
    lgfhat :: Array{Complex{Float64},2}, FFT::FFTW_Type, IFFT::IFFTW_Type, L_buffer :: Array{Float64,2} )
  M = size(w,1)
  N = size(w,2)
  #padded_w :: Array{Float64,2} =
  #[ w zeros(M,N-1); zeros(M-1,N) zeros(M-1,N-1) ]
  dp.w_buffer_small[1:M, 1:N] = w
  A_mul_B!(dp.s_buffer_small_2, FFT, dp.w_buffer_small)
  dp.s_buffer_small_3 .= dp.s_buffer_small_2.*lgfhat
  A_mul_B!(dp.s_buffer_small, IFFT, dp.s_buffer_small_3)
  L_buffer .= dp.s_buffer_small[M:end,N:end];

  nothing
end

function apply_lgf_big(w :: Array{Float64,2},
    dp :: Depot,
    lgfhat :: Array{Complex{Float64},2}, FFT::FFTW_Type, IFFT::IFFTW_Type, L_buffer :: Array{Float64,2})

  M = size(w,1)+2
  N = size(w,2)+2

  dp.w_buffer[2:M-1, 2:N-1] = w
  A_mul_B!(dp.s_buffer_2, FFT, dp.w_buffer)
  dp.s_buffer_3 .= dp.s_buffer_2.*lgfhat
  A_mul_B!(dp.s_buffer, IFFT, dp.s_buffer_3)
  L_buffer .= dp.s_buffer[M:end,N:end];

  nothing
end


function schur_comp(mn, ec::AbstractArray, opr :: Function)
  nb = size(ec,1)
  f = zeros(nb,nb)
  y = zeros(Float64, mn[1], mn[2])
  @simd for k=1:nb
      x = reshape( full(ec[k,:]),mn[1],mn[2])
      @fastmath opr(x,y)
      f[:,k] = -ec*(reshape(y,mn[1]*mn[2],1))
  end
  f
end


# Test DST
function dst(a :: Array{Float64,2})
    n = size(a,1);
    m = size(a,2)
    # Pad or truncate a if necessary

    y=zeros(2*(n+1),m)
    y[2:n+1,:]=a
    y[n+3:2*(n+1),:]=-flipdim(a, 1)
    yy=fft(y,1);
    b=yy[2:n+1,:]/(-2*im)

    real(b);
end

function dst2d(a, FFT)
    n = size(a,1);
    m = size(a,2)
    # Pad or truncate a if necessary

    y=zeros(2*(n+1), 2*(m+1))
    y[2:n+1,2:m+1]=a
    y[n+3:2*(n+1), m+3:2*(m+1)]=flipdim(flipdim(a,1),2)
    y[n+3:2*(n+1), 2:m+1]=-flipdim(a, 1)
    y[2: n + 1, m + 3:2 * (m + 1)]=-flipdim(a, 2)
    yy=FFT*y;
    b=yy[2:n+1,2:m+1]/(-4)

    real(b);
end
# LocalWords:  Nordmark
