function init_ns(job)
    # set up function handles
    m = job.dims[1]
    n = job.dims[2]
    println()
    println("-------------- FFT_plans -----------------------------------------")
    tic()
    tmp = rfft(zeros(Float64, 2*m+3,2*n+3))
    const FFT_big  = plan_rfft( zeros(Float64, 2*m+3,2*n+3), flags = FFTW.MEASURE)#plan_fft(zeros(2*m+3,2*n+3), flags=FFTW.MEASURE)
    const IFFT_big = plan_irfft(tmp, 2*m+3 , flags = FFTW.MEASURE)#plan_ifft(zeros(2*m+3,2*n+3), flags=FFTW.MEASURE)

    tmp = rfft(zeros(Float64, 2*m-1,2*n-1))
    const FFT      = plan_rfft( zeros(Float64, 2*m-1,2*n-1), flags = FFTW.MEASURE) #plan_fft(zeros(2*m-1,2*n-1), flags=FFTW.MEASURE)
    const IFFT     = plan_irfft(tmp, 2*m-1 , flags = FFTW.MEASURE)#plan_ifft(zeros(2*m-1,2*n-1), flags=FFTW.MEASURE)

    #FFT_dst2d  = plan_fft(zeros(2*(m+1),2*(n+1)), flags=FFTW.MEASURE)

    const R2R = plan_r2r(zeros(Float64, m,n), FFTW.RODFT00, flags=FFTW.MEASURE)
    const IR2R = plan_r2r(zeros(Float64, m,n), FFTW.RODFT00, flags=FFTW.MEASURE)
    toc()
    println("-------------- FFT_plans completed -------------------------------")
    println()
    println("-------------- building g_times ----------------------------------")
    tic()
    solver = Solver()
    solver.dt = job.dt

    const solver.g_times = w ::Array{Float64} ->
    apply_lgf(w, job.ghat, FFT, IFFT)#(IFFT * ( (FFT*([ w spzeros(size(w,1),size(w,2)-1); spzeros(size(w,1)-1,2 * size(w,2)-1)]) .*job.ghat )))[size(w,1):end,size(w,2):end] #apply_lgf(w, job.ghat, FFT, IFFT)
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
        const solver.ainv_times   = w :: Array{Float64,2} -> 1/(4 * (n+1) * (m+1)) * IR2R * (lam .* (R2R*w))  #4/( (n+1) * (m+1)) * dst2d(lam.*dst2d(w, FFT_dst2d),FFT_dst2d) #r2r(lam.*r2r(w, FFTW.RODFT00), FFTW.RODFT00) #4/( (n+1) * (m+1)) * dst( (dst( lam.*(dst( (dst(w))')') ))')'; #r2r(lam.*r2r(w, FFTW.RODFT01), FFTW.RODFT10) #DST * ( lam .* (DST * w) );
        const solver.la_inv_times = w :: Array{Float64,2} ->  solver.g_times(solver.ainv_times(w))
        job.z = schur_comp(job.dims, job.ec, solver.la_inv_times)
        const job.z = Symmetric(job.z)

        solver.zinv_times = w :: Array{Float64,2} -> -job.z\w
        job.cno = cond(job.z)
        solver.visc = w :: Array{Float64,2} -> (0.5*job.dt/job.r)*lap(w)
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
    wcrossx_x = - (yx' - job.center[2])
    wcrossx_y = (xy' - job.center[1])

    # nonlinear terms
    vel_a =
    (w :: Array{Float64,2},t :: Float64 ,u :: Array{Float64,2},v :: Array{Float64,2})->
    velocity(w, t, job.ghatbig,
                job.frame_rot, job.frame_xvel, job.frame_yvel,
                wcrossx_x, wcrossx_y,FFT_big,IFFT_big, u, v)

    solver.r1 =
    (w :: Array{Float64,2},t :: Float64,u :: Array{Float64,2},v :: Array{Float64,2}) ->
    -ctnonlin(w,t,vel_a,u,v)


    # Last bit
    solver.r2 = (t :: Float64) -> reshape([ -(job.frame_rot(t)* wcrossx_xb + job.frame_xvel(t)*wcrossx_xb);
                         -(job.frame_rot(t)* wcrossx_yb + job.frame_yvel(t)*wcrossx_yb) ], size(job.bdy,1)*2,1)

    solver.step = (soln :: Soln,slvr::Solver) -> cn_ab2(soln, slvr)
    solver.speed = (soln :: Soln,frame) -> vel_out(soln,3,job.ghatbig,
            job.frame_rot, job.frame_xvel, job.frame_yvel,
            wcrossx_x, wcrossx_y,FFT_big,IFFT_big)
    solver.cfl = (soln::Soln) -> solver.dt * maximum(maximum(vel_out(soln,3,job.ghatbig,
            job.frame_rot, job.frame_xvel, job.frame_yvel,
            wcrossx_x,wcrossx_y, FFT_big,IFFT_big) ))
    solver;
end

function apply_lgf(w :: Array{Float64,2}, lgfhat :: Array{Complex64,2}, FFT::FFTW_Type, IFFT::IFFTW_Type)
  M = size(w,1)
  N = size(w,2)
  padded_w =
  [ w zeros(M,N-1); zeros(M-1,N) zeros(M-1,N-1) ]
  # Tranform, multiply, inverse tranform [convolution]
  s_tmp :: Array{Float64,2} = (IFFT * ( (FFT*padded_w) .*lgfhat ))#real(ifft(fft(padded_w).*lgfhat)) #IFFT * ( (FFT*padded_w) .*lgfhat )
  return s_tmp[M:end,N:end];
end

function apply_lgf_big(w :: Array{Float64,2}, lgfhat :: Array{Complex64,2}, FFT::FFTW_Type, IFFT::IFFTW_Type)
  M = size(w,1)+2
  N = size(w,2)+2
  w2 = [ zeros(1,N); [ zeros(M-2,1) w zeros(M-2,1)]; zeros(1,N)]
  padded_w =
  [ w2 zeros(M,N-1); zeros(M-1,N) zeros(M-1,N-1) ]
  # Tranform, multiply, inverse tranform [convolution]
  s_tmp :: Array{Float64,2} = (IFFT * ( (FFT*padded_w) .*lgfhat ))#real(ifft(fft(padded_w).*lgfhat)) #IFFT * ( (FFT*padded_w) .*lgfhat )
  return s_tmp[M:end,N:end];
end


function schur_comp(mn, ec, opr :: Function)
  nb = size(ec,1)
  f = zeros(nb,nb)

  @simd for k=1:nb
      x = ( reshape( full(ec[k,:]),mn[1],mn[2]))
      @fastmath y = opr(x)
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