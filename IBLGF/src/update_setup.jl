function update_setup(parms, obj)
    tic()
    job = Job()

    # discrete delta function()
    ddf = getfield(Ddfs, Symbol(parms.ddf))
    job.ddf = ddf
    pts = obj.xy
    job.pts = pts;

    # Mesh
    m = floor( (parms.axlim[2]-parms.axlim[1])*parms.ppl +1)
    n = floor( (parms.axlim[4]-parms.axlim[3])*parms.ppl +1)
    job.dims = [m, n]
    job.domain = [parms.axlim[1] (parms.axlim[1]+(m-1)/parms.ppl) parms.axlim[3] (parms.axlim[3]+(n-1)/parms.ppl)]
    job.center = [ (parms.orig[1]-parms.axlim[1])*parms.ppl+1  [parms.orig[2]-parms.axlim[3]]*parms.ppl+1  ]
    xg,yg = meshgrid(1:m,1:n);  # the grid()
    job.x = ((xg-1.0+job.domain[1]*parms.ppl)/parms.ppl)';
    job.y = ((yg-1.0+job.domain[3]*parms.ppl)/parms.ppl)';
    job.r = parms.rey/parms.ppl

    xb = (pts[1,:]-parms.axlim[1])*parms.ppl+ones(size(pts[1,:]))
    yb = (pts[2,:]-parms.axlim[3])*parms.ppl+ones(size(pts[1,:]))
    job.bdy = [xb yb]

    # Inital soln
    job.soln = Soln()
    job.soln.t = parms.hzn[1];

    # Time stepping
    job.dt, job.nstp = get_dt(parms)
    job.it = 0
    job.checkpoint = 1
    job.times = parms.hzn[1]:parms.hzn[2]:parms.hzn[3]

    job.frame_rot  = t:: Float64-> parms.mot.r(t/parms.ppl)/parms.ppl;
    job.frame_xvel = t:: Float64-> parms.mot.x(t/parms.ppl);
    job.frame_yvel = t:: Float64-> parms.mot.y(t/parms.ppl);

    const job.g = setup_lap(job.dims);
    gbig = setup_lap(job.dims + [2, 2]);
    g_p = setup_lap(job.dims + [1, 1]);

    const job.ghat = transform_lgf(job.g);
    const job.ghatbig = transform_lgf(gbig);
    const job.ghat_p= transform_lgf(g_p);

    # Compute the regularlization matrices
    const job.etx, job.ety = et_mat(job.dims,job.bdy,job.ddf);
    const job.ec = ec_mat(job.dims,job.etx,job.ety);

    toc()
    return job
end
