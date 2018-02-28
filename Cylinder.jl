module Cylinder
    # Add path
    export job_run
    push!(LOAD_PATH, pwd())
    FFTW.set_num_threads(Sys.CPU_CORES * 2)


    # Modules
    using Util
    using Ddfs
    using Lgf
    using IBLGF
    using FFTW


    # DataTypes
    mutable struct Mot
        x :: Function
        y :: Function
        r :: Function
        Mot() = new()
    end

    mutable struct Tstep
        auto; cfl
        Tstep() = new()
    end

    mutable struct Parms
        solver;buffer;axlim;ppl;orig;rey;tstep;tstart;tend;tinc;ibph;ddf;hzn
        base :: Function
        ic :: Function
        mot :: Mot
        Parms() = new()
    end

    function job_run(ppl)
        # ------------------------------------------------------------------------------
        box = 3
        rey = 10
        cfl = .2
        ibph = 2
        ddf = "yang3"
        # ------------------------------------------------------------------------------
        parms = Parms()

        parms.solver = 0         # which solver to use
        parms.buffer = 1         # 1 = store checkpoints in memory(); 0 = write to separate files
        parms.axlim = [-1 1 -1 1] * box         # domain size (body coords)
        parms.ppl = ppl         # resolution [grid points per unit length = 1/h];
        parms.orig = [0  0]         # location of origin in body coords
        parms.tstart = 0         # starting time
        parms.tend = 1;         # ending time
        parms.tinc = 0.1         # time increment
        parms.mot = Mot()
        parms.mot.x = t -> 0 # zeros(size(t))         # motion [velocity] of body in x dir()
        parms.mot.y = t -> 0 # zeros(size(t))         # motion [velocity] of body in y dir()
        parms.mot.r = t -> 1 # ones(size(t))         # motion [rotation rate] of body in theta dir()
        parms.rey = rey;         # Reynolds number based on L
        parms.tstep = Tstep()
        parms.tstep.auto = 1;         # automatically choose dt
        parms.tstep.cfl = cfl         #
        parms.ibph = ibph         # IB point spacing [rel. to grid()]
        parms.base = (x,y) -> zeros(size(x))         # base solution/initial condition
        parms.ddf = ddf         # discrete delta function()
        parms.ic = (x,y) -> zeros(size(x))         # base solution/initial condition
        parms.hzn = [0, 0.1, 0.1]

        # ------------------------------------------------------------------------------
        # Circle for now
        xobj = [ 0.0 0.0 0.0 ; 0.0 1.0 0.0]
        obj = Geo()
        obj = Util.create_circle(xobj, parms.ibph./parms.ppl)

        # ------------------------------------------------------------------------------
        job = update_setup(parms,obj);
        solver_static = init_ns(job);
        job = time_step(parms,job,solver_static)
        #x,y,u = post_proc[parms,job,solver,job.soln,"v [lab]",0]

        # ------------------------------------------------------------------------------
        parms, job, solver_static;
    end

end
