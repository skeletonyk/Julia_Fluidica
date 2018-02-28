function time_step(parms, job, solver)
#
    tic()
    job.soln.w = zeros(job.dims[1],job.dims[2])
    job.soln.u = zeros(job.dims[1]+2,job.dims[2]+1)
    job.soln.v = zeros(job.dims[1]+1,job.dims[2]+2)
    job.soln.r1 = zeros(job.dims[1],job.dims[2])
    job.soln.t = 0

    first::Bool = true
    for k=1:job.nstp
        w = deepcopy(job.soln.w)
        t = job.soln.t
        r1old = deepcopy(job.soln.r1)
        @time @fastmath job.soln.w, job.soln.t, job.soln.f, job.soln.r1  = solver.step(w, t, solver, first, r1old)
        job.it = job.it + 1
        cfl:: Float64 = solver.cfl(job.soln)

        if cfl > 1.5
            error("The CFL has become too large and Fluidica cannot continue().  Consider choosing a lower time step so that CFL < 1.5")
            return
        end
        first = false
    end
    toc()
    #job.runtime = toc()
    job.checkpoint = job.checkpoint+1
    #job.data[job.checkpoint] = job.soln

    return job
end
