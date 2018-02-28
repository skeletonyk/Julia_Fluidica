function time_step(parms,job,solver)
#
    tic()
    job.soln.w = zeros(job.dims[1],job.dims[2])
    job.soln.u = zeros(job.dims[1]+2,job.dims[2]+1)
    job.soln.v = zeros(job.dims[1]+1,job.dims[2]+2)
    job.soln.t = 0
    for k=1:job.nstp
        @time @fastmath job.soln = solver.step(job.soln,solver)
        job.it = job.it + 1
        cfl = solver.cfl(job.soln)

        if isnan(cfl) || cfl > 1.5
            error("The CFL has become too large and Fluidica cannot continue().  Consider choosing a lower time step so that CFL < 1.5")
            return
        end
    end
    toc()
    #job.runtime = toc()
    job.checkpoint = job.checkpoint+1
    #job.data[job.checkpoint] = job.soln

    return job
end
