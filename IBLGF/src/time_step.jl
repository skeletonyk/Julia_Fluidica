function time_step(parms, job, solver)
#
    # Inital soln
    soln = Soln()
    soln.t = parms.hzn[1];
    soln.w = zeros(job.dims[1],job.dims[2])
    soln.u = zeros(job.dims[1]+2,job.dims[2]+1)
    soln.v = zeros(job.dims[1]+1,job.dims[2]+2)
    soln.r1 = zeros(job.dims[1],job.dims[2])
    soln.t = 0
    soln.it = 0

    first::Bool = true
    tic()
    for k=1:100#job.nstp
        #w = copy(soln.w)
        #t = soln.t

        r1old = copy(soln.r1)
        #@time #soln.w, soln.t, soln.f, soln.r1  =
        @time solver.step(soln, solver, first)
        soln.it = soln.it + 1
        #cfl = solver.cfl(soln)

        #if cfl > 1.5
        #    error("The CFL has become too large and Fluidica cannot continue().  Consider choosing a lower time step so that CFL < 1.5")
        #    return
        #end
        first = false
    end
    toc()
    #job.runtime = toc()
    #soln.checkpoint = soln.checkpoint+1
    #job.data[job.checkpoint] = job.soln

    return soln
end
