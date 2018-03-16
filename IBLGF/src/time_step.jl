function time_step(parms, job, solver)
#
    # Inital soln
    soln = Soln()
    soln.t = parms.hzn[1];
    soln.w = zeros(Float64, job.dims[1],job.dims[2])
    soln.u = zeros(Float64, job.dims[1]+2,job.dims[2]+1)
    soln.v = zeros(Float64, job.dims[1]+1,job.dims[2]+2)
    soln.r1 = zeros(Float64, job.dims[1],job.dims[2])
    soln.t = 0.0
    soln.it = 0
    r1_buffer = similar(soln.w)
    rhs_buffer = similar(soln.w)
    wstar_buffer = similar(soln.w)
    L_inv_buffer = similar(soln.w)
    L_inv_buffer_big = zeros(Float64, job.dims[1]+2,job.dims[2]+2)

    first::Bool = true

    tic()
    for k=1:100#job.nstp
        #w = copy(soln.w)
        #t = soln.t

        r1old = copy(soln.r1)
        #@time #soln.w, soln.t, soln.f, soln.r1  =
        solver.step(soln, solver, first, r1_buffer, rhs_buffer, wstar_buffer, L_inv_buffer, L_inv_buffer_big)
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
