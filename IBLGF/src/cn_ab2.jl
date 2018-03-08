function cn_ab2(soln:: Soln, h :: Solver_static, first::Bool)
# Split Crank-Nicholson Adams-Bashforth scheme for DAE (including nested soln.
# of KKT system())
#
# the equations to be integrated are
#
# dw/dt + B' f = r1[w,t] + nu L w       [1]
#   B L^{-1} w = r2[t]                  (2)
#
# s is the solution with components w,f, and t.  w and f on input are
# the differential and algebraic variables, respectively, at time t.
#
# snew is the solution with same components at time t + h.dt
#
# nu L w is to be integrated with CN
# r1 is to be integrated with AB2
# the B' f term uses implicit Euler and yields a first-order-accurate f.
# Eq [2] is enforced at each time level.

# h is a structure with function handles for the various operators
# appearing in [1] and [2].  It is setup outside this routine, and needs
# to be refreshed whenever h.dt changes.  The needed function handles are:
#
# h.r1[w,t]         # computes first part of rhs of [1]
# h.r2[t]           # computes rhs of [2]
# h.b_times[w]      # applies B to data size(w)
# h.bt_times[f]     # applies B' to fdata size(f)
# h.g_times[w]      # applies L^{-1} to data size(w)
# h.zinv_times[f]   # applies Z^{-1}, where Z = B L^{-1} A^{-1} B'.  to data size(f)
# h.ainv_times[w]   # applies A^{-1}, where A = I - 0.5*h.dt * nu L

# note that the whole trick is to do have function handles that
# apply these operators as efficiently as possible!!!

    #snew = Soln()
    #snew.u = s.u;
    #snew.v = s.v;
    #w :: Array{Float64,2} = (s.w)
    #t :: Float64 = s.t

    #print("non-linear ->")
    #@time
    r1 :: Array{Float64,2} = h.r1(soln.w, soln.t);  # compute current nonlinear term

    rhs :: Array{Float64,2} = similar(r1)

    if !first
        rhs = h.dt.*(1.5 .* r1 - 0.5 .* soln.r1);  # add old nonlinear term
    else
        rhs = h.dt.*r1;  # use Euler step for first timestep
    end
    #snew.r1old = r1;  # save current nonlinear term for next step
    rhs += soln.w + h.visc(soln.w)

    #print("a_inv      ->")
    #@time
    wstar :: Array{Float64,2} = h.ainv_times(rhs)
    #print("g_times    ->")
    #@time
    tmp :: Array{Float64,2} = h.b_times(h.g_times(wstar))-h.r2(soln.t)
    #print("z_inv      ->")
    #@time
    soln.f = h.zinv_times(tmp )./h.dt

    soln.w = wstar - h.ainv_times( h.dt*h.bt_times(soln.f ))
    soln.t = soln.t + h.dt

    #w,t,f,r1

end
