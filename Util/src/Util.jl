module Util
    using Ddfs
    using Base.repeat

    export create_circle
    export meshgrid
    export get_dt
    export Geo

    include("create_object.jl")

    # Types


    mutable struct Geo
        xyorig :: Array{Float64,2}
        xy :: Array{Float64,2}
        Geo() = new()
    end

    # Simple ones ------------------------------
    function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
        m, n = length(vy), length(vx)

        gx :: Array{Float64,2} = reshape(repeat(vx, inner = m, outer = 1), m, n)
        gy :: Array{Float64,2} = reshape(repeat(vy, inner = 1, outer = n), m, n)

        return gx, gy
    end

    function get_dt(parms)
        if parms.tstep.auto == 1
            # try to find a good time step
            vmax = abs(parms.axlim[1]) * sqrt(2) #analyze_mot(parms)
            dt = parms.tstep.cfl/vmax
        else
            dt = parms.tstep.cfl
        end
        println(dt)
        nstp = max(1,floor(parms.hzn[3]*parms.ppl / dt)); # want the time increment to be an integer no. of time steps
        dt = parms.hzn[3]*parms.ppl/nstp
        println(floor(parms.hzn[3]*parms.ppl /dt))
        nstp = convert(Int64, floor(parms.hzn[3]*parms.ppl /dt))

        dt, nstp
    end

# Simple ones ------------------------------
end
