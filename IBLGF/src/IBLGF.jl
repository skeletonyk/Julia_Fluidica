module IBLGF
    using Util
    using FFTW
    using Ddfs
    using Lgf

    import Util.meshgrid

    mutable struct Soln
        f :: Array{Float64,2}
        w :: Array{Float64,2}
        r1old :: Array{Float64,2}
        t :: Float64
        u :: Array{Float64,2}
        v :: Array{Float64,2}
        Soln() = new()
    end

    mutable struct Job
        ddf :: Function
        pts :: Array{Float64,2}
        dims :: Vector
        domain
        center
        x
        y
        r
        bdy :: Array{Float64,2}
        soln :: Soln
        dt :: Float64
        nstp ::Int64
        it
        times
        g :: Array{Float64,2}
        ghat :: Array{Complex64,2}
        ghatbig :: Array{Complex64,2}
        ghat_p :: Array{Complex64,2}
        etx :: AbstractSparseMatrix
        ety :: AbstractSparseMatrix
        ec  :: AbstractSparseMatrix
        t   :: Array{Float64,2}
        z   :: Array{Float64,2}
        frame_rot :: Function
        frame_xvel :: Function
        frame_yvel :: Function
        cno
        runtime
        checkpoint
        Job() = new()
    end

    mutable struct Solver
        dt :: Float64
        g_times :: Function
        tinv_times :: Function
        ainv_times :: Function
        la_inv_times :: Function
        zinv_times :: Function
        visc :: Function
        b_times :: Function
        bt_times :: Function
        r1 :: Function
        r2 :: Function
        step :: Function
        speed :: Function
        cfl :: Function
        vel_a :: Function
        Solver() = new()
    end

    FFTW_Type = Base.DFT.FFTW.rFFTWPlan{Float64,-1,false,2}
    IFFTW_Type = Base.DFT.ScaledPlan{Complex{Float64},Base.DFT.FFTW.rFFTWPlan{Complex{Float64},1,false,2},Float64}

    include("et_mat.jl")
    include("ec_mat.jl")
    include("update_setup.jl")
    include("init_ns.jl")
    include("differential_operators.jl")
    include("velocity.jl")
    include("time_step.jl")
    include("cn_ab2.jl")

    export update_setup
    export init_ns
    export time_step
    export cn_ab2
    export Job
    export Soln
    export Solver


end
