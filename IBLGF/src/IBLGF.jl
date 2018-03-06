module IBLGF
    using Util
    using Devectorize
    using FFTW
    using Ddfs
    using Lgf

    import Util.meshgrid

    mutable struct Soln
        f :: Array{Float64,2}
        w :: Array{Float64,2}
        r1 :: Array{Float64,2}
        t :: Float64
        u :: Array{Float64,2}
        v :: Array{Float64,2}
        it :: Int64
        Soln() = new()
    end

    struct Job_static
        ddf :: Function
        pts :: Array{Float64,2}
        dims :: Array{Int64,1}
        domain :: Array{Float64,2}
        center :: Array{Float64,2}
        x  :: Array{Float64,2}
        y :: Array{Float64,2}
        r :: Float64
        bdy :: Array{Float64,2}
        dt :: Float64
        nstp ::Int64
        times :: AbstractVector
        g :: Array{Float64,2}
        ghat :: Array{Complex64,2}
        ghatbig :: Array{Complex64,2}
        ghat_p :: Array{Complex64,2}
        etx :: AbstractSparseMatrix
        ety :: AbstractSparseMatrix
        ec  :: AbstractSparseMatrix
        #t   :: Array{Float64,2}
        #z   ::Base.LinAlg.BunchKaufman{Float64,Array{Float64,2}}
        frame_rot :: Function
        frame_xvel :: Function
        frame_yvel :: Function
        cno :: Float64
        #runtime
        checkpoint :: Int64
    end
    mutable struct Job
        ddf :: Function
        pts :: Array{Float64,2}
        dims :: Array{Int64,1}
        domain :: Array{Float64,2}
        center :: Array{Float64,2}
        x  :: Array{Float64,2}
        y :: Array{Float64,2}
        r :: Float64
        bdy :: Array{Float64,2}
        dt :: Float64
        nstp ::Int64
        times :: AbstractVector
        g :: Array{Float64,2}
        ghat :: Array{Complex64,2}
        ghatbig :: Array{Complex64,2}
        ghat_p :: Array{Complex64,2}
        etx :: AbstractSparseMatrix
        ety :: AbstractSparseMatrix
        ec  :: AbstractSparseMatrix
        #t   :: Array{Float64,2}
        #z   ::Base.LinAlg.BunchKaufman{Float64,Array{Float64,2}}
        frame_rot :: Function
        frame_xvel :: Function
        frame_yvel :: Function
        cno :: Float64
        #runtime
        checkpoint :: Int64
        Job() = new()
    end


    mutable struct Solver
        dt :: Float64
        g_times :: Function
        #tinv_times :: Function
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

    struct Solver_static
        dt :: Float64
        g_times :: Function
        #tinv_times :: Function
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
    export Solver_static


end
