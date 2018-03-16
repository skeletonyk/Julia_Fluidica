mutable struct Depot
    u :: Array{Float64,2}
    v :: Array{Float64,2}
    u_buffer :: Array{Float64,2}
    v_buffer :: Array{Float64,2}
    ufr :: Array{Float64,2}
    vfr :: Array{Float64,2}
    w_buffer :: Array{Float64,2}
    w_lap :: Array{Float64,2}
    t1 :: Array{Float64,2}
    t2 :: Array{Float64,2}
    s_buffer :: Array{Float64,2}
    s_buffer_2 :: Array{Complex{Float64},2}
    s_buffer_3 :: Array{Complex{Float64},2}
    w_buffer_small :: Array{Float64,2}
    s_buffer_small :: Array{Float64,2}
    s_buffer_small_2 ::Array{Complex{Float64},2}
    s_buffer_small_3 :: Array{Complex{Float64},2}
    wx_out :: Array{Float64,2}
    wy_out :: Array{Float64,2}
    avg_u :: Array{Float64,2}
    avg_v :: Array{Float64,2}
    Depot() = new()
end

function depot_ini(m :: Int64, n :: Int64)
    bf_dpt = Depot()

    bf_dpt.w_buffer = zeros(Float64, 2*m+3, 2*n+3)
    bf_dpt.t1 = zeros(Float64,m,n)
    bf_dpt.t2 = zeros(Float64,m,n)
    bf_dpt.w_lap = similar(bf_dpt.t1)
    bf_dpt.s_buffer = zeros(Float64, 2*m+3, 2*n+3)
    bf_dpt.s_buffer_2 = rfft(bf_dpt.s_buffer)
    bf_dpt.s_buffer_3 = similar(bf_dpt.s_buffer_2)
    bf_dpt.w_buffer_small = zeros(Float64, 2*m-1, 2*n-1)
    bf_dpt.s_buffer_small = similar(bf_dpt.w_buffer_small) #zeros(Float64, 2*m-1, 2*n-1)
    bf_dpt.s_buffer_small_2 = rfft(bf_dpt.s_buffer_small)
    bf_dpt.s_buffer_small_3 = similar(bf_dpt.s_buffer_small_2)

    bf_dpt.u = zeros(Float64, m+2, n+1)
    bf_dpt.v = zeros(Float64, m+1, n+2)
    bf_dpt.ufr = zeros(Float64, m+2, n+1)
    bf_dpt.vfr = zeros(Float64, m+1, n+2)
    bf_dpt.u_buffer = similar(bf_dpt.ufr)
    bf_dpt.v_buffer = similar(bf_dpt.vfr)

    bf_dpt.wx_out = zeros(m, n+1)
    bf_dpt.wy_out = zeros(m+1, n)
    bf_dpt.avg_u = zeros(Float64, m+1, n)
    bf_dpt.avg_v = zeros(Float64, m, n+1)


    return bf_dpt
end
