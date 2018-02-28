function curltx(u :: Array{Float64,2})
    # - d/dy
    # 1/12/17 [TC] changed from original version to take array argument

    m = size(u,1)
    n = size(u,2)

    return (-u[1:m,2:n]+u[1:m,1:n-1])
end

function curlty(v :: Array{Float64,2})
    #  d/dx
    # 1/12/17 [TC] changed from original version to take array argument

    m = size(v,1)
    n = size(v,2)

    return (v[2:m,1:n]-v[1:m-1,1:n])
end

function curl(s :: Array{Float64,2})
    # given a streamfunction field, s, compute u,v [on cell edges]
    # s[M,N] --> ( u[M,N-1],v[M-1,N]
    m = size(s,1)
    n = size(s,2)

    u :: Array{Float64,2} = s[1:m,2:n]-s[1:m,1:n-1]
    v :: Array{Float64,2} = s[1:m-1,1:n]-s[2:m,1:n]

    u, v
end
function lap(f :: Array{Float64,2})
    # discrete Laplacian with zero Dirichlet BC

    m = size(f,1)
    n = size(f,2)

    fext = [ zeros(1,n+2); [ zeros(m,1) f zeros(m,1)]; zeros(1,n+2)]
    a = -4*fext[2:m+1,2:n+1] +
            fext[1:m,2:n+1]+fext[3:m+2,2:n+1] +
            fext[2:m+1,1:n]+fext[2:m+1,3:n+2]

    return a
end
