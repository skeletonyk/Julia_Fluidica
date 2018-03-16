function curltx(u :: Array{Float64,2}, u_out ::Array{Float64,2})
    # - d/dy
    # 1/12/17 [TC] changed from original version to take array argument

    #m = size(u,1)
    #n = size(u,2)
    for i=1:size(u,1), j=1:size(u,2)-1
        u_out[i,j] = -u[i,j+1] + u[i,j]
    end
    nothing
    #return (-u[1:m,2:n]+u[1:m,1:n-1])
end

function curlty(v :: Array{Float64,2}, v_out :: Array{Float64,2})
    #  d/dx
    # 1/12/17 [TC] changed from original version to take array argument

    #m = size(v,1)
    #n = size(v,2)
    for i = 1:size(v,1)-1, j=1:size(v,2)
        v_out[i,j] = v[i + 1, j] - v[i, j]
    end
    nothing

    #return (v[2:m,1:n]-v[1:m-1,1:n])
end

function curl(s :: Array{Float64,2}, u_out :: Array{Float64,2}, v_out :: Array{Float64,2})
    # given a streamfunction field, s, compute u,v [on cell edges]
    for i=1:size(s,1), j=1:size(s,2)-1
        u_out[i,j] = s[i,j+1] - s[i,j]
    end
    for i=1:size(s,1)-1, j=1:size(s,2)
        v_out[i,j] = -s[i+1,j] + s[i,j]
    end
    nothing
end
function lap(f :: Array{Float64,2}, f_out :: Array{Float64,2})
    # discrete Laplacian with zero Dirichlet BC

    m = size(f,1)
    n = size(f,2)

    fext :: Array{Float64,2} = [ zeros(1,n+2); [ zeros(m,1) f zeros(m,1)]; zeros(1,n+2)]

    [-4*fext[i+1,j+1] + fext[i,j+1] + fext[i+2,j+1] + fext[i+1,j] + fext[i+1,j+2] for i = 1 : m, j = 1 : n]

    for i = 2 : m-1, j = 2 : n-1
        f_out[i,j] = -4*f[i,j] + f[i-1,j] + fext[i+1,j] + fext[i,j-1] + fext[i,j+1]
    end

    for j = 2 : n-1
        f_out[1,j] = -4*f[1,j] + f[2,j] + fext[i,j-1] + fext[i,j+1]
    end

    #-4*fext[2:m+1,2:n+1] + fext[1:m,2:n+1]+fext[3:m+2,2:n+1] + fext[2:m+1,1:n]+fext[2:m+1,3:n+2]

end
