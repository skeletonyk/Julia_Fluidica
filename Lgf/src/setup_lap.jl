using MAT

function setup_lap(mn)
    m = maximum(mn)
    n = minimum(mn)
    g = zeros(m,n)
    if m <= 1000 && n <= 1000
        # use precomputed lookup table to save time
        #disp("Using precomputed table for LGF")
        fname = "Lgf/src/lgf_lap_table.mat"
        if isfile(fname)
            file = matopen(fname)
            gsave = read(file, "gsave")
            @inbounds for i=1:m
                @inbounds for j=1:minimum([i,n])
                    g[i,j] = gsave[i,j]
                end
            end
            @inbounds for i=1:n
                @inbounds for j=i+1:n
                    g[i,j] = g[j,i]
                end
            end
        else
            error("Could not find precomputed lookup table for Lattice Greens Function")
        end
    else   # build g
        println("Building table for LGF")
        for i=1:m
            @fastmath @inbounds @simd for j=1:minimum([i,n])
                g[i,j] = ker_lap(i-1,j-1)
            end
        end
        for i=1:n
            @fastmath @inbounds @simd for j=i+1:n
                g[i,j] = g[j,i]
            end
        end
    end

    # transpose if necessary
    if (mn[1]< mn[2])
        g = g'
    end
    return g
end

function transform_lgf(lgf :: Array{Float64,2})
  # padded fft of LGF
  lgf_large =
  [ flipdim(flipdim(lgf[2:end,2:end],1),2) flipdim(lgf[2:end,:],1);
    flipdim(lgf[:,2:end],2) lgf ]
  # Tranform augmented table
  tmp = rfft(lgf_large)
  return tmp
end
