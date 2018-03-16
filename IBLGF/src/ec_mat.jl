function ec_mat(mn, etx, ety)

    # altered 1/12/17 [TC].
    # changed functionality of curltx and curlty so they can be reused in
    # nonline terms

    nb = size(etx,2)
    ctxetx = spzeros(mn[1]*mn[2],nb)
    ctyety = spzeros(mn[1]*mn[2],nb)
    tmp = zeros(mn[1], mn[2])

    for k=1:nb
        curltx(reshape(full(etx[:,k]),mn[1],mn[2]+1),tmp)
        ctxetx[:,k] = reshape(tmp, prod(mn), 1)

        curlty(reshape(full(ety[:,k]),mn[1]+1,mn[2]),tmp)
        ctyety[:,k] = reshape(tmp, prod(mn), 1)
    end

    return [ctxetx' ; ctyety']
end
