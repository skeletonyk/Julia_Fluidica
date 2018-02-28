function ec_mat(mn, etx, ety)

    # altered 1/12/17 [TC].
    # changed functionality of curltx and curlty so they can be reused in
    # nonline terms

    nb = size(etx,2)
    ctxetx = spzeros(mn[1]*mn[2],nb)
    ctyety = spzeros(mn[1]*mn[2],nb)

    for k=1:nb
        ctxetx[:,k] = reshape( curltx(reshape(full(etx[:,k]),mn[1],mn[2]+1)), prod(mn), 1)
        ctyety[:,k] = reshape( curlty(reshape(full(ety[:,k]),mn[1]+1,mn[2])), prod(mn), 1)
    end

    return [ctxetx' ; ctyety']
end
