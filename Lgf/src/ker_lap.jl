using QuadGK

function ker_lap(i,j )
    fnc = t ->( real((1 - ( (t-sqrt(im))./(t+sqrt(im)) ).^(j+abs(i)) .* ( (t+sqrt(-im))./(t - sqrt(-im)) ).^(j-abs(i)) )) ./t )
    v,~ = quadgk(fnc,1e-14,1.0,reltol=sqrt(eps))
    v = v/2/pi
    return v
end
