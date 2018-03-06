module Ddfs
export yang3

function yang3(r)
    f = zeros(size(r))
    rr = abs.(r)

    pts = find( x -> x<=1 , rr)
    f[pts] = 17/48 + sqrt(3) * pi/108 + rr[pts]/4 - rr[pts].^2/4 + (1-2*rr[pts])/16 .* sqrt.(-12*rr[pts].^2+12*rr[pts]+1) -
    sqrt(3)/12 .* asin.(sqrt(3)/2*(2*rr[pts]-1))

    pts = find(x -> x > 1 && x <= 2, rr )
    f[pts] = 55/48-sqrt(3)*pi/108-13*rr[pts]/12+rr[pts].^2/4+(2*rr[pts]-3)/48.*sqrt.(-12*rr[pts].^2+36*rr[pts]-23) +
    sqrt(3)/36 * asin.(sqrt(3)/2*(2*rr[pts]-3))

    return f
end
end
