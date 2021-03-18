function horner(x, coeffs)
    s = coeffs[end]
    for k in length(coeffs)-1:-1:1
        s = coeffs[k] + x * s
    end
    return s
end


# Returns the real dilogarithm of a real number of type `Float64`.
# Author: Alexander Voigt
# License: MIT
function li2(x::Float64)::Float64
    cp = (
        1.0706105563309304277e+0,
       -4.5353562730201404017e+0,
        7.4819657596286408905e+0,
       -6.0516124315132409155e+0,
        2.4733515209909815443e+0,
       -4.6937565143754629578e-1,
        3.1608910440687221695e-2,
       -2.4630612614645039828e-4
    )
    cq = (
        1.0000000000000000000e+0,
       -4.5355682121856044935e+0,
        8.1790029773247428573e+0,
       -7.4634190853767468810e+0,
        3.6245392503925290187e+0,
       -8.9936784740041174897e-1,
        9.8554565816757007266e-2,
       -3.2116618742475189569e-3
    )

    # transform to [0, 1/2)
    (y, r, s) = if x < -1.0
        l = log(1.0 - x)
        (1.0/(1.0 - x), -pi*pi/6.0 + l*(0.5*l - log(-x)), 1.0)
    elseif x == -1.0
        return -pi*pi/12.
    elseif x < 0.0
        l = log1p(-x)
        (x/(x - 1.0), -0.5*l*l, -1.0)
    elseif x == 0.0
        return 0.0
    elseif x < 0.5
        (x, 0.0, 1.0)
    elseif x < 1.
        (1.0 - x, pi*pi/6.0 - log(x)*log(1.0 - x), -1.0)
    elseif x == 1.
        return pi*pi/6.0
    elseif x < 2.
        l = log(x)
        (1.0 - 1.0/x, pi*pi/6.0 - l*(log(1.0 - 1.0/x) + 0.5*l), 1.0)
    else
        l = log(x)
        (1.0/x, pi*pi/3.0 - 0.5*l*l, -1.0)
    end

    z = y - 0.25
    p = horner(z, cp)
    q = horner(z, cq)

    return r + s*y*p/q
end
