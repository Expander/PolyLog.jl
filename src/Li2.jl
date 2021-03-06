"""
    li2(x::Float64)::Float64

Returns the real dilogarithm of a real number `x` of type `Float64`.

Author: Alexander Voigt

License: MIT

# Example
```julia
li2(1.0)
```
"""
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
        (x/(x - 1.0), -0.5*log1p(-x)^2, -1.0)
    elseif x == 0.0
        return 0.0
    elseif x < 0.5
        (x, 0.0, 1.0)
    elseif x < 1.0
        (1.0 - x, pi*pi/6.0 - log(x)*log(1.0 - x), -1.0)
    elseif x == 1.0
        return pi*pi/6.0
    elseif x < 2.0
        l = log(x)
        (1.0 - 1.0/x, pi*pi/6.0 - l*(log(1.0 - 1.0/x) + 0.5*l), 1.0)
    else
        (1.0/x, pi*pi/3.0 - 0.5*log(x)^2, -1.0)
    end

    z = y - 0.25
    z2 = z*z
    z4 = z2*z2

    p = cp[1] + z * cp[2] + z2 * (cp[3] + z * cp[4]) +
        z4 * (cp[5] + z * cp[6] + z2 * (cp[7] + z * cp[8]))
    q = cq[1] + z * cq[2] + z2 * (cq[3] + z * cq[4]) +
        z4 * (cq[5] + z * cq[6] + z2 * (cq[7] + z * cq[8]))

    r + s*y*p/q
end


"""
    li2(z::ComplexF64)::ComplexF64

Returns the complex dilogarithm of a complex number `z` of type `ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```julia
li2(1.0 + 1.0im)
```
"""
function li2(z::ComplexF64)::ComplexF64
    clog(z) = 0.5*log(abs2(z)) + angle(z)*1.0im

    rz = real(z)
    iz = imag(z)

    if iz == 0.0
        if rz <= 1.0
            return li2(rz)
        else # rz > 1.
            return li2(rz) - pi*log(rz)*1.0im
        end
    end

    nz = abs2(z)

    if nz < eps(Float64)
        return z*(1.0 + 0.25*z)
    end

    (u::ComplexF64, rest::ComplexF64, sgn::Float64) = if rz <= 0.5
        if nz > 1.0
            (-clog(1.0 - 1.0 / z), -0.5 * clog(-z)^2 - pi^2 / 6.0, -1.0)
        else # nz <= 1.
            (-clog(1.0 - z), 0.0 + 0.0im, 1.0)
        end
    else # rz > 0.5
        if nz <= 2.0*rz
            l = -clog(z)
            (l, l * clog(1.0 - z) + pi^2 / 6.0, -1.0)
        else # nz > 2.0*rz
            (-clog(1.0 - 1.0 / z), -0.5 * clog(-z)^2 - pi^2 / 6.0, -1.0)
        end
    end

    bf = (
        - 1.0/4.0,
          1.0/36.0,
        - 1.0/3600.0,
          1.0/211680.0,
        - 1.0/10886400.0,
          1.0/526901760.0,
        - 4.0647616451442255e-11,
          8.9216910204564526e-13,
        - 1.9939295860721076e-14,
          4.5189800296199182e-16
    )

    u2::ComplexF64 = u*u
    u4::ComplexF64 = u2*u2

    rest + sgn * (
        u +
        u2 * (bf[1] +
        u  * (bf[2] +
        u2 * (
            bf[3] +
            u2*bf[4] +
            u4*(bf[5] + u2*bf[6]) +
            u4*u4*(bf[7] + u2*bf[8] + u4*(bf[9] + u2*bf[10]))
        )))
    )
end
