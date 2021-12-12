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
        0.9999999999999999502e+0,
       -2.6883926818565423430e+0,
        2.6477222699473109692e+0,
       -1.1538559607887416355e+0,
        2.0886077795020607837e-1,
       -1.0859777134152463084e-2
    )
    cq = (
        1.0000000000000000000e+0,
       -2.9383926818565635485e+0,
        3.2712093293018635389e+0,
       -1.7076702173954289421e+0,
        4.1596017228400603836e-1,
       -3.9801343754084482956e-2,
        8.2743668974466659035e-4
    )

    # transform to [0, 1/2]
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

    y2 = y*y
    y4 = y2*y2

    p = cp[1] + y * cp[2] + y2 * (cp[3] + y * cp[4]) +
        y4 * (cp[5] + y * cp[6])
    q = cq[1] + y * cq[2] + y2 * (cq[3] + y * cq[4]) +
        y4 * (cq[5] + y * cq[6] + y2 * cq[7])

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
