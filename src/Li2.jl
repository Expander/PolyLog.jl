"""
    reli2(x::Real)

Returns the real dilogarithm ``\\Re[\\operatorname{Li}_2(x)]`` of a
real number ``x`` of type `Real`.

Implemented as rational function approximation with a maximum error of
`5e-17` [[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> reli2(1.0)
1.6449340668482264
```
"""
reli2(x::Real) = _reli2(float(x))

_reli2(x::Float16) = oftype(x, _reli2(Float32(x)))

_reli2(x::Float32) = oftype(x, _reli2(Float64(x)))

function _reli2(x::Float64)::Float64
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
    (y, rest, sgn) = if x < -1.0
        l = log(1.0 - x)
        (1.0/(1.0 - x), -zeta2 + l*(0.5*l - log(-x)), 1.0)
    elseif x == -1.0
        return -0.5*zeta2
    elseif x < 0.0
        (x/(x - 1.0), -0.5*log1p(-x)^2, -1.0)
    elseif x == 0.0
        return 0.0
    elseif x < 0.5
        (x, 0.0, 1.0)
    elseif x < 1.0
        (1.0 - x, zeta2 - log(x)*log(1.0 - x), -1.0)
    elseif x == 1.0
        return zeta2
    elseif x < 2.0
        l = log(x)
        (1.0 - 1.0/x, zeta2 - l*(log(1.0 - 1.0/x) + 0.5*l), 1.0)
    else
        (1.0/x, 2.0*zeta2 - 0.5*log(x)^2, -1.0)
    end

    y2 = y*y
    y4 = y2*y2

    p = cp[1] + y * cp[2] + y2 * (cp[3] + y * cp[4]) +
        y4 * (cp[5] + y * cp[6])
    q = cq[1] + y * cq[2] + y2 * (cq[3] + y * cq[4]) +
        y4 * (cq[5] + y * cq[6] + y2 * cq[7])

    rest + sgn*y*p/q
end


"""
    li2(z::Complex)

Returns the complex dilogarithm ``\\operatorname{Li}_2(z)`` of a
complex number ``z`` of type `Complex`.

If only real arguments ``z\\in\\mathbb{R}`` are considered and one is
interested only in the real part of the dilogarithm,
``\\Re[\\operatorname{Li}_2(z)]``, refer to the function
[`reli2`](@ref), which may be a faster alternative.

See also [`reli2`](@ref).

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li2(1.0 + 1.0im)
0.616850275068085 + 1.4603621167531196im
```
"""
li2(z::Complex) = _li2(float(z))

li2(z::Real) = li2(Complex(z))

_li2(z::ComplexF16) = oftype(z, _li2(ComplexF32(z)))

_li2(z::ComplexF32) = oftype(z, _li2(ComplexF64(z)))

function _li2(z::ComplexF64)::ComplexF64
    clog(z) = 0.5*log(abs2(z)) + angle(z)*1.0im

    rz = real(z)
    iz = imag(z)

    if iz == 0.0
        if rz <= 1.0
            return reli2(rz)
        else # rz > 1.
            return reli2(rz) - pi*log(rz)*1.0im
        end
    end

    nz = abs2(z)

    if nz < eps(Float64)
        return z*(1.0 + 0.25*z)
    end

    (u, rest, sgn) = if rz <= 0.5
        if nz > 1.0
            (-clog(1.0 - inv(z)), -0.5*clog(-z)^2 - zeta2, -1.0)
        else # nz <= 1.
            (-clog(1.0 - z), 0.0 + 0.0im, 1.0)
        end
    else # rz > 0.5
        if nz <= 2.0*rz
            l = -clog(z)
            (l, l*clog(1.0 - z) + zeta2, -1.0)
        else # nz > 2.0*rz
            (-clog(1.0 - inv(z)), -0.5*clog(-z)^2 - zeta2, -1.0)
        end
    end

    B = (
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

    u2 = u*u
    u4 = u2*u2

    rest + sgn*(
        u +
        u2*(B[1] +
        u *(B[2] +
        u2*(
            B[3] +
            u2*B[4] +
            u4*(B[5] + u2*B[6]) +
            u4*u4*(B[7] + u2*B[8] + u4*(B[9] + u2*B[10]))
        )))
    )
end
