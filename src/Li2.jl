# rational function approximation of Re[Li2(x)] for x in [0, 1/2]
function reli2_approx(x::Float32)::Float32
    cp = (1.00000020f0, -0.780790946f0, 0.0648256871f0)
    cq = (1.00000000f0, -1.03077545f0, 0.211216710f0)

    x2 = x*x

    p = cp[1] + x * cp[2] + x2 * cp[3]
    q = cq[1] + x * cq[2] + x2 * cq[3]

    x*p/q
end


# rational function approximation of Re[Li2(x)] for x in [0, 1/2]
function reli2_approx(x::Float64)::Float64
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

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    x*p/q
end


# approximation of Re[Li2(x)] for x in [0, 1/2]
function reli2_approx(x::BigFloat)::BigFloat
    sum = x
    xn = x*x

    for k in 2:typemax(Int64)
        term = xn/oftype(x, k)^2
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        xn *= x
    end

    sum
end


# series expansion of Li2(z) for |z| <= 1 and Re(z) <= 0.5
# in terms of u = -log(1-z)
function li2_approx(u::ComplexF32)::ComplexF32
    B = (-1.0f0/4, 1.0f0/36, -1.0f0/3600, 1.0f0/211680)
    u2 = u*u
    u + u2*(B[1] + u*(B[2] + u2*(B[3] + u2*B[4])))
end


# series expansion of Li2(z) for |z| <= 1 and Re(z) <= 0.5
# in terms of u = -log(1-z)
function li2_approx(u::ComplexF64)::ComplexF64
    B = (
        - 1.0/4,
          1.0/36,
        - 1.0/3600,
          1.0/211680,
        - 1.0/10886400,
          1.0/526901760,
        - 4.0647616451442255e-11,
          8.9216910204564526e-13,
        - 1.9939295860721076e-14,
          4.5189800296199182e-16
    )

    u2 = u*u
    u4 = u2*u2

    u +
    u2*(B[1] +
    u *(B[2] +
    u2*(
        B[3] +
        u2*B[4] +
        u4*(B[5] + u2*B[6]) +
        u4*u4*(B[7] + u2*B[8] + u4*(B[9] + u2*B[10]))
    )))
end


# series expansion of Li2(z) for |z| <= 1 and Re(z) <= 0.5
function li2_approx(z::Complex{T})::Complex{T} where T
    if abs2(z) < (99/100)^2
        li2_approx_taylor(z)
    else
        li2_approx_unity(-clog1p(-z))
    end
end


# Taylor series expansion of Li2(z) for |z| < 1
function li2_approx_taylor(z::Complex{T})::Complex{T} where T
    sum = z
    zn = z*z

    for k in 2:typemax(Int64)
        term = zn/convert(T, k)^2
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        zn *= z
    end

    sum
end


# series expansion of Li2(z) for |z| <= 1 and Re(z) <= 0.5
# in terms of u = -log(1-z)
function li2_approx_unity(u::Complex{T})::Complex{T} where T
    u2 = u*u
    c0 = inv(4*big(pi)^2)
    p = u2*c0^2
    sum = one(T)/72

    for n in 2:typemax(Int64)
        old_sum = sum
        sgn = iseven(n) ? -1 : 1
        sum += sgn*p/(2*n + 1)*zeta(2*n, T)
        sum == old_sum && break
        p *= u2*c0
    end

    u + u2*(-one(T)/4 + 2*u*sum)
end


"""
    reli2(x::Real)

Returns the real dilogarithm ``\\Re[\\operatorname{Li}_2(x)]`` of a
real number ``x`` of type `Real`.

For ``x`` of type `Float16`, `Float32` or `Float64` the implementation
is an adaptation of
[[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> reli2(1.0)
1.6449340668482264

julia> reli2(BigFloat("1.0"))
1.644934066848226436472415166646025189218949901206798437735558229370007470403202
```
"""
reli2(x::Real) = _reli2(float(x))

_reli2(x::Float16) = oftype(x, _reli2(Float32(x)))

function _reli2(x::T)::T where T
    # transform to [0, 1/2]
    if x < -one(x)
        l = log(one(x) - x)
        reli2_approx(inv(one(x) - x)) - zeta2(typeof(x)) + l*(one(x)/2*l - log(-x))
    elseif x == -one(x)
        -one(x)/2*zeta2(typeof(x))
    elseif x < zero(x)
        -reli2_approx(x/(x - one(x))) - one(x)/2*log1p(-x)^2
    elseif iszero(x)
        zero(x)
    elseif x < one(x)/2
        reli2_approx(x)
    elseif x == one(x)/2
        oftype(x, pi)^2/12 - log(oftype(x, 2))^2/2
    elseif x < one(x)
        -reli2_approx(one(x) - x) + zeta2(typeof(x)) - log(x)*log1p(-x)
    elseif x == one(x)
        zeta2(typeof(x))
    elseif x < 2*one(x)
        l = log(x)
        reli2_approx(one(x) - inv(x)) + zeta2(typeof(x)) - l*(log(one(x) - inv(x)) + one(x)/2*l)
    else
        -reli2_approx(inv(x)) + 2*zeta2(typeof(x)) - one(x)/2*log(x)^2
    end
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
0.6168502750680851 + 1.4603621167531196im

julia> li2(BigFloat(1) + 1im)
0.6168502750680849136771556874922594459571062129525494141508343360137528014012052 + 1.460362116753119547679775739491787597608795299373993707847946932920340157070418im
```
"""
li2(z::Complex) = _li2(float(z))

li2(z::Real) = li2(Complex(z))

_li2(z::ComplexF16) = oftype(z, _li2(ComplexF32(z)))

function _li2(z::Complex{T})::Complex{T} where T
    rz, iz = reim(z)

    if iszero(iz)
        if rz <= one(T)
            complex(reli2(rz))
        else # Re(z) > 1
            complex(reli2(rz), -convert(T, pi)*log(rz))
        end
    else
        nz = abs2(z)

        if nz < eps(T)
            z*(one(T) + one(T)/4*z)
        else
            if rz <= one(T)/2
                if nz > one(T)
                    -li2_approx(-clog1p(-inv(z))) - one(T)/2*clog(-z)^2 - zeta2(T)
                else # |z|^2 <= 1
                    li2_approx(-clog1p(-z))
                end
            else # Re(z) > 1/2
                if nz <= 2*rz
                    l = -clog(z)
                    -li2_approx(l) + l*clog1p(-z) + zeta2(T)
                else # |z|^2 > 2*Re(z)
                    -li2_approx(-clog1p(-inv(z))) - one(T)/2*clog(-z)^2 - zeta2(T)
                end
            end
        end
    end
end

function _li2(z::Complex{BigFloat})::Complex{BigFloat}
    rz, iz = reim(z)

    if iszero(iz)
        if rz <= one(BigFloat)
            complex(reli2(rz))
        else # Re(z) > 1
            complex(reli2(rz), -convert(BigFloat, pi)*log(rz))
        end
    else
        nz = abs2(z)

        if nz < eps(BigFloat)
            z*(one(BigFloat) + one(BigFloat)/4*z)
        else
            if rz <= one(BigFloat)/2
                if nz > one(BigFloat)
                    -li2_approx(inv(z)) - one(BigFloat)/2*clog(-z)^2 - zeta2(BigFloat)
                else # |z|^2 <= 1
                    li2_approx(z)
                end
            else # Re(z) > 1/2
                if nz <= 2*rz
                    -li2_approx(one(BigFloat) - z) - clog(z)*clog1p(-z) + zeta2(BigFloat)
                else # |z|^2 > 2*Re(z)
                    -li2_approx(inv(z)) - one(BigFloat)/2*clog(-z)^2 - zeta2(BigFloat)
                end
            end
        end
    end
end
