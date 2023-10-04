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


# rational function approximation of Re[Li2(x)] for x in [0, 1/2]
function reli2_approx(x::Float32)::Float32
    cp = (1.00000020f0, -0.780790946f0, 0.0648256871f0)
    cq = (1.00000000f0, -1.03077545f0, 0.211216710f0)

    x2 = x*x

    p = cp[1] + x * cp[2] + x2 * cp[3]
    q = cq[1] + x * cq[2] + x2 * cq[3]

    x*p/q
end


# approximation of Re[Li2(x)] for x in [0, 1/2]
# todo(alex): benchmark this routine
function reli2_approx(x::BigFloat)::BigFloat
    u = -log1p(-x)
    u2 = u*u

    p = u2 # powers of u2
    sum = one(x)/36

    for n in 2:typemax(Int64)
        old_sum = sum
        sum += p*bernoulli(2*n)/fac(2*n + 1)
        sum == old_sum && break
        p *= u2
    end

    u + u2*(-one(x)/4 + u*sum)
end


# Bernoulli number
function bernoulli(n)
    A = Vector{Rational{BigInt}}(undef, n + 1)
    for m in 0:n
        A[m + 1] = 1//(m + 1)
        for j in m:-1:1
            A[j] = j*(A[j] - A[j + 1])
        end
    end
    A[1]
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

function _reli2(x::Float32)::Float32
    # transform to [0, 1/2]
    if x < -1.0f0
        l = log(1.0f0 - x)
        reli2_approx(1.0f0/(1.0f0 - x)) - zeta2F32 + l*(0.5f0*l - log(-x))
    elseif x == -1.0f0
        -0.5f0*zeta2F32
    elseif x < 0.0f0
        -reli2_approx(x/(x - 1.0f0)) - 0.5f0*log1p(-x)^2
    elseif x == 0.0f0
        0.0f0
    elseif x < 0.5f0
        reli2_approx(x)
    elseif x < 1.0f0
        -reli2_approx(1.0f0 - x) + zeta2F32 - log(x)*log1p(-x)
    elseif x == 1.0f0
        zeta2F32
    elseif x < 2.0f0
        l = log(x)
        reli2_approx(1.0f0 - 1.0f0/x) + zeta2F32 - l*(log(1.0f0 - 1.0f0/x) + 0.5f0*l)
    else
        -reli2_approx(1.0f0/x) + 2.0f0*zeta2F32 - 0.5f0*log(x)^2
    end
end

function _reli2(x::Float64)::Float64
    # transform to [0, 1/2]
    if x < -1.0
        l = log(1.0 - x)
        reli2_approx(1.0/(1.0 - x)) - zeta2 + l*(0.5*l - log(-x))
    elseif x == -1.0
        -0.5*zeta2
    elseif x < 0.0
        -reli2_approx(x/(x - 1.0)) - 0.5*log1p(-x)^2
    elseif x == 0.0
        0.0
    elseif x < 0.5
        reli2_approx(x)
    elseif x < 1.0
        -reli2_approx(1.0 - x) + zeta2 - log(x)*log1p(-x)
    elseif x == 1.0
        zeta2
    elseif x < 2.0
        l = log(x)
        reli2_approx(1.0 - 1.0/x) + zeta2 - l*(log(1.0 - 1.0/x) + 0.5*l)
    else
        -reli2_approx(1.0/x) + 2.0*zeta2 - 0.5*log(x)^2
    end
end

function _reli2(x::BigFloat)::BigFloat
    # transform to [0, 1/2]
    if x < -one(x)
        l = log(one(x) - x)
        reli2_approx(inv(one(x) - x)) - zeta(2, typeof(x)) + l*(one(x)/2*l - log(-x))
    elseif x == -one(x)
        -one(x)/2*zeta(2, typeof(x))
    elseif x < zero(x)
        -reli2_approx(x/(x - one(x))) - one(x)/2*log1p(-x)^2
    elseif iszero(x)
        zero(x)
    elseif x < one(x)/2
        reli2_approx(x)
    elseif x == one(x)/2
        BigFloat(pi)^2/12 - log(big(2))^2/2
    elseif x < one(x)
        -reli2_approx(one(x) - x) + zeta(2, typeof(x)) - log(x)*log1p(-x)
    elseif x == one(x)
        zeta(2, typeof(x))
    elseif x < 2*one(x)
        l = log(x)
        reli2_approx(one(x) - inv(x)) + zeta(2, typeof(x)) - l*(log(one(x) - inv(x)) + one(x)/2*l)
    else
        -reli2_approx(inv(x)) + 2*zeta(2, typeof(x)) - one(x)/2*log(x)^2
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
0.616850275068085 + 1.4603621167531196im
```
"""
li2(z::Complex) = _li2(float(z))

li2(z::Real) = li2(Complex(z))

_li2(z::ComplexF16) = oftype(z, _li2(ComplexF32(z)))

function _li2(z::ComplexF32)::ComplexF32
    clog(z) = log(abs(z)) + angle(z)*1.0f0im

    rz, iz = reim(z)

    if iz == 0.0f0
        if rz <= 1.0f0
            reli2(rz)
        else # Re(z) > 1
            reli2(rz) - pi*clog(rz)*1.0f0im
        end
    else
        nz = abs2(z)

        if nz < eps(Float32)
            z*(1.0f0 + 0.25f0*z)
        else
            if rz <= 0.5f0
                if nz > 1.0f0
                    -li2_approx(-clog(1.0f0 - inv(z))) - 0.5f0*clog(-z)^2 - zeta2F32
                else # |z|^2 <= 1
                    li2_approx(-clog(1.0f0 - z))
                end
            else # Re(z) > 1/2
                if nz <= 2.0f0*rz
                    l = -clog(z)
                    -li2_approx(l) + l*clog(1.0f0 - z) + zeta2F32
                else # |z|^2 > 2*Re(z)
                    -li2_approx(-clog(1.0f0 - inv(z))) - 0.5f0*clog(-z)^2 - zeta2F32
                end
            end
        end
    end
end

function _li2(z::ComplexF64)::ComplexF64
    clog(z) = 0.5*log(abs2(z)) + angle(z)*1.0im

    rz, iz = reim(z)

    if iz == 0.0
        if rz <= 1.0
            reli2(rz)
        else # Re(z) > 1
            reli2(rz) - pi*log(rz)*1.0im
        end
    else
        nz = abs2(z)

        if nz < eps(Float64)
            z*(1.0 + 0.25*z)
        else
            if rz <= 0.5
                if nz > 1.0
                    -li2_approx(-clog(1.0 - inv(z))) - 0.5*clog(-z)^2 - zeta2
                else # |z|^2 <= 1
                    li2_approx(-clog(1.0 - z))
                end
            else # Re(z) > 1/2
                if nz <= 2.0*rz
                    l = -clog(z)
                    -li2_approx(l) + l*clog(1.0 - z) + zeta2
                else # |z|^2 > 2*Re(z)
                    -li2_approx(-clog(1.0 - inv(z))) - 0.5*clog(-z)^2 - zeta2
                end
            end
        end
    end
end
