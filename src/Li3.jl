# Li3(x) for x in [-1,0]
function li3_neg(x::Float64)::Float64
    cp = (
        0.9999999999999999795e+0, -2.0281801754117129576e+0,
        1.4364029887561718540e+0, -4.2240680435713030268e-1,
        4.7296746450884096877e-2, -1.3453536579918419568e-3
    )
    cq = (
        1.0000000000000000000e+0, -2.1531801754117049035e+0,
        1.6685134736461140517e+0, -5.6684857464584544310e-1,
        8.1999463370623961084e-2, -4.0756048502924149389e-3,
        3.4316398489103212699e-5
    )

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    x*p/q
end

# Li3(x) for x in [0,1/2]
function li3_pos(x::Float64)::Float64
    cp = (
        0.9999999999999999893e+0, -2.5224717303769789628e+0,
        2.3204919140887894133e+0, -9.3980973288965037869e-1,
        1.5728950200990509052e-1, -7.5485193983677071129e-3
    )
    cq = (
        1.0000000000000000000e+0, -2.6474717303769836244e+0,
        2.6143888433492184741e+0, -1.1841788297857667038e+0,
        2.4184938524793651120e-1, -1.8220900115898156346e-2,
        2.4927971540017376759e-4
    )

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    x*p/q
end

"""
    li3(x::Real)

Returns the real trilogarithm ``\\Re[\\operatorname{Li}_3(x)]`` of a
real number ``x`` of type `Real`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li3(1.0)
1.2020569031595942
```
"""
li3(x::Real) = _li3(float(x))

_li3(x::Float16) = oftype(x, _li3(Float32(x)))

_li3(x::Float32) = oftype(x, _li3(Float64(x)))

function _li3(x::Float64)::Float64
    # transformation to [-1,0] and [0,1/2]
    (neg, pos, sgn, rest) = if x < -1.0
        l = log(-x)
        (li3_neg(inv(x)), 0.0, 1.0, -l*(zeta2 + 1/6*l*l))
    elseif x == -1.0
        return -0.75*zeta3
    elseif x < 0.0
        (li3_neg(x), 0.0, 1.0, 0.0)
    elseif x == 0.0
        return 0.0
    elseif x < 0.5
        (0.0, li3_pos(x), 1.0, 0.0)
    elseif x == 0.5
        return 0.53721319360804020
    elseif x < 1.0
        l = log(x)
        (li3_neg((x - 1.0)/x), li3_pos(1.0 - x), -1.0, zeta3 + l*(zeta2 + l*(-0.5*log(1.0 - x) + 1/6*l)))
    elseif x == 1.0
        return zeta3
    elseif x < 2.0
        l = log(x)
        (li3_neg(1.0 - x), li3_pos((x - 1.0)/x), -1.0, zeta3 + l*(zeta2 + l*(-0.5*log(x - 1.0) + 1/6*l)))
    else # x >= 2.0
        l = log(x)
        (0.0, li3_pos(inv(x)), 1.0, l*(2*zeta2 - 1/6*l*l))
    end

    rest + sgn*(neg + pos)
end

"""
    li3(z::Complex)

Returns the complex trilogarithm ``\\operatorname{Li}_3(z)`` of a
complex number ``z`` of type `Complex`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li3(1.0 + 1.0im)
0.8711588834109378 + 1.267083441888924im
```
"""
li3(z::Complex) = _li3(float(z))

_li3(z::ComplexF16) = oftype(z, _li3(ComplexF32(z)))

_li3(z::ComplexF32) = oftype(z, _li3(ComplexF64(z)))

function _li3(z::ComplexF64)::ComplexF64
    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return zeta3 + 0.0im
        end
        if real(z) == -1.0
            return -0.75*zeta3 + 0.0im
        end
        if real(z) == 0.5
            return 0.53721319360804020 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        C = (
            -3.4722222222222222e-03, 1.1574074074074074e-05,
            -9.8418997228521038e-08, 1.1482216343327454e-09,
            -1.5815724990809166e-11, 2.4195009792525152e-13,
            -3.9828977769894877e-15
        )

        u = lnz + pz*im
        u2 = u*u
        u4 = u2*u2
        u8 = u4*u4
        c0 = zeta3 + u*(zeta2 - 1/12*u2)
        c1 = 0.25*(3.0 - 2.0*clog(-u))

        c0 +
        c1*u2 +
        u4*(C[1] + u2*C[2]) +
        u8*(C[3] + u2*C[4] + u4*(C[5] + u2*C[6])) +
        u8*u8*C[7]
    else
        B = (
            1.0, -3.0/8.0, 17.0/216.0, -5.0/576.0,
            1.2962962962962963e-04,  8.1018518518518519e-05,
           -3.4193571608537595e-06, -1.3286564625850340e-06,
            8.6608717561098513e-08,  2.5260875955320400e-08,
           -2.1446944683640648e-09, -5.1401106220129789e-10,
            5.2495821146008294e-11,  1.0887754406636318e-11,
           -1.2779396094493695e-12, -2.3698241773087452e-13,
            3.1043578879654623e-14,  5.2617586299125061e-15
        )

        (u, rest) = if nz <= 1.0
            (-clog(1.0 - z), 0.0 + 0.0im)
        else # nz > 1.0
            arg = pz > 0.0 ? pz - pi : pz + pi
            lmz = lnz + arg*im # clog(z)
            (-clog(1.0 - inv(z)), -lmz*(1/6*lmz*lmz + zeta2))
        end

        u2 = u*u
        u4 = u2*u2
        u8 = u4*u4

        rest +
        u*B[1] +
        u2*(B[2] + u*B[3]) +
        u4*(B[4] + u*B[5] + u2*(B[6] + u*B[7])) +
        u8*(B[8] + u*B[9] + u2*(B[10] + u*B[11]) +
            u4*(B[12] + u*B[13] + u2*(B[14] + u*B[15]))) +
        u8*u8*(B[16] + u*B[17] + u2*B[18])
    end
end
