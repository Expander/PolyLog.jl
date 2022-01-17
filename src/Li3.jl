# Li3(x) for x in [-1,0]
function li3_neg(x::Float64)::Float64
    cp = (
        0.9999999999999999795, -2.0281801754117129576,
        1.4364029887561718540, -0.42240680435713030268,
        0.047296746450884096877, -0.0013453536579918419568
    )
    cq = (
        1.0000000000000000000, -2.1531801754117049035,
        1.6685134736461140517, -0.56684857464584544310,
        0.081999463370623961084, -0.0040756048502924149389,
        0.000034316398489103212699
    )

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    return x*p/q
end

# Li3(x) for x in [0,1/2]
function li3_pos(x::Float64)::Float64
    cp = (
        0.9999999999999999893, -2.5224717303769789628,
        2.3204919140887894133, -0.93980973288965037869,
        0.15728950200990509052, -0.0075485193983677071129
    )
    cq = (
        1.0000000000000000000, -2.6474717303769836244,
        2.6143888433492184741, -1.1841788297857667038,
        0.24184938524793651120, -0.018220900115898156346,
        0.00024927971540017376759
    )

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    return x*p/q
end

"""
    li3(z::Float64)::Float64

Returns the real trilogarithm of a real number `z` of type `Float64`.

Author: Alexander Voigt

License: MIT

# Example
```julia
li3(1.0)
```
"""
function li3(x::Float64)::Float64
    z2::Float64 = 1.6449340668482264
    z3::Float64 = 1.2020569031595943

    # transformation so that y in [-1,0] and z in [0,1/2]
    (neg, pos, s, r) = if x < -1.0
        l = log(-x)
        (li3_neg(inv(x)), 0.0, 1.0, -l*(z2 + 1/6*l*l))
    elseif x == -1.0
        return -0.75*z3 # todo: necessary?
    elseif x < 0.0
        (li3_neg(x), 0.0, 1.0, 0.0)
    elseif x == 0.0
        return 0.0
    elseif x < 0.5
        (0.0, li3_pos(x), 1.0, 0.0)
    elseif x == 0.5
        return 0.53721319360804020 # todo necessary?
    elseif x < 1.0
        l = log(x)
        (li3_neg((x - 1.0)/x), li3_pos(1.0 - x), -1.0, z3 + l*(z2 + l*(-0.5*log(1.0 - x) + 1/6*l)))
    elseif x == 1.0
        return z3
    elseif x < 2.0
        l = log(x)
        (li3_neg(1.0 - x), li3_pos((x - 1.0)/x), -1.0, z3 + l*(z2 + l*(-0.5*log(x - 1.0) + 1/6*l)))
    else # x >= 2.0
        l = log(x)
        (0.0, li3_pos(inv(x)), 1.0, l*(2*z2 - 1/6*l*l))
    end

    r + s*(neg + pos)
end

"""
    li3(z::ComplexF64)::ComplexF64

Returns the complex trilogarithm of a complex number `z` of type `ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```julia
li3(1.0 + 1.0im)
```
"""
function li3(z::ComplexF64)::ComplexF64
    function clog(z)
        az::Float64 = angle(z)
        return 0.5*log(abs2(z)) + (imag(z) == 0.0 && az < 0.0 ? -az : az)*1.0im
    end

    z2::Float64 = 1.6449340668482264
    z3::Float64 = 1.2020569031595943
    bf = (
        1.0, -3.0/8.0, 17.0/216.0, -5.0/576.0,
        1.2962962962962963e-04,  8.1018518518518519e-05,
       -3.4193571608537595e-06, -1.3286564625850340e-06,
        8.6608717561098513e-08,  2.5260875955320400e-08,
       -2.1446944683640648e-09, -5.1401106220129789e-10,
        5.2495821146008294e-11,  1.0887754406636318e-11,
       -1.2779396094493695e-12, -2.3698241773087452e-13,
        3.1043578879654623e-14,  5.2617586299125061e-15
    )

    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return z3 + 0.0im
        end
        if real(z) == -1.0
            return -0.75*z3 + 0.0im
        end
        if real(z) == 0.5
            return 0.53721319360804020 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        v::ComplexF64  = lnz + pz*im
        v2::ComplexF64 = v*v
        v4::ComplexF64 = v2*v2
        v8::ComplexF64 = v4*v4
        c0::ComplexF64 = z3 + v*(z2 - v2/12.0)
        c1::ComplexF64 = 0.25 * (3.0 - 2.0*clog(-v))

        cs = (
            -3.4722222222222222e-03, 1.1574074074074074e-05,
            -9.8418997228521038e-08, 1.1482216343327454e-09,
            -1.5815724990809166e-11, 2.4195009792525152e-13,
            -3.9828977769894877e-15
        )

        return c0 +
            c1*v2 +
            v4*(cs[1] + v2*cs[2]) +
            v8*(cs[3] + v2*cs[4] + v4*(cs[5] + v2*cs[6])) +
            v8*v8*cs[7]
    end

    (u::ComplexF64, rest::ComplexF64) = if nz <= 1.0
        (-clog(1.0 - z), 0.0 + 0.0im)
    else # nz > 1.0
        arg::Float64 = pz > 0.0 ? pz - pi : pz + pi
        lmz::ComplexF64 = lnz + arg*im # clog(z)
        (-clog(1.0 - 1.0/z), -lmz*(lmz*lmz/6.0 + z2))
    end

    u2::ComplexF64 = u*u
    u4::ComplexF64 = u2*u2
    u8::ComplexF64 = u4*u4

    rest +
    u*bf[1] +
    u2*(bf[2] + u*bf[3]) +
    u4*(bf[4] + u*bf[5] + u2*(bf[6] + u*bf[7])) +
    u8*(bf[8] + u*bf[9] + u2*(bf[10] + u*bf[11]) +
        u4*(bf[12] + u*bf[13] + u2*(bf[14] + u*bf[15]))) +
    u8*u8*(bf[16] + u*bf[17] + u2*bf[18])
end
