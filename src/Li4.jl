"""
    li4(x::Float64)::Float64

Returns the real 4th order polylogarithm
``\\Re[\\operatorname{Li}_4(x)]`` of a real number ``x`` of type
`Float64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li4(1.0)
1.0823232337111381
```
"""
function li4(x::Float64)::Float64
    z2::Float64 = pi^2/6
    z3::Float64 = 1.2020569031595943
    z4::Float64 = pi^4/90

    # transform so that y in [-1,1]
    (y, r, s) = if x < -1.0
        l = log(-x)^2
        (inv(x), -7/4*z4 + l*(-z2/2 - 1/24*l), -1.0)
    elseif x == -1.0
        return -7/8*z4
    elseif x < 1.0
        (x, 0.0, 1.0)
    elseif x == 1.0
        return z4
    else # x > 1.0
        l = log(x)^2
        (inv(x), 2*z4 + l*(z2 - 1/24*l), -1.0)
    end

    if y < 0.0
        cp = (
            0.9999999999999999952e+0, -1.8532099956062184217e+0,
            1.1937642574034898249e+0, -3.1817912243893560382e-1,
            3.2268284189261624841e-2, -8.3773570305913850724e-4
        )
        cq = (
            1.0000000000000000000e+0, -1.9157099956062165688e+0,
            1.3011504531166486419e+0, -3.7975653506939627186e-1,
            4.5822723996558783670e-2, -1.8023912938765272341e-3,
            1.0199621542882314929e-5
        )
        y2 = y*y
        y4 = y2*y2
        p = cp[1] + y * cp[2] + y2 * (cp[3] + y * cp[4]) +
            y4 * (cp[5] + y * cp[6])
        q = cq[1] + y * cq[2] + y2 * (cq[3] + y * cq[4]) +
            y4 * (cq[5] + y * cq[6] + y2 * cq[7])
        r + s*y*p/q
    elseif y < 0.5
        cp = (
            1.0000000000000000414e+0, -2.0588072418045364525e+0,
            1.4713328756794826579e+0, -4.2608608613069811474e-1,
            4.2975084278851543150e-2, -6.8314031819918920802e-4
        )
        cq = (
            1.0000000000000000000e+0, -2.1213072418045207223e+0,
            1.5915688992789175941e+0, -5.0327641401677265813e-1,
            6.1467217495127095177e-2, -1.9061294280193280330e-3
        )
        y2 = y*y
        y4 = y2*y2
        p = cp[1] + y * cp[2] + y2 * (cp[3] + y * cp[4]) +
            y4 * (cp[5] + y * cp[6])
        q = cq[1] + y * cq[2] + y2 * (cq[3] + y * cq[4]) +
            y4 * (cq[5] + y * cq[6])
        r + s*y*p/q
    elseif y < 0.9
        cp = (
            3.1862424089379644171e-9, 9.9999994512820440983e-1,
           -3.7747775757314436877e+0, 5.6768601579106503485e+0,
           -4.2980010983965359383e+0, 1.6926918997627241534e+0,
           -3.1696542409634207715e-1, 2.0368233561036489459e-2
        )
        cq = (
            1.0000000000000000000e+0, -3.8372780020703093770e+0,
            5.9043463224806406145e+0, -4.6235611716349019621e+0,
            1.9221731606152274856e+0, -3.9773415425970717428e-1,
            3.2672161963833850649e-2, -4.5557334754282730111e-4
        )
        p = evalpoly(y, cp)
        q = evalpoly(y, cq)
        r + s*p/q
    else # y <= 1.0
        l = log(y)
        r + s*(z4 + l*(l*(pi^2/12 + l*(11/36 + l*(-1/48 + l*(-1/1440 + l^2/604800)) - log(abs(l))/6)) + z3))
    end
end

"""
    li4(z::ComplexF64)::ComplexF64

Returns the complex 4th order polylogarithm
``\\operatorname{Li}_4(z)`` of a complex number ``z`` of type
`ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li4(1.0 + 1.0im)
0.9593189135784193 + 1.1380391966769827im
```
"""
function li4(z::ComplexF64)::ComplexF64
    function clog(z)
        az::Float64 = angle(z)
        return 0.5*log(abs2(z)) + (imag(z) == 0.0 && az < 0.0 ? -az : az)*1.0im
    end

    z4::Float64 = 1.0823232337111382
    bf = (
        1.0                   , -7.0/16.0              ,
        1.1651234567901235e-01, -1.9820601851851852e-02,
        1.9279320987654321e-03, -3.1057098765432099e-05,
       -1.5624009114857835e-05,  8.4851235467732066e-07,
        2.2909616603189711e-07, -2.1832614218526917e-08,
       -3.8828248791720156e-09,  5.4462921032203321e-10,
        6.9608052106827254e-11, -1.3375737686445215e-11,
       -1.2784852685266572e-12,  3.2605628580248922e-13,
        2.3647571168618257e-14, -7.9231351220311617e-15
    )

    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return z4 + 0.0im
        end
        if real(z) == -1.0
            return -7.0/8.0*z4 + 0.0im
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
        c1::Float64 = 1.2020569031595943 # zeta(3)
        c2::Float64 = 0.82246703342411322
        c3::ComplexF64 = (11.0/6.0 - clog(-v))/6.0
        c4::Float64 = -1.0/48.0

        cs = (
            -6.9444444444444444e-04, 1.6534391534391534e-06,
            -1.0935444136502338e-08, 1.0438378493934049e-10,
            -1.2165942300622435e-12, 1.6130006528350101e-14,
            -2.3428810452879340e-16
        )

        return z4 + v2 * (c2 + v2 * c4) +
            v * (
                c1 +
                c3*v2 +
                v4*(cs[1] + v2*cs[2]) +
                v8*(cs[3] + v2*cs[4] + v4*(cs[5] + v2*cs[6])) +
                v8*v8*cs[7]
            )
    end

    (u::ComplexF64, rest::ComplexF64, sgn::Float64) = if nz <= 1.0
        (-clog(1.0 - z), 0.0 + 0.0im, 1.0)
    else # nz > 1.0
        arg::Float64 = pz > 0.0 ? pz - pi : pz + pi
        lmz::ComplexF64 = lnz + arg*im # clog(z)
        lmz2::ComplexF64 = lmz*lmz
        (-clog(1.0 - 1.0/z), 1.0/360.0*(-7.0*pi^4 + lmz2*(-30.0*pi^2 - 15.0*lmz2)), -1.0)
    end

    u2::ComplexF64 = u*u
    u4::ComplexF64 = u2*u2
    u8::ComplexF64 = u4*u4

    rest + sgn * (
       u*bf[1] +
       u2*(bf[2] + u*bf[3]) +
       u4*(bf[4] + u*bf[5] + u2*(bf[6] + u*bf[7])) +
       u8*(bf[8] + u*bf[9] + u2*(bf[10] + u*bf[11]) +
           u4*(bf[12] + u*bf[13] + u2*(bf[14] + u*bf[15]))) +
       u8*u8*(bf[16] + u*bf[17] + u2*bf[18])
    )
end
