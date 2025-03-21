# Re[Li4(x)] for x in [-1,0]
function reli4_neg(x::Float64)::Float64
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

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    x*p/q
end

# Re[Li4(x)] for x in [0,1/2]
function reli4_half(x::Float64)::Float64
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

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6])

    x*p/q
end

# Re[Li4(x)] for x in [1/2,8/10]
function reli4_mid(x::Float64)::Float64
    cp = (
        3.2009826406098890447e-9, 9.9999994634837574160e-1,
       -2.9144851228299341318e+0, 3.1891031447462342009e+0,
       -1.6009125158511117090e+0, 3.5397747039432351193e-1,
       -2.5230024124741454735e-2
    )
    cq = (
        1.0000000000000000000e+0, -2.9769855248411488460e+0,
        3.3628208295110572579e+0, -1.7782471949702788393e+0,
        4.3364007973198649921e-1, -3.9535592340362510549e-2,
        5.7373431535336755591e-4
    )

    x2 = x*x
    x4 = x2*x2

    p = cp[1] + x * cp[2] + x2 * (cp[3] + x * cp[4]) +
        x4 * (cp[5] + x * cp[6] + x2 * cp[7])
    q = cq[1] + x * cq[2] + x2 * (cq[3] + x * cq[4]) +
        x4 * (cq[5] + x * cq[6] + x2 * cq[7])

    p/q
end

# Re[Li4(x)] for x in [8/10,1]
function reli4_one(x::Float64)::Float64
    l = log(x)
    l2 = l*l

    ZETA4_F64 +
    l*(ZETA3_F64 +
    l*(ZETA2_F64/2 +
    l*(11/36 - 1/6*log(-l) +
    l*(-1/48 +
    l*(-1/1440 +
    l2*(1/604800 - 1/91445760*l2))))))
end

"""
    reli4(x::Real)

Returns the real 4th order polylogarithm
``\\Re[\\operatorname{Li}_4(x)]`` of a real number ``x`` of type
`Real`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> reli4(1.0)
1.0823232337111381
```
"""
reli4(x::Real) = _reli4(float(x))

_reli4(x::Float16) = oftype(x, _reli4(Float32(x)))

_reli4(x::Float32) = oftype(x, _reli4(Float64(x)))

function _reli4(x::Float64)::Float64
    # transform to [-1,1]
    (x, rest, sgn) = if x < -1.0
        l = log(-x)^2
        (inv(x), -7/4*ZETA4_F64 + l*(-ZETA2_F64/2 - 1/24*l), -1.0)
    elseif x == -1.0
        return -7/8*ZETA4_F64
    elseif x < 1.0
        (x, 0.0, 1.0)
    elseif x == 1.0
        return ZETA4_F64
    else # x > 1.0
        l = log(x)^2
        (inv(x), 2*ZETA4_F64 + l*(ZETA2_F64 - 1/24*l), -1.0)
    end

    if x < 0.0
        rest + sgn*reli4_neg(x)
    elseif iszero(x)
        x
    elseif x < 0.5
        rest + sgn*reli4_half(x)
    elseif x < 0.8
        rest + sgn*reli4_mid(x)
    else # x <= 1.0
        rest + sgn*reli4_one(x)
    end
end

"""
    li4(z::Complex)

Returns the complex 4th order polylogarithm
``\\operatorname{Li}_4(z)`` of a complex number ``z`` of type
`Complex`.

If only real arguments ``z\\in\\mathbb{R}`` are considered and one is
interested only in the real part of the 4th order polylogarithm,
``\\Re[\\operatorname{Li}_4(z)]``, refer to the function
[`reli4`](@ref), which may be a faster alternative.

See also [`reli4`](@ref).

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li4(1.0 + 1.0im)
0.9593189135784194 + 1.138039196676983im
```
"""
li4(z::Complex) = _li4(float(z))

li4(z::Real) = li4(Complex(z))

_li4(z::ComplexF16) = oftype(z, _li4(ComplexF32(z)))

_li4(z::ComplexF32) = oftype(z, _li4(ComplexF64(z)))

function _li4(z::ComplexF64)::ComplexF64
    if iszero(imag(z))
        if real(z) <= 1.0
            Complex(reli4(real(z)), imag(z))
        else
            Complex(reli4(real(z)), -pi/6*log(real(z))^3)
        end
    else
        nz::Float64  = abs(z)
        pz::Float64  = angle(z)
        lnz::Float64 = log(nz)

        if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
            C = (
                -6.9444444444444444e-04, 1.6534391534391534e-06,
                -1.0935444136502338e-08, 1.0438378493934049e-10,
                -1.2165942300622435e-12, 1.6130006528350101e-14,
                -2.3428810452879340e-16
            )

            u  = lnz + pz*im
            u2 = u*u
            u4 = u2*u2
            u8 = u4*u4
            c1 = ZETA3_F64
            c2 = 0.82246703342411322
            c3 = 1/6*(11/6 - clog(-u))
            c4 = -1/48

            ZETA4_F64 + u2*(c2 + u2*c4) + u*(
                c1 +
                c3*u2 +
                u4*(C[1] + u2*C[2]) +
                u8*(C[3] + u2*C[4] + u4*(C[5] + u2*C[6])) +
                u8*u8*C[7]
            )
        else
            B = (
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

            (u, rest, sgn) = if nz <= 1.0
                (-clog1p(-z), 0.0 + 0.0im, 1.0)
            else # |z|^2 > 1
                arg = pz > 0.0 ? pz - pi : pz + pi
                lmz = lnz + arg*im # clog(z)
                lmz2 = lmz*lmz
                (-clog1p(-inv(z)), -7/4*ZETA4_F64 + lmz2*(-0.5*ZETA2_F64 - 1/24*lmz2), -1.0)
            end

            u2 = u*u
            u4 = u2*u2
            u8 = u4*u4

            rest + sgn*(
                u*B[1] +
                u2*(B[2] + u*B[3]) +
                u4*(B[4] + u*B[5] + u2*(B[6] + u*B[7])) +
                u8*(B[8] + u*B[9] + u2*(B[10] + u*B[11]) +
                u4*(B[12] + u*B[13] + u2*(B[14] + u*B[15]))) +
                u8*u8*(B[16] + u*B[17] + u2*B[18])
            )
        end
    end
end
