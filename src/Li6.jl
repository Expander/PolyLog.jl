"""
    li6(z::ComplexF64)::ComplexF64

Returns the complex 6th order polylogarithm
``\\operatorname{Li}_6(z)`` of a complex number ``z`` of type
`ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li6(1.0 + 1.0im)
0.996149796835317 + 1.0335544477237482im
```
"""
function li6(z::ComplexF64)::ComplexF64
    function clog(z)
        az::Float64 = angle(z)
        return 0.5*log(abs2(z)) + (imag(z) == 0.0 && az < 0.0 ? -az : az)*1.0im
    end

    bf = (
        1.0                   , -31.0/64.0             ,
        1.5241340877914952e-01, -3.4365555877057613e-02,
        5.7174797239368999e-03, -6.8180453746570645e-04,
        4.9960361948734493e-05, -4.9166051196039048e-07,
       -3.0632975161302164e-07,  1.4414599270849095e-08,
        3.7272438230924107e-09, -3.7300867345487607e-10,
       -5.1246526816085832e-11,  9.0541930956636683e-12,
        6.7381882615512517e-13, -2.1215831150303135e-13,
       -6.8408811719011698e-15,  4.8691178462005581e-15
    )

    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return zeta6 + 0.0im
        end
        if real(z) == -1.0
            return -31.0/32.0*zeta6 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        v::ComplexF64  = lnz + pz*im
        v2::ComplexF64 = v*v
        c1::Float64 = 1.0369277551433699 # zeta(5)
        c2::Float64 = 0.54116161685556910
        c3::Float64 = 0.20034281719326571
        c4::Float64 = 0.068538919452009435
        c5::ComplexF64 = (137.0/60.0 - clog(-v))/120.0
        c6::Float64 = -1.0/1440.0

        cs = (
            -1.6534391534391534e-05, 2.2964432686654909e-08,
            -9.9413128513657614e-11, 6.6912682653423394e-13,
            -5.7933058574392549e-15
        )

        return zeta6 + v * c1 +
            v2 * (c2 + v * c3 +
            v2 * (c4 + v * c5 +
            v2 * (c6 +
            v * (cs[1] +
            v2 * (cs[2] +
            v2 * (cs[3] +
            v2 * (cs[4] +
            v2 * (cs[5]))))))))
    end

    (u::ComplexF64, rest::ComplexF64, sgn::Float64) = if nz <= 1.0
        (-clog(1.0 - z), 0.0 + 0.0im, 1.0)
    else # nz > 1.0
        arg::Float64 = pz > 0.0 ? pz - pi : pz + pi
        lmz::ComplexF64 = lnz + arg*im # clog(z)
        lmz2::ComplexF64 = lmz*lmz
        (-clog(1.0 - 1.0/z), -31/16*zeta6 + lmz2*(-7/8*zeta4 + lmz2*(-1/24*zeta2 - 1/720*lmz2)), -1.0)
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
