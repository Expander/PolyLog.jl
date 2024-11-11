"""
    li6(z::Complex)

Returns the complex 6th order polylogarithm
``\\operatorname{Li}_6(z)`` of a complex number ``z`` of type
`Complex`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li6(1.0 + 1.0im)
0.996149796835317 + 1.0335544477237482im
```
"""
li6(z::Complex) = _li6(float(z))

li6(z::Real) = li6(Complex(z))

_li6(z::ComplexF16) = oftype(z, _li6(ComplexF32(z)))

_li6(z::ComplexF32) = oftype(z, _li6(ComplexF64(z)))

function _li6(z::ComplexF64)::ComplexF64
    if iszero(imag(z))
        if iszero(real(z))
            return z
        end
        if real(z) == 1.0
            return complex(ZETA6_F64)
        end
        if real(z) == -1.0
            return complex(-31/32*ZETA6_F64)
        end
    end

    nz::Float64  = abs(z)
    pz::Float64  = angle(z)
    lnz::Float64 = log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        C = (
            -1.6534391534391534e-05, 2.2964432686654909e-08,
            -9.9413128513657614e-11, 6.6912682653423394e-13,
            -5.7933058574392549e-15
        )

        u  = lnz + pz*im
        u2 = u*u
        u4 = u2*u2
        c1 = ZETA5_F64
        c2 = 0.54116161685556910
        c3 = 0.20034281719326571
        c4 = 0.068538919452009435
        c5 = 1/120*(137/60 - clog(-u))
        c6 = -1/1440

        ZETA6_F64 + u2*c2 + u4*(c4 + u2*c6) +
        u*(c1 + u2*c3 + u4*(c5 + u2*C[1]) +
           u4*u4*(C[2] + u2*C[3] + u4*(C[4] + u2*C[5])))
    else
        B = (
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

        (u, rest, sgn) = if nz <= 1.0
            (-clog1p(-z), 0.0 + 0.0im, 1.0)
        else # nz > 1.0
            arg = pz > 0.0 ? pz - pi : pz + pi
            lmz = lnz + arg*im # clog(z)
            lmz2 = lmz*lmz
            (-clog1p(-inv(z)), -31/16*ZETA6_F64 + lmz2*(-7/8*ZETA4_F64 + lmz2*(-1/24*ZETA2_F64 - 1/720*lmz2)), -1.0)
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

li6(::Missing) = missing
