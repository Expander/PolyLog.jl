"""
    li5(z::Complex)

Returns the complex 5th order polylogarithm
``\\operatorname{Li}_5(z)`` of a complex number ``z`` of type
`ComplexF`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li5(1.0 + 1.0im)
0.9874666591701124 + 1.068441607107422im
```
"""
li5(z::Complex) = _li5(float(z))

_li5(z::ComplexF16) = oftype(z, _li5(ComplexF32(z)))

_li5(z::ComplexF32) = oftype(z, _li5(ComplexF64(z)))

function _li5(z::ComplexF64)::ComplexF64
    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return ZETA5_F64 + 0.0im
        end
        if real(z) == -1.0
            return -15/16*ZETA5_F64 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        C = (
            -1.1574074074074074e-04, 2.0667989417989418e-07,
            -1.0935444136502338e-09, 8.6986487449450412e-12,
            -8.6899587861588824e-14, 1.0081254080218813e-15
        )

        u  = lnz + pz*im
        u2 = u*u
        u4 = u2*u2
        u8 = u4*u4
        c1 = ZETA4_F64
        c2 = 0.5*ZETA3_F64
        c3 = 0.27415567780803774
        c4 = 1/24*(25/12 - clog(-u))
        c5 = -1/240

        ZETA5_F64 + u*(c1 + u2*(c3 + u2*c5)) +
        u2*c2 + u4*(c4 + u2*C[1]) +
        u8*(C[2] + u2*C[3] + u4*(C[4] + u2*C[5])) +
        u8*u8*C[6]
    else
        B = (
            1.0                   , -15.0/32.0             ,
            1.3953189300411523e-01, -2.8633777006172840e-02,
            4.0317412551440329e-03, -3.3985018004115226e-04,
            4.5445184621617666e-06,  2.3916808048569012e-06,
           -1.2762692600122747e-07, -3.1628984306505932e-08,
            3.2848118445335192e-09,  4.7613713995660579e-10,
           -8.0846898171909830e-11, -7.2387648587737207e-12,
            1.9439760115173968e-12,  1.0256978405977236e-13,
           -4.6180551009884830e-14, -1.1535857196470580e-15,
            1.0903545401333394e-15
        )

        (u, rest) = if nz <= 1.0
            (-clog(1.0 - z), 0.0 + 0.0im)
        else # nz > 1.0
            arg = pz > 0.0 ? pz - pi : pz + pi
            lmz = lnz + arg*im # clog(z)
            lmz2 = lmz*lmz
            (-clog(1.0 - inv(z)), lmz*(-7/4*ZETA4_F64 - lmz2*(1/6*ZETA2_F64 + 1/120*lmz2)))
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
        u8*u8*(B[16] + u*B[17] + u2*(B[18] + u*B[19]))
    end
end
