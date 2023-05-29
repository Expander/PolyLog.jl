PolyLog.jl
==========

[![test](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml/badge.svg)](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml)

The PolyLog.jl package provides Julia implementations of real and
complex polylogarithms, including the real and complex dilogarithm and
trilogarithm.


Example
-------

```.jl
using PolyLog

# real polylogarithms for real arguments
reli1(1.0)          # Re[Li_1(x)]
reli2(1.0)          # Re[Li_2(x)] (dilogarithm)
reli3(1.0)          # Re[Li_3(x)] (trilogarithm)
reli4(1.0)          # Re[Li_4(x)]
reli(10, 1.0)       # Re[Li_n(x)] for all integers n (here: n = 10)
reli(1//2, 1.0)     # Re[Li_r(x)] for all half-integer rationals r=n/2 (here: n = 1)

# complex polylogarithms for real or complex arguments and integer order
li0(1.0 + 1.0im)      # Li_0(z)
li1(1.0 + 1.0im)      # Li_1(z)
li2(1.0 + 1.0im)      # Li_2(z) (dilogarithm)
li3(1.0 + 1.0im)      # Li_3(z) (trilogarithm)
li4(1.0 + 1.0im)      # Li_4(z)
li5(1.0 + 1.0im)      # Li_5(z)
li6(1.0 + 1.0im)      # Li_6(z)
li(10, 1.0 + 1.0im)   # Li_n(z) for all integers n (here: n = 10)
li(1//2, 1.0 + 1.0im) # Li_r(x) for all half-integer rationals r=n/2 (here: n = 1)
```


Documentation
-------------

[https://docs.juliahub.com/PolyLog/](https://docs.juliahub.com/PolyLog/)


Notes
-----

The implementation of the real dilogarithm is an adaptation of
[[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Julia.

The implementation of the general n-th order polylogarithm is an
adaptation of [[arXiv:2010.09860](https://arxiv.org/abs/2010.09860)].


Copying
-------

PolyLog.jl is licenced under the MIT License.


Links
-----

Refer to the package
[Polylogarithms.jl](https://github.com/mroughan/Polylogarithms.jl) for
a Julia implementation of polylogarithms of arbitrary complex order.
