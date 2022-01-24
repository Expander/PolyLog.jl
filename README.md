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

li0(1.0)         # real polylogarithm Re[Li_0(x)] of order 0
li1(1.0)         # real polylogarithm Re[Li_1(x)] of order 1
li2(1.0)         # real polylogarithm Re[Li_2(x)] of order 2 (dilogarithm)
li3(1.0)         # real polylogarithm Re[Li_3(x)] of order 3 (trilogarithm)
li4(1.0)         # real polylogarithm Re[Li_4(x)] of order 4
li(10, 1.0)      # real polylogarithm Re[Li_n(x)] of integer order n >= 0 (here: n = 10)

li2(1.0 + 1.0im) # complex polylogarithm Li_2(z) of order 2 (dilogarithm)
li3(1.0 + 1.0im) # complex polylogarithm Li_3(z) of order 3 (trilogarithm)
li4(1.0 + 1.0im) # complex polylogarithm Li_4(z) of order 4
li5(1.0 + 1.0im) # complex polylogarithm Li_5(z) of order 5
li6(1.0 + 1.0im) # complex polylogarithm Li_6(z) of order 6
```


Documentation
-------------

[https://docs.juliahub.com/PolyLog/](https://docs.juliahub.com/PolyLog/)


Notes
-----

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Julia.


Citation
--------

~~~.bibtex
@software{PolyLog.jl,
    author       = {{Alexander Voigt}},
    title        = {{PolyLog.jl}},
    year         = {2021},
    version      = {1.7.0},
    url          = {https://github.com/Expander/PolyLog.jl},
    note         = {[License: MIT]}
}
~~~


Copying
-------

PolyLog.jl is licenced under the MIT License.


Links
-----

Refer to the package
[Polylogarithms.jl](https://github.com/mroughan/Polylogarithms.jl) for
a Julia implementation of polylogarithms of arbitrary order.
