PolyLog.jl
==========

[![test](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml/badge.svg)](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml)

The PolyLog.jl package provides Julia implementations of real and
complex polylogarithms, including the real and complex dilogarithm and
trilogarithm.


Example
-------

```.jl
include("src/PolyLog.jl")

println(li2(1.0))
println(li2(1.0 + 1.0im))
println(li3(1.0 + 1.0im))
println(li4(1.0 + 1.0im))
```


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
    year         = "2021",
    version      = {1.2.0},
    url          = {https://github.com/Expander/PolyLog.jl},
    note         = "[License: MIT]"
}
~~~


Copying
-------

PolyLog.jl is licenced under the MIT License.


Links
-----

Refer to the package
[Polylogarithms.jl](https://github.com/mroughan/Polylogarithms.jl) for
a Julia implementation of polylogarithm of arbitrary order.
