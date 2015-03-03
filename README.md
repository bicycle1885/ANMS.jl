# ANMS - Adaptive Nelder-Mead Simplex Optimization

[![Build Status](https://travis-ci.org/bicycle1885/ANMS.jl.svg?branch=master)](https://travis-ci.org/bicycle1885/ANMS.jl)

Nelder-Mead algorithm with adaptive parameters based on Fuchang Gao and Lixing Han (2010), Springer US. "Implementing the Nelder-Mead simplex algorithm with adaptive parameters" ([doi:10.1007/s10589-010-9329-3](http://link.springer.com/article/10.1007/s10589-010-9329-3)).

This algorithm would be faster than the standard Nelder-Mead algorithm when optimizing a high dimensional objective function.

## API

`ANMS` exports a function:

```julia
nelder_mead(f::Function, x₀::Vector{Float64}; iteration::Int=1_000_000, ftol::Float64=1.0e-8, xtol::Float64=1.0e-8)
```

where

* `f::Function`: objective function
* `x₀::Vector{Float64}`: starting point (one vertex of a simplex).

Finally, this function returns a tuple of `(minimizer::Vector{Float64}, function-value::Float64)`.


## Example

```julia
using ANMS

a = 1.0
b = 100.0
# rosenbrock function takes its minimum value 0 at (a, a^2)
f = x -> (a - x[1])^2 + b*(x[2] - x[1]^2)^2
x₀ = [0., 0.]
x, fmin = nelder_mead(f, x₀)
@show x fmin
```

Output:

    x => [0.9999999976422846,0.9999999951910932]
    fmin => 6.432600122573516e-18
