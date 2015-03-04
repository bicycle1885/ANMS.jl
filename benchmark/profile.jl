using ANMS

quadratic(x) = dot(x, x)

let
    n = 200
    x₀ = ones(n)
    nelder_mead(quadratic, x₀)
    @profile nelder_mead(quadratic, x₀)
    Profile.print()
end
