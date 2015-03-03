using ANMS

quadratic(x) = dot(x, x)

let
    x₀ = ones(2)
    nelder_mead(quadratic, x₀)

    for n in [2, 5, 10, 50, 100, 200]
        x₀ = ones(n)
        @time x, fmin = nelder_mead(quadratic, x₀)
        @assert fmin < 1.0e-5
    end
end
