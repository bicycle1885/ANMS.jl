using ANMS
using Base.Test

const tol = 1e-8

let
    # rosenbrock function
    a = 1.0
    b = 100.0
    f = x -> (a - x[1])^2 + b*(x[2] - x[1]^2)^2
    x₀ = [0., 0.]
    x, fmin = nelder_mead(f, x₀)
    @test norm(x .- [a, a]) < tol
    @test abs(fmin - 0.0) < tol
end

let
    # quadratic function
    f = x -> dot(x, x)
    for n in 2:10
        x₀ = ones(n)
        x, fmin = nelder_mead(f, x₀)
        @test norm(x .- zeros(n)) < tol
        @test abs(fmin - 0.0) < tol
    end
end

let
    run(`julia example.jl` |> DevNull)
end
