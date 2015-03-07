using ANMS

a = 1.0
b = 100.0
# rosenbrock function takes its minimum value 0 at (a, a^2)
f = x -> (a - x[1])^2 + b*(x[2] - x[1]^2)^2
x₀ = [0., 0.]
x, fmin = nelder_mead(f, x₀)
@show x fmin

@assert norm(x - [a, a]) < 1.0e-8
@assert fmin < 1.0e-8
