using ConstrainedLeastSquares
using Random
using Test

const CLSQ = ConstrainedLeastSquares
const Z = [
    (:a, [1.0, 1.9], [0.1, 0.1]),
    (:b, 1.0, 0.1),
    (:c, [2.9, 5.0], [0.1, 0.05])
]

@testset "creation of parameter index" begin
    target = Dict([:a => 1:2, :b => 3, :c => 4:5, :a_1 => 1, :a_2 => 2, :c_1 => 4, :c_2 => 5])
    @test CLSQ._build_index(Z) == target
    β = [
        (:a, [1.0, 1.9]),
        (:b, 1.0),
        (:c, [2.9, 5.0])
    ]
    @test CLSQ._build_index(β) == target
end

@testset "creation of parameter array" begin
    target = [1.0, 1.9, 1.0, 2.9, 5.0] # values defined from Z
    @test CLSQ._build_var_array(Z, 2) == target
    target = [0.1, 0.1, 0.1, 0.1, 0.05]  # uncertainties defined from Z
    @test CLSQ._build_var_array(Z, 3) == target
end

@testset "set covariance" begin
    z = [
        (:a, 1.5, 0.1),
        (:b, [2., 3.], [0.1, 0.14]),
        (:c, 4.0, 0.2)
    ]
    β = [(:A, 1.5)] # not needed here
    inp = LSQInput(β, z)
    
    @test_throws ArgumentError covariance!(inp, :a, :b, 10.0)
    covariance!(inp, :a, :c, 10.0)
    @test inp.Σ[1, 4] == inp.Σ[4, 1] == 10.0
    covariance!(inp, :b_1, :b_2, 20.0)
    @test inp.Σ[2, 3] == inp.Σ[3, 2] == 20.0
end

@testset "regression test with linear model" begin
    Random.seed!(42)
    # toy data
    N = 4
    xs = collect(Float64, 1:N) .+ randn(N)/10
    ys = xs * 2.2 .+ 3.1 .+ randn(N)/10

    # measured parameters
    z = [
        (:x, xs, 0.05xs),
        (:y, ys, 0.05ys)
    ]

    # unknown parameters
    β₀ = [
        (:a,  1.0),
        (:b,  1.0)
    ]

    linear_input = LSQInput(β₀, z)
    linear_model = @withinput linear_input begin
        return @. a*x + b - y
    end

    output = solve(linear_model, linear_input)

    @test output.χ²   ≈ 0.20962691792604748
    @test output.β[1] ≈ 2.113191779601563
    @test output.β[2] ≈ 3.3110566163534476
end