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
    @test uncertainty_β(output)[1] ≈ 0.21165344125510738
    @test uncertainty_β(output)[2] ≈ 0.419366917354832
end

@testset "regression test with example from Lars Nielsen 2001" begin
    # see Lars Nielsen, 2001: Evaluation of Measurements by the Method of Least Squares
    # measured parameters
    z = [
        (:ms, 199.988816, 0.000010), # known sum of weights in g
        (:mr, 199.999043, 0.000008), # reference weight in g
        (:ρr, 7833.01, 0.29), # densities in kg/m³
        (:ρ,  7965.76, 0.71), 
        (:a, 1.1950, 0.0035), # air density
        (:I, [199.988617 # scale indication in "div"
            199.988608
            174.992133
            175.009992
            150.013558
            149.980675
            125.002083
            124.984217
            100.005650
            99.982925
            74.986433
            75.004325
            50.007892
            49.974992
            24.978533
            24.996417],
        ones(16)*2.3e-5),
        (:Ir, [199.998867, 199.998875], ones(2)*2.3e-5)
    ]

    # unknown parameters
    β₀ = [
        (:f,  1.0), # g per div
        (:A,  1.0), # non-linearity
        (:m1, 1.0), # mass 1 to 4
        (:m2, 1.0),
        (:m3, 1.0),
        (:m4, 1.0)
    ]

    mass_input = LSQInput(β₀, z);
    mass_model = @withinput mass_input begin
        mass_combinations = [
            m1 + m2 + m3 + m4 # I1
            m1 + m2 + m3 + m4 # I2
            m1 + m2 + m3      # I3
            m1 + m2 + m4      # I4
            m1 + m2           # I5
            m1 + m3 + m4      # I6
            m1 + m4           # I7
            m1 + m3           # I8
            m1                # I9
            m2 + m3 + m4      # I10
            m2 + m3           # I11
            m2 + m4           # I12
            m2                # I13
            m3 + m4           # I14
            m3                # I15
            m4                # I16
        ]
        return [
            @. mass_combinations*(1 - (a - 1.2)*(1/ρ - 1/8000)) - f*(I + A*I^2);
            @. mr*(1 - (a - 1.2)*(1/ρr - 1/8000)) - f*(Ir + A*Ir^2);
            ms - (m1 + m2 + m3 + m4)
        ]
    end
    output = solve(mass_model, mass_input; tol=1e-9);
    
    @test output.χ²   ≈ 8.070939437514037
    @test output.β[1] ≈ 1.0000018199321477 # estimated f
    @test output.β[3] ≈ 100.00577273693919 # estimated m1
    @test uncertainty_β(output)[1] ≈ 1.9030092538139585e-7 # estimated std. uncertainty f
    @test uncertainty_β(output)[3] ≈ 1.0989125630912822e-5 # estimated std. uncertainty m1
end