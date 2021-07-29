# ConstrainedLeastSquares

This package implements a general least squares method described in the paper *Evaluation of Measurements by the Method of Least Squares* by 
Lars Nielsen (2001).

The input of the least squares method is a set of `m` measured parameters `z₁ … zₘ` and their respective uncertainties and covariances (covariance matrix `Σ`), a guess for a set of `k` unknown parameters `β₁ … βₖ` and a set of `n` constraints which define the measurement model:

```
f₁(β, ζ) = 0
f₂(β, ζ) = 0
    …
fₙ(β, ζ) = 0
```

The chi-square function to minimize is:

```
χ²(ζ, z) = (z - ζ)ᵀ Σ⁻¹(z - ζ)
```

`ζ` is a set of `m` parameters that are optimized such that the constraint functions are satisfied while the deviation of `ζ` to `z`, weighted by the uncertainty and covariances of `z`, is minimal.
The output of the method is an estimate of optimal parameters `β'` and `ζ'` and their uncertainties and covariance matrices.

## Installation

The package is an early prototype and not registered in the Julia registry. Install it in Julia's package mode *via* (press `]` in the Julia REPL):

```julia
add https://github.com/nluetts/ConstrainedLeastSquares.jl
```

## Usage

### Straight line fit

Fitting a straight line may be done as follows:

```julia
using ConstrainedLeastSquares
using Random

Random.seed!(42)

# Some toy data:
N = 4
xs = collect(Float64, 1:N) .+ randn(N)/10
ys = xs * 2.2 .+ 3.1 .+ randn(N)/10

```
Parameters are defined as vector of tuples:

```julia
# measured parameters: symbol, value(s), uncertainty/uncertainties
z = [
    (:x, xs, 0.05xs),
    (:y, ys, 0.05ys)
]

# unknown parameters: symbol, guess(es)
β₀ = [
    (:a,  1.0),
    (:b,  1.0)
]
```
Input parameters are compiled into an `LSQInput` object:
```julia
linear_input = LSQInput(β₀, z)
```

Covariance can be set for individual input parameters.
The ith element of an array parameter is accessed via `_i`.

```julia
set_covariance!(linear_input, :y_3, :y_2, 0.12)
set_covariance!(linear_input, :y_3, :y_4, 0.14)
```

The `@withinput` macro declares a model for a specific input dataset.
The model expression is wrapped in an `begin ... end` block.
It has to return the values of the constraint functions.

```julia
linear_model = @withinput linear_input begin
    return @. a*x + b - y
end
```

The least squares problem is solved with the `solve()` function:

```julia
output = solve(linear_model, linear_input);
```

```
iteration 1: χ² = 0.31650639079113685
iteration 2: χ² = 0.23051235592939925
iteration 3: χ² = 0.2305741386176354
iteration 4: χ² = 0.2305741992074792
iteration 5: χ² = 0.23057419925831718
```

To get the estimates of the unknown parameters,
simply access the field `β` (similar for `ζ`, `Σζ` and `χ²`):


```julia
output.β
```

```
2-element Vector{Float64}:
 2.1226812851114447
 3.2970009738685717
```

Standard uncertainties are retrieved *via*:

```julia
uncertainty_β(output)
```

```
2-element Vector{Float64}:
 0.23619689830895196
 0.42492929738409446
```

To print the correlation matrix, do:

```julia
print_corr_β(output)
```

```
correlation of β parameters
┌───┬───────────┬───────────┐
│   │         a │         b │
├───┼───────────┼───────────┤
│ a │       1.0 │ -0.848899 │
│ b │ -0.848899 │       1.0 │
└───┴───────────┴───────────┘
```

Similarly, `print_cov_β` displays the covariance matrix.

The matrices are retrieved with the `covariance_β()` and `correlation_β()` functions, e.g.:

```julia
covariance_β(output)
```

```
2×2 Named Matrix{Float64}
A ╲ B │         :a          :b
──────┼───────────────────────
:a    │   0.055789  -0.0852015
:b    │ -0.0852015    0.180565
```

### Example from Lars Nielsen's paper

The first example from Lars Nielsen's paper (*Calibration of an analytical balance*) can be formulated as follows:

```julia
using ConstrainedLeastSquares

# measured parameters
z = [
    (:ms, 199.988816, 0.000010), # known sum of weights in g
    (:mr, 199.999043, 0.000008), # reference weight in g
    (:ρr, 7833.01, 0.29), # densities in kg/m³
    (:ρ,  7965.76, 0.71), 
    (:a, 1.1950, 0.0035), # air density
    (:I, [199.988617      # scale indication in "div"
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
    (:m1, 1.0), # unknown mass 1 to 4
    (:m2, 1.0),
    (:m3, 1.0),
    (:m4, 1.0)
]

mass_input = LSQInput(β₀, z)

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

output = solve(mass_model, mass_input; tol=1e-9, iter_max=10);
```

```
iteration 1: χ² = 0.0003222308915394452
iteration 2: χ² = 8.07786170097814
iteration 3: χ² = 8.07093712244527
iteration 4: χ² = 8.070939442718299
iteration 5: χ² = 8.070939434553656
iteration 6: χ² = 8.070939437514037
```

```julia
output.β
```

```
6-element Vector{Float64}:
   1.0000018199321477
  -4.2474329421029445e-9
 100.00577273693919
  50.00796185700379
  24.97860013758237
  24.996475404873838
```