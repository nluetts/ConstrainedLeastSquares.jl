module ConstrainedLeastSquares

using Espresso:      subs
using LinearAlgebra: Diagonal
using NamedArrays:   NamedArray
using PrettyTables:  PrettyTables
using Zygote:        jacobian

include("input.jl")
include("output.jl")
include("model.jl")
include("solve.jl")

export
    LSQInput,
    LSQOutput,
    @withinput_debug,
    @withinput,
    correlation_β,
    covariance_β,
    covariance!,
    names_z,
    names_β,
    print_corr_β,
    print_cov_β,
    solve,
    uncertainty_β

end # module
