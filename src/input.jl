"""
Build index from input parameters, e.g.

```julia
z = [
    (:a, [1.0, 1.9], [0.1, 0.1]),
    (:b, 1.0, 0.1),
    (:c, [2.9, 5.0], [0.1, 0.05])
]
```

is turned into 

```julia
Dict([
    :a => 1:2,
    :b => 3,
    :c => 4:5,
    :a_1 => 1,
    :a_2 => 2,
    :c_1 => 4,
    :c_2 => 5
])
```

indicating at which index the parameters are stored in an array
generated with `_build_var_array()`.
"""
function _build_index(params)
    i = 1
    indices = []
    for (k, v...) in params
        j = i + length(v[1]) - 1
        if i == j
            push!(indices, k => i)
        else
            push!(indices, k => i:j)
            for (m, n) in enumerate(i:j)
                # store also the location of indexed variables
                push!(indices, Symbol(k, "_$m") => n)
            end
        end
        i += (j - i + 1)
    end
    Dict(indices)
end

"""
Build array from input parameters, e.g.

    ```julia
    z = [
        (:a, [1.0, 1.9], [0.1, 0.1]),
        (:b, 1.0, 0.1),
        (:c, [2.9, 5.0], [0.1, 0.05])
    ]
    ```
    
    yields
    
    ```julia
    # read values at position 2 of tuples in `z`
    _build_var_array(z, 2) == [1.0, 1.9, 1.0, 2.9, 5.0]
    # read uncertainties at position 3 of tuples in `z`
    _build_var_array(z, 3) == [0.1, 0.1, 0.1, 0.1, 0.05]
    ```
"""
_build_var_array(var_dict, n) = reduce((acc, elm) -> append!(acc, elm[n]), values(var_dict); init=Float64[])


struct LSQInput
    β::Vector{Float64}
    z::Vector{Float64}
    Σ::Matrix{Float64}
    βindex::Dict{Symbol, Union{Int, UnitRange{Int}}}
    zindex::Dict{Symbol, Union{Int, UnitRange{Int}}}
end

"""
    function LSQInput(β, z)

Compile the measured and guessed parameters for a constrained least squares problem.

`β` is a vector of pairs, each holding a symbol and initial guess for (an) unknown parameter(s).
The guess can be a scalar or vector of floats.

`z` is a vector of 3-tuples, each holding a symbol, the measured value and standard uncertainty
of a measured parameter. Values and uncertainties can be floats or vectors of floats.

# Example

For a linear fit, one could provide:

```julia
# measured parameters
z = [
    (:x, [1.0, 1.9, 3.2, 3.8], [0.1,  0.1, 0.2, 0.3]),
    (:y, [2.9, 5.0, 6.9, 9.1], [0.1, 0.05, 0.1, 0.2])
]

# unknown parameters, a = slope, b = offset
β₀ = [
    (:a,  1.0),
    (:b,  1.0)
]

input = LSQInput(z, β₀)
```
"""
function LSQInput(β, z)
    uz = _build_var_array(z, 3) # third element in tuple = uncertainty    
    return LSQInput(
        _build_var_array(β, 2),
        _build_var_array(z, 2), # second element in tuple = value
        Diagonal(uz.^2) |> Matrix,
        _build_index(β),
        _build_index(z)
    )
end


"""
    function covariance!(input::LSQInput, x::Symbol, y::Symbol, value::Float64)

Set covariance of two input parameters stored in `input`.

If parameter is part of an array, use underscore to indicate which 
covariance shall be set, e.g.

```julia
function covariance!(input, :x_1, :x_2, 0.1)
```

sets the covariance of the first and second element in
the array parameter `x`.
"""
function set_covariance!(input::LSQInput, x::Symbol, y::Symbol, value::Float64)
    i = input.zindex[x]
    j = input.zindex[y]
    
    (!isa(i, Integer) || !isa(j, Integer)) && ArgumentError("Covariance can only be set on individual variables.")
      
    input.Σ[i, j] = input.Σ[j, i] = value

    return nothing
end