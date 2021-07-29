"""
Take data in `input` and replace variable symbols in `modelex`
expression with the arrays β and ζ that hold the unknown parameters and
estimates of measured parameters, respectively.
"""
function _makemodel(input::LSQInput, modelex::Expr; debug=false)

    modelname = gensym("model")

    # substitute symbols in model expression with entries in beta and zeta,
    # according to mapping in input
    subssyms = let
        substitutions = [k => :(ζ[$v]) for (k, v) in input.zindex]
        append!(substitutions, [k => :(β[$v]) for (k, v) in input.βindex])
        Dict(substitutions)
    end
    modelex = subs(modelex, subssyms)

    modelcode = quote
        function $modelname(β::Vector{Float64}, ζ::Vector{Float64})
            $modelex
        end
    end

    # print generated model to identify problems
    if debug
        println(modelcode)
    end

    eval(modelcode)
end

"""
    yourmodel = @withinput input begin
        # your model expression
    end

Declare a model using the data in `input`.

All variables that are not assigned to must be present in
the input data, i.e.

```
y = x + 1
```

will throw an error when the model is evaluated if `x` is not included
in the input data.

**The model has to return the values of the constraints**, e.g.
if two constraints are present:

```
y - a*x = 0
a - 2b  = 0
```

The model has to return:

```
return [
    y - a*x
    a - 2b
]
```

See README and tests for further examples.
"""
macro withinput(input::Symbol, modelex::Expr)
    
    modelex = Expr(:quote, modelex)
    esc(quote
        ConstrainedLeastSquares._makemodel($input, $modelex)
    end)
end

"""
Like `@withinput`, but prints the generated model.
"""
macro withinput_debug(input::Symbol, modelex::Expr)
    
    modelex = Expr(:quote, modelex)
    quote
        makemodel($input, $modelex; debug=true)
    end
end