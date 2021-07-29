"""
Build D matrix used in LSQ solver algorithm.
"""
function _buildD(func, β, ζ, Σ)

    Nf = func(β, ζ) |> length
    Nβ = length(β)
    Nζ = length(ζ)

    ∇fβ, ∇fζ = jacobian(func, β, ζ)

    return [
        zeros(Nβ, Nβ) zeros(Nβ, Nζ)      transpose(∇fβ);
        zeros(Nζ, Nβ)        inv(Σ)      transpose(∇fζ);
                  ∇fβ           ∇fζ       zeros(Nf, Nf)
    ]
end

"""
Apply solver algorithm to model `func` and data `input`.
"""
function solve(func::Function, input::LSQInput; iter_max=100, tol=1e-8)
    
    βₗ, z, Σ = input.β, input.z, input.Σ
    ζₗ = copy(z)
    Nf = func(βₗ, ζₗ) |> length
    Nβ = length(βₗ)
    Nζ = length(ζₗ)

    vₗ = [βₗ; ζₗ; zeros(Nf)] # iteratively updated vector
    
    # previous chi-square
    χ²ₗ₋₁ = χ² = nothing
    
    for l in 1:iter_max
        D = _buildD(func, βₗ, ζₗ, Σ)
        f = func(βₗ, ζₗ)
        vₗ += D \ [zeros(Nβ); Σ\(z - ζₗ); -f]
        βₗ = vₗ[1:Nβ]
        ζₗ = vₗ[Nβ+1:Nβ+Nζ]
        χ² = transpose(z - ζₗ)*inv(Σ)*(z - ζₗ)
        println("iteration $l: χ² = $(χ²)")
        if !isnothing(χ²ₗ₋₁) && abs((χ²ₗ₋₁ - χ²)/χ²) < tol
            break
        else
            χ²ₗ₋₁ = χ²
        end
    end
    
    invD = inv(_buildD(func, βₗ, ζₗ, Σ))
    Σβ = invD[1:Nβ, 1:Nβ]
    Σζ = invD[Nβ+1:Nβ+Nζ, Nβ+1:Nβ+Nζ]
    Σβζ= invD[1:Nβ, Nβ+1:Nβ+Nζ]
        
    return LSQOutput(βₗ, ζₗ, Σβ, Σζ, Σβζ, χ², Nf - Nβ, input.βindex, input.zindex)
end