struct LSQOutput
    β::Vector{Float64}
    ζ::Vector{Float64}
    Σβ::Matrix{Float64}
    Σζ::Matrix{Float64}
    Σβζ::Matrix{Float64}
    χ²::Float64
    dof::Int
    βindex::Dict{Symbol, Union{Int, UnitRange{Int}}}
    zindex::Dict{Symbol, Union{Int, UnitRange{Int}}}
end

_names_ordered(index) = map(x->x[1], sort(collect(index), by=x->x[2]))

"""Return names of measured parameters `z`."""
names_z(output::LSQOutput) = _names_ordered(output.zindex)

"""Return names of unknown parameters `β`."""
names_β(output::LSQOutput) = _names_ordered(output.βindex)

"""Return standard uncertainty of unknown parameters `β`."""
uncertainty_β(output::LSQOutput) = [sqrt(output.Σβ[i, i]) for i in 1:length(output.β)]

"""Return covariance matrix of unknown parameters `β`."""
function covariance_β(output::LSQOutput)
    names = names_β(output)
    return NamedArray(output.Σβ, (names, names))
end

"""Return correlation matrix of unknown parameters `β`."""
function correlation_β(output::LSQOutput)
    cov = covariance_β(output)
    N = size(cov)[1]
    corr = similar(cov)
    for i in 1:N
        for j in 1:N
            corr[i, j] = cov[i, j]/√(cov[i, i]*cov[j, j])
        end
    end
    return corr
end

"""Display covariance matrix of unknown parameters `β`."""
function print_cov_β(output::LSQOutput)
    pretty_table(
        covariance_β(output);
        row_names=names_β(output),
        header=names_β(output),
        title="covariance of β parameters"
    )
end

"""Display correlation matrix of unknown parameters `β`."""
function print_corr_β(output::LSQOutput)
    pretty_table(
        correlation_β(output);
        row_names=names_β(output),
        header=names_β(output),
        title="correlation of β parameters"
    )
end