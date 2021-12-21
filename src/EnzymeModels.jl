module EnzymeModels

using ReactionNetworks: AbstractReactionNetwork
import ReactionNetworks.stoichiometric_matrix
export Enzyme, EnzymeMM

struct Enzyme <: AbstractReactionNetwork{4,3}
    k1::Float64
    k2::Float64
    k3::Float64
end
function stoichiometric_matrix(::Enzyme)
    return [
        [-1, -1, 1, 0],
        [1, 1, -1, 0],
        [0, 1, -1, 1],
    ]
end
function (rn::Enzyme)(v::AbstractVector{Float64}, x::AbstractVector{Int})
    v[1] = rn.k1*x[1]*x[2]
    v[2] = rn.k2*x[3]
    v[3] = rn.k3*x[3]
    return nothing
end

struct EnzymeMM <: AbstractReactionNetwork{2,1}
    k1::Float64
    k2::Float64
    k3::Float64
    e0::Int64
end
function stoichiometric_matrix(::EnzymeMM)
    return [
        [-1, 1],
    ]
end
function (rn::EnzymeMM)(v::AbstractVector{Float64}, x::AbstractVector{Int})
    K = (rn.k2 + rn.k3)/rn.k1
    v[1] = rn.k3 * min(rn.e0, x[1]) * x[1] / (K + x[1])
    return nothing
end

end
