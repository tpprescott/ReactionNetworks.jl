module ReactionNetworks

using Random
export ReactionChannels, AbstractReactionNetwork, simulate

struct ReactionChannel
    X::Vector{Float64}
    function ReactionChannel()
        new(randexp(1))
    end
end

function Base.reset(r::ReactionChannel) 
    empty!(r.X)
    push!(r.X, randexp())
    return nothing
end
function Base.iterate(r::ReactionChannel, i::Int=1)
    length(r.X)==i && push!(r.X, randexp())
    r.X[i], i+1
end
Base.IteratorSize(::Type{ReactionChannel}) = Base.IsInfinite()
Base.eltype(::Type{ReactionChannel}) = Float64

abstract type AbstractReactionNetwork{N, M} end
function stoichiometric_matrix(::AbstractReactionNetwork) end
function reaction_rates!(v::AbstractVector{Float64}, RN::AbstractReactionNetwork{N,M}, x::AbstractVector{Int}) where {N,M} where T
    RN(v, x)
    return nothing
end
function reaction_rates(RN::AbstractReactionNetwork{N,M}, x::AbstractVector{Int}) where {N,M} where T
    v = zeros(M)
    reaction_rates!(v, RN, x)
    return v
end

ReactionChannels{M} = NTuple{M, Union{Nothing,ReactionChannel}}
ReactionChannels{M}() where M = tuple(Iterators.repeated(nothing, M)...)

function simulate(
    RN::AbstractReactionNetwork{N,M},
    x0::AbstractVector{Int},
    tspan::Tuple{Float64,Float64},
    channels::ReactionChannels{M} = ReactionChannels{M}(),
) where {N,M}
    
    length(x0)==N || error("Incorrect dimensions")
    S = stoichiometric_matrix(RN)
    t0, t1 = tspan
    t0<t1 || error("Need t0<t1 in tspan")

    x = copy(x0)
    t = t0

    xsav = Vector{Int}()
    tsav = Vector{Float64}()
    csav = tuple(((isnothing(chn) ? ReactionChannel() : chn) for chn in channels)...)

    distance_generators = broadcast(Iterators.Stateful, csav)
    distance_to_next_jump = [first(dg) for dg in distance_generators]
    speeds = reaction_rates(RN, x)
    times_to_next_jump = distance_to_next_jump ./ speeds
        
    while true
        push!(tsav, t)
        append!(xsav, x)
        
        dt, jumper_idx = findmin(times_to_next_jump)
        t += dt
        t>t1 && break
        @. distance_to_next_jump -= dt*speeds
        distance_to_next_jump[jumper_idx] = first(distance_generators[jumper_idx])

        x .+= S[jumper_idx]
        
        reaction_rates!(speeds, RN, x)
        @. times_to_next_jump = distance_to_next_jump / speeds
    end
    
    return tsav, reshape(xsav,N,:), csav
end

include("EnzymeModels.jl")
export EnzymeModels

end # module
