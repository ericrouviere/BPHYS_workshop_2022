
#struct Network
#    seq::Vector{Int8} # sequence
#    r::Vector{Float64} # structure
#    SM::Matrix{Float64} # spring matrix
#    H::Matrix{Float64} # hessian
#    B::Matrix{Float64} # hessian
#end


struct Network
    seq::Vector{Int8} # sequence
    r::Vector{Float64} # structure
    SM::Matrix{Float64} # spring matrix
    H::Matrix{Float64} # hessian
    B::Matrix{Float64} # hessian
    b::Vector{Float64}
    c::Vector{Float64}
end



struct Pokes
    sites::Vector{Vector{Int}}
    strains::Vector{Vector{Float64}}
end


struct EvoParams
    # stores the evolutionary parameters for evolving populations.
    P::Int # Population size
    a::Float64 # selection threshold eg a=0.3 means remove bottom 30%
    N::Int # number of generations
    μ::Float64 # mutation rate
    τ::Int # enviromental timescale
    sampleLastHalf::Bool # return sequences randomly selected from last half of simulation
    numSeqs2Keep::Int # number of sequences to return
    diffSeqs::Bool # start evolvePop with different sequences or not
end

struct NetParams
    # Stores network parameters for building elastic network.
    W::Int # width of network
    L::Int # Length of network
    offset::Float64 # disorder in structure
    numTypes::Int # number of possible nodes types
    k_min::Float64 # maximum stiffness of a bond
    k_max::Float64 # minimum stiffness of a bond
    binaryK::Bool # is the spring constant table K binary?
    numSoft::Int # number of soft interaction in spring constant table K
end





# aliases for types
const FloatMatrix = AbstractArray{<:AbstractFloat, 2}
const FloatVector = AbstractArray{<:AbstractFloat, 1}
const IntMatrix = AbstractArray{<:Integer,2}
const IntVector = AbstractArray{<:Integer,1}
const NumMatrix = AbstractArray{<:Real,2}
const NumVector = AbstractArray{<:Real,1}
