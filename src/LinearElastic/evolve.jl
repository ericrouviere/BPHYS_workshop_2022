function evolveMC(Q::Network,
                  K::FloatMatrix,
                  pokes::Pokes,
                  assay::String,
                  N::Int,
                  ζ::Real;
                  ΔT=1)
    # evolve using monte carlo evolution
    # WARNING this code is slow due to heavy allocation.

    Q = deepcopy(Q)
    Q_wt = deepcopy(Q)
    ϕ = computeFitness(Q, pokes, assay)
    times, fitnesses, networks = Int[], Float64[], Network[]
    ΔT == 0 && (ΔT = N-1) 
    for t in 0:N-1
        networkCopy!(Q_wt, Q)
        mutateNode!(Q, K)
        ϕ_mut = computeFitness(Q, pokes, assay)
        # selection
        if rand() < exp(ζ*(ϕ_mut-ϕ))
            ϕ = ϕ_mut
        else
            networkCopy!(Q,Q_wt)
        end
        # recording:
        if mod(t, ΔT) == 0
            push!(times, t)
            push!(fitnesses, ϕ)
            push!(networks, deepcopy(Q))
        end 
    end
    return times, fitnesses, networks
end



function evolvePop!(nets::Vector{Network}, # population of networks
                   K::FloatMatrix, # stiffness interaction table
                   pokesList::Vector{Pokes}, # Strains and Sites for ligand perturbation
                   assayList::Vector{String},
                   P::Int, # Population size
                   a::Number, # selection threshold eg a=0.3 means remove bottom 30%
                   N::Int, # number of generations
                   μ::Number, # mutation rate
                   τ::Int; # enviromental timescale
                   sampleLastHalf::Bool=true,
                   numSeqs2Keep::Int=length(nets)) 

    # Evolve a popluation of sequences under a fluctuating selection
    @assert 0 <= μ <= 1
    @assert 0 <= a <= 1

    P, numNodes = length(nets), length(nets[1].seq)
    τ ==0 && (τ=N)
    env, oldenv, epo = 1, 1, 1
    mutatedRecord = Int[]
    numWinners = Int(floor((1-a)*P) + ((1-a)<1) * 1)
    # setting sequences accounting
    sgens, numSeqs2Keep = samplingDetails(N, P, numSeqs2Keep)
    seqs = Matrix{eltype(nets[1].seq)}(undef, numNodes, numSeqs2Keep)
    fits = zeros(P)
    for t in 0:N-1 # evolution loop
        env, oldenv, epo = updateEnvironment(t, τ, env, oldenv, epo) 
        updateFitness!(fits, nets, pokesList, assayList, env, oldenv, mutatedRecord)
        applySection!(nets, fits, P, numWinners)
        mutatedRecord = mutateAtRate!.(nets, Ref(K), μ)
        sampleSequences!(seqs, nets, sgens, sampleLastHalf, t)
    end

    if !sampleLastHalf
        seqs = hcat([copy(nets[i].seq) for i in eachindex(nets)]...)
    end
    return seqs
end


function evolvePop(r::FloatVector,
                   A::NumMatrix,
                   K::FloatMatrix, # stiffness interaction table
                   pokesList::Vector{Pokes}, # Strains and Sites for ligand perturbation
                   assayList::Vector{String},
                   ep::EvoParams,
                   seed::UInt)

    Random.seed!(seed)
    nets0 = buildPopulation(r,A,K,ep.P; ep.diffSeqs)
    seqs = evolvePop!(nets0, K, pokesList, assayList, ep.P,
                           ep.a, ep.N, ep.μ, ep.τ; ep.sampleLastHalf, ep.numSeqs2Keep)
    return seqs
end


function evolvePop(np::NetParams,
                   ep::EvoParams,
                   actSiteNodes::IntVector,
                   assayList::Vector{String},
                   seed::UInt)
        r, A, K = buildStructAdjMatTable(np)
        pokesList = makePokesListFlux(r, actSiteNodes)
        seqs = evolvePop(r, A, K, pokesList, assayList, ep, seed)
    return seqs, r, A, K, pokesList
end


function samplingDetails(N,P,numSeqs2Keep)
    # save up to 1 sequence per generation
    
    numGen2Sample = Int(floor(N / 2))
    shift = Int(ceil(N/2))
    if numSeqs2Keep > numGen2Sample
        println("You are trying to keep too many sequences\n"*
                "Setting numSeqs2Keep = numGen2Sample")
        numSeqs2Keep = numGen2Sample
    end
    # the generations to sample 1 sequence
    sgens = randperm(numGen2Sample)[1:numSeqs2Keep] .+ shift .- 1
    sort!(sgens)
    return sgens, numSeqs2Keep
end


function sampleSequences!(seqs, nets, sgens, sampleLastHalf, t)
    # sample sequences from the last half of evolution.
    if sampleLastHalf && t in sgens
        P = length(nets)
        netInd = rand(1:P)
        i = findfirst( x -> x == t, sgens)
        seqs[:,i] .= nets[netInd].seq
    end
end

function applySection!(nets, fits, P, numWinners)
    rankedIndicies = sortperm(fits, rev=true)
    winners = rankedIndicies[1:numWinners]
    losers = rankedIndicies[numWinners+1:end]
    replicates = rand(winners, P-numWinners)
    networkCopy!(nets, losers, replicates)
    @views fits[losers] .= fits[replicates]
end



function updateFitness!(fits, nets, pokesList, assayList, env, oldenv, mutatedRecord)
    # update fitnesses based on the environ
    pokes, assay = pokesList[env], assayList[env]
    if env == oldenv
        fits[mutatedRecord] .= computeFitness.(nets[mutatedRecord], [pokes], assay)
    else
        fits .= computeFitness.(nets, [pokes], assay)
    end
    return nothing
end

function updateEnvironment(t, τ, env, oldenv, epo)
    # this is bad code but it does the job.
    if t % τ == 0
     #   println("τ == 0 -------------------------- ")
        oldenv = env
        env = 3
        epo = 3 - oldenv
    else
        oldenv = env
        env = epo
    end
    return env, oldenv, epo
end

