
function evolvePopEnsemble(actSiteNodes::IntVector,
                           assayList::Vector{String},
                           np::NetParams,
                           ep::EvoParams,
                           seeds::Vector{UInt},
                           numReps::Int;
                           parallel::Bool=true)

    if parallel
        out = pmap(x -> evolvePop(np, ep, actSiteNodes, assayList, x), seeds)
    else
        out = map(x -> evolvePop(np, ep, actSiteNodes, assayList, x), seeds)
    end

    seqsList = getindex.(out, 1)
    structsList = getindex.(out, 2)
    adjMatList = getindex.(out, 3)
    tablesList = getindex.(out, 4)
    perturbationsList = getindex.(out, 5)
    return seqsList, structsList, adjMatList, tablesList, perturbationsList
end


function computeEnergyDMSEns(seqsList, structsList, adjMatList, tablesList, perturbationsList)

    function kernel(seqs, r, A, K, pokesList)
        numNodes, popSize = size(seqs)
        numTypes = size(K, 1)
        pokes = pokesList[3]
        dms = Array{Float32, 4}(undef,numTypes, numNodes, 2, popSize)
        for j in 1:popSize
            seq = seqs[:,j]
            Q = buildNetwork(seq, r, A, K)
            dms[:,:,:,j] .= Float32.(computeEnergyDMS(Q, K, pokes)) 
        end
        return dms
    end
    
    numReps = length(seqsList)
    numNodes, popSize = size(seqsList[1])
    numTypes = size(tablesList[1],1)
    megadms = Array{Float32, 5}(undef, numTypes, numNodes, 2, popSize, numReps)
    out = pmap( kernel, seqsList, structsList, adjMatList, tablesList, perturbationsList)
    for k in 1:numReps; megadms[:,:,:,:,k] .= out[k] end
    return megadms
end

function countMutSensitivityEns(dms; x=1e-2)
    @assert ndims(dms) == 5
    numTypes, numNodes, numStrains, popSize, numReps = size(dms)
    numSensMutsArray = Matrix{Int64}(undef, popSize, numReps)
    numSensPosArray = similar(numSensMutsArray)
    for i in 1:popSize, j in 1:numReps
        numSensMuts, numSensPos = countMutSensitivity(dms[:,:,:,i,j]; x)
        numSensMutsArray[i,j] = numSensMuts
        numSensPosArray[i,j] = numSensPos
    end
    return numSensMutsArray, numSensPosArray
end















