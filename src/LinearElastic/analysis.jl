

function computeDMS(Q::Network, K::FloatMatrix, pokes::Pokes, assay::String)
    # compute deep mutational scan of network
    Q = deepcopy(Q)
    Q_wt = deepcopy(Q)
    numTypes = size(K,1)
    numNodes = length(Q.seq)
    ϕ_wt = computeFitness(Q_wt,pokes,assay)
    dms = zeros(numTypes, numNodes)
    for nodeIndex in 1:numNodes

        type_wt = Q_wt.seq[nodeIndex]
        for newType in 1:numTypes
            newType == type_wt && continue
            mutateNode!(Q, K; nodeIndex, newType)
            ϕ_mut = computeFitness(Q, pokes, assay) 
            dms[newType, nodeIndex]  = ϕ_mut - ϕ_wt
            networkCopy!(Q, Q_wt)
        end
    end
    return dms, ϕ_wt
end


#function countMutSensitivity(dms::FloatMatrix, ϕ_wt::Number; x::Number=1e-2)
#    # count the number of mutatant that change the fitness by at least x.
#    # count the number of positions that have at least 1 mutation that 
#    # changes fitness by at least x.
#
#    rel_dms =  abs.(dms ./ ϕ_wt)
#    max_dms =  maximum(rel_dms, dims=1)
#    numSensMuts = sum(rel_dms .> x)
#    numSensPositions = sum( max_dms .> x)
#    return numSensMuts, numSensPositions
#end


function countMutSensitivity(dms::Array{<:AbstractFloat,3}; x::Number=1e-2)
    # count the number of mutatant that change the fitness by at least x.
    # count the number of positions that have at least 1 mutation that 
    # changes fitness by at least x.
    dms = abs.(dms)
    max_dms = maximum(dms, dims=(1,3))
    numSensMuts = sum(dms .> x)
    numSensPos = sum( max_dms .> x)
    return  numSensMuts, numSensPos
end



function computeEnergyDMS(Q::Network, K::FloatMatrix, pokes::Pokes)
    # compute deep mutational scan of network, measuring the energies, not fitness.

    Q_wt = deepcopy(Q)
    numTypes = size(K,1)
    numNodes = length(Q.seq)
    energies_wt = log.(10,computeEnergiesFast(Q, pokes))
    dms = zeros(numTypes, numNodes, length(energies_wt))

    for nodeIndex in 1:numNodes
        type_wt = Q_wt.seq[nodeIndex]
        for newType in 1:numTypes
            newType == type_wt && continue
            mutateNode!(Q, K; nodeIndex, newType)
            dms[newType, nodeIndex, :]  = log.(10, computeEnergiesFast(Q, pokes)) .- energies_wt
            networkCopy!(Q, Q_wt)
        end
    end
    return dms
end




























