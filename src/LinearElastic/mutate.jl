

function updateSpringMatrix!(Q::Network, K::FloatMatrix, nodeIndex::Int, newType::Int)
    for i in 1:size(Q.SM,1)
        if Q.SM[i,nodeIndex] != 0.0
            otherNodeType = Q.seq[i]
            k = K[newType, otherNodeType]
            Q.SM[i,nodeIndex] = Q.SM[nodeIndex,i] = k
        end
    end
    return nothing
end


function mutateNode!(Q, K;
                     nodeIndex=rand(1:length(Q.seq)),
                     newType=rand(setdiff(1:size(K,1), Q.seq[nodeIndex]))) # allocates memory
    # Mutates 1 node of the network to a new type and updates the physics,
    # Defaults to random mutation if nodeIndex and newType not given. 

    # update Sequence
    Q.seq[nodeIndex] = newType

    # update Spring Matrix
    updateSpringMatrix!(Q, K, nodeIndex, newType)

    # update Hessian
    updateHessian!(Q, nodeIndex) 
    return nothing
end


function mutateAtRate!(Q::Network,
                       K::FloatMatrix,
                       μ::Number)

    # Mutate each position of sequence with probability μ.
    didMutate = false
    numTypes = size(K,1)
    for i in eachindex(Q.seq)
        if rand() < μ
            didMutate = true
            newType=randNewType(numTypes, Q.seq[i])
            Q.seq[i] = newType
            updateSpringMatrix!(Q, K, i, newType)
            updateHessian!(Q, i) 
        end
    end
    return didMutate
end


function randNewType(numTypes, oldType)
    a = rand( 1:(numTypes-1) )
    newType = ((oldType-1 + a) % numTypes) +1
    return newType
end

