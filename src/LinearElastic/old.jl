#function computeHessian(r::Vector{Float64},
#                        SM::Matrix{Float64})
#    # compute Hessian matrix of Hamiltonian
#    # returns hessian
#    N = size(SM,1)
#    H = zeros(2N, 2N)
#    for i in 1:N
#        sxx = 0.0
#        syy = 0.0 
#        sxy = 0.0
#        i2 = 2i
#        xi = r[i2-1]
#        yi = r[i2]
#        for j in 1:N
#            k = SM[i,j]
#            if k > 0 && i != j
#                j2 = 2j
#                dx = xi - r[j2-1]
#                dy = yi - r[j2]
#                c = -k / (dx^2 + dy^2)
#                Exx = c * dx^2 
#                Eyy = c * dy^2 
#                Exy = c * dx*dy 
#                sxx -= Exx
#                syy -= Eyy
#                sxy -= Exy
#                H[i2-1,j2-1] = Exx
#                H[i2,j2] = Eyy
#                H[i2,j2-1] = Exy
#                H[i2-1,j2] = Exy
#            end
#        end
#        H[i2-1,i2-1] = sxx
#        H[i2,i2] = syy
#        H[i2,i2-1] = sxy
#        H[i2-1,i2] = sxy
#    end
#    return H
#end


#function computeHessian!(H::Matrix,
#                         r::Vector,
#                         SM::Matrix)
#    # compute the hessian matrix of the network
#    N = size(SM,1)
#    for i in 1:N
#        sxx = 0.0
#        syy = 0.0 
#        sxy = 0.0
#        i2 = 2i
#        xi = r[i2-1]
#        yi = r[i2]
#        for j in 1:N
#            
#            j2 = 2j
#            # fill H with zeros
#            H[i2-1,j2-1] = 0.0
#            H[i2,j2-1] = 0.0
#            H[i2-1,j2] = 0.0
#            H[i2,j2] = 0.0
#
#            k = SM[i,j]
#            if k > 0 && i != j
#                dx = xi - r[j2-1]
#                dy = yi - r[j2]
#                c = -k / (dx^2 + dy^2)
#                Exx = c * dx^2 
#                Eyy = c * dy^2 
#                Exy = c * dx*dy 
#                sxx -= Exx
#                syy -= Eyy
#                sxy -= Exy
#                H[i2-1,j2-1] = Exx
#                H[i2,j2-1] = Exy
#                H[i2-1,j2] = Exy
#                H[i2,j2] = Eyy
#            end
#        end
#        H[i2-1,i2-1] = sxx
#        H[i2,i2-1] = sxy
#        H[i2-1,i2] = sxy
#        H[i2,i2] = syy
#    end
#    return nothing
#end

#function sampleSequences!(seqBucket, nets, sampleLastHalf, t, gen2StartSampling, numSeqs2KeepPerGen)
#    # sample sequence from the last half of evolution.
#    if sampleLastHalf && t >= gen2StartSampling
#        inds = randperm(P)[1:numSeqs2KeepPerGen]
#        # [copy(nets[i].seq) for i in inds]
#        selectedSeqs = [copy(nets[i].seq) for i in inds]
#        append!(seqBucket, selectedSeqs)
#    end
#end


#function samplingDetails(N, P, numSeqs2Keep)
#    # some book keeping for sampling sequence from the population
#    # in evolvePop.
#    numGen2Sample = Int(floor(N / 2))
#    maxSeqsSampleSize = P * numGen2Sample
#    if numSeqs2Keep > maxSeqsSampleSize
#        println("You are trying to keep too many sequences\n"*
#                "Setting numSeqs2Keep = maxSeqsSampleSize")
#        numSeqs2Keep = maxSeqsSampleSize
#    end
#    numSeqs2KeepPerGen = Int(ceil(numSeqs2Keep / numGen2Sample))
#    gen2StartSampling = Int(ceil(N/2))
#    return numGen2Sample,maxSeqsSampleSize,numSeqs2KeepPerGen, gen2StartSampling,numSeqs2Keep
#end

