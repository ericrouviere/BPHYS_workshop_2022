# convert between postions as a vector r = [x1,y1,x2,y2,...]
# to position as a matrix xy = [x1 y1; x1 y1;..]
r2xy(x::Vector) = Array(reshape(x, 2, :)')

# inverse of r2xy
xy2r(xy::Matrix) = xy'[:]

function networkCopy!(Q_des::Network, Q_src::Network)
    # copy one network into another without allocating memory.  
    Q_des.seq .= Q_src.seq
    Q_des.r .= Q_src.r
    Q_des.SM .= Q_src.SM
    Q_des.H .= Q_src.H
    return nothing
end

function networkCopy!(nets::Vector{Network}, des_inds::IntVector, src_inds::IntVector)
    @assert length(des_inds) == length(src_inds)
    for i in eachindex(des_inds)
        j = src_inds[i]
        k = des_inds[i]
        networkCopy!(nets[k], nets[j])
    end
    return nothing
end

function pop2Seqs(nets::Vector{Network})
    Q = nets[1]
    numNodes = length(Q.seq)
    P = length(nets)
    seqs = Matrix{eltype(Q.seq)}(undef, numNodes, P)
    for j in 1:P
        seqs[:,j] .= nets[j].seq
    end
    return seqs
end

function convertSettings2NetParams(S::Dict)
    @assert isodd(S["L"])
    @assert S["offset"] >= 0
    @assert S["numTypes"] >= 2
    @assert S["k_min"] <= S["k_max"]
    @assert 0 <= S["numSoft"] <= S["numTypes"]
    return NetParams(S["W"], S["L"], Float64(S["offset"]), S["numTypes"],
                     Float64(S["k_min"]), Float64(S["k_max"]),
                     S["binaryK"], S["numSoft"])
end

function convertSettings2EvoParams(S::Dict)
    @assert 0 <= S["μ"] <= 1
    @assert 0 <= S["a"] <= 1
    return EvoParams(S["P"], Float64(S["a"]),S["N"], Float64(S["μ"]), S["τ"],
                     S["sampleLastHalf"],S["numSeqs2Keep"],S["diffSeqs"])
end

function compressPokesList(pokesList::Vector{Pokes})
# this function is temporary. will change later
    strains = pokesList[3].strains
    sites = pokesList[3].sites
    return sites, strains
end

function expandPokesList(compressedPokesList)
    # this just aweful code
    sites = compressedPokesList[1]
    strains = compressedPokesList[2]
    pokes1 = Pokes(deepcopy([sites[1]]), deepcopy( [strains[1]] ) )
    pokes2 = Pokes(deepcopy([sites[2]]), deepcopy( [strains[2]] ) )
    pokes3 = Pokes(deepcopy(sites), deepcopy(strains))
    return [pokes1, pokes2, pokes3]
end

function saveEvoData(filename, seqsList, structsList, adjMatList, tablesList, perturbationsList, seeds)
    # save the data from evolvePopEnsemble
    sequences = cat(seqsList..., dims=3)
    structures = cat(structsList..., dims=2)
    tables = cat(tablesList..., dims=3)
    A = adjMatList[1]
    perturbations = compressPokesList.(perturbationsList)
    
    jldopen(filename,"w") do file
        write(file, "sequences", sequences)
        write(file, "structures", structures)
        write(file, "A", A)
        write(file, "tables", tables)
        write(file, "perturbations", perturbations)
        write(file, "seeds", seeds)
    end
end

function getSeqs(nets::Vector{Network}, inds)
    return [copy(nets[i].seq) for i in inds]
end


































