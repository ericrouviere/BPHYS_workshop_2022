
function buildStructure(W::Int, L::Int, offset::Number)

    @assert isodd(L)
    numLayers = Int((L-1)/2)
    stride = 2W-1
    numNodes = stride * numLayers + W - 1
    
    v1 = [cos(pi/6), sin(pi/6)] # basis vector 1 for triangular matrix
    v2 = [0,1] # basis vector 2 for triangular matrix
    V = [v1 v2] # change of basis matrix. From lattice basis to standard basis

    layer = (V * [zeros(W-1) 0:W-2 ; ones(W) -1:W-2]')'
    layer_x = layer[:,1]
    layer_y = layer[:,2]

    px = zeros(stride*numLayers + W-1 )
    py = zeros(stride*numLayers + W -1)
    for i in 1:numLayers    
        px[(i-1)*stride+1:(i)*stride] .= layer_x .+ (i-1)*2*sin(pi/3)
        py[(i-1)*stride+1:(i)*stride] .= layer_y
    end

    px[end-W+2:end] .= layer_x[1:W-1] .+ (numLayers)*2*sin(pi/3)
    py[end-W+2:end] .= layer_y[1:W-1]

    # build Adjacency Matrix
    A = zeros(numNodes,numNodes);
    for i in 2:numNodes
        for j in 1:i-1
            dx = px[i] - px[j]
            dy = py[i] - py[j]
            dist = sqrt(dx^2 + dy^2)
            if dist < 1.1
                A[i,j] = A[j,i] = 1.
            end
        end
    end

    x = copy(px) .+ offset .* randn(numNodes)
    y = copy(py) .+ offset .* randn(numNodes)
    r = xy2r([x y])
    return r, A
end


function buildSpringTable(numTypes::Int, k_min::Number, k_max::Number; binary=true, numSoft=0)
    # generate a Matrix K_ij where i and j are nodes types and
    # and K_ij is the spring constant that connect nodes of type i
    # and type j.

    if binary # make all but numSoft entries stiff
        K = k_max*ones(numTypes,numTypes)
        inds = findall(LowerTriangular(K) .> 0)
        rp = randperm(length(inds))
        softInds = inds[rp[1:numSoft]]
        K[softInds] .= k_min
    else # pull K from log uniform distribution
        log_diff = log(10, k_max) - log(10, k_min)
        K = 10 .^ (rand(numTypes, numTypes) .* log_diff .+ log(10,k_min))
    end
    K = Float64.(Symmetric(K, :L))
    return K
end

function seq2Springs(seq, K, A)
    # convert a sequence into a matrix that stores
    # the spring constant between node i and node j.
    SM = zeros(size(A))
    for j in axes(A,2), i in axes(A,1)
        if A[i,j] != 0
            SM[i,j] = K[seq[i], seq[j]]
        end
    end
    return SM
end

function randSeq(r, K)
    # generate a random sequence of node types.
    numNodes = Int(length(r)/2)
    numTypes = size(K, 1)
    seq = rand(1:numTypes, numNodes)
    return seq
end

function buildStructAdjMatTable(W::Int, L::Int, offset::Number, numTypes::Int,
                      k_min::Number, k_max::Number, binaryK::Bool, numSoft::Int)
    
    r, A = buildStructure(W,L,offset);
    K = buildSpringTable(numTypes, k_min, k_max; binary=binaryK, numSoft)
    return r, A, K
end

buildStructAdjMatTable(np::NetParams) = buildStructAdjMatTable(np.W, np.L, np.offset,
                             np.numTypes, np.k_min, np.k_max, np.binaryK, np.numSoft)

function buildNetwork(seq::IntVector, r::FloatVector, A::NumMatrix, K::FloatMatrix)
    SM = seq2Springs(seq, K, A);
    H = computeHessian(r, SM)
    b = Vector{eltype(H)}(undef,length(r))
    #c = similar(b)
    @assert 2*length(seq) == length(r) == 2*size(SM,1) == size(H,1) 
    
    Q = Network(seq, copy(r), SM, H, similar(H), b, similar(b))
    return Q
end

function buildNetwork(r::FloatVector, A::NumMatrix, K::FloatMatrix)
    seq = randSeq(r, K)
    return buildNetwork(seq, r, A, K)
end

function buildNetwork(W::Int, L::Int, offset::Number, numTypes::Int,
                      k_min::Number, k_max::Number, binaryK::Bool, numSoft::Int)
    # build network, Spring table and Adjacency matrix
    r,A,K = buildStructAdjMatTable(W, L, offset, numTypes, k_min, k_max, binaryK, numSoft)
    Q = buildNetwork(r, A, K)
    return Q, K, A
end


function buildPopulation(r::FloatVector, A::NumMatrix, K::FloatMatrix, P::Int; diffSeqs::Bool=true)
    if diffSeqs
        nets = [buildNetwork(r, A, K) for i in 1:P]
    else
        Q = buildNetwork(r, A, K)
        nets = [deepcopy(Q) for i in 1:P]
    end
    return nets
end

function buildPopulation(W::Int, L::Int, offset::Number, numTypes::Int,
                      k_min::Number, k_max::Number, binaryK::Bool, numSoft::Int,
                      P::Int; same=true)
    
    r, A = buildStructure(W,L,offset);
    K = buildSpringTable(numTypes, k_min, k_max; binary=binaryK, numSoft)
    return buildPopulation(r,A,K,P; same)
end

function makePokesListFlux(r::FloatVector, actSiteNodes::IntVector; numStrains=2)
    # this  is a temporary function that is very specific. VERY BAD.
    @assert numStrains == 2 # this is temporary

    site = sort([actSiteNodes * 2; actSiteNodes * 2 .- 1])
    sites1 = [copy(site)]
    sites2 = [copy(site)]
    sites3 = [copy(site), copy(site)]
    
    strains = createNormalStrains(numStrains, r, actSiteNodes)
    strains1 = [copy(strains[1])]
    strains2 = [copy(strains[2])]
    strains3 = strains
    
    pokesList =[ Pokes(sites1, strains1),
                 Pokes(sites2, strains2),
                 Pokes(sites3, strains3)]
    return pokesList
end


