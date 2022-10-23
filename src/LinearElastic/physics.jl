

function buildStructureMatrix(r::FloatVector, A::NumMatrix)
    # construct the structure matrix T    
    numBonds = Int(sum(A)/2)
    numDims = length(r)
    T = zeros(numBonds, numDims)
    k = 0
    AA = UpperTriangular(A)
    @inbounds for i in axes(AA,1), j in axes(AA,2)
        if AA[i,j] > 0
            k +=1
            j2 = 2j
            i2 = 2i
            dx = r[i2-1] - r[j2-1]
            dy = r[i2] - r[j2]
            d_recip = 1 / sqrt(dx^2 + dy^2)
            T[k, 2*(j-1) + 1] = -dx * d_recip
            T[k, 2*(j-1) + 2] = -dy * d_recip 
            T[k, 2*(i-1) + 1] = dx * d_recip
            T[k, 2*(i-1) + 2] = dy * d_recip
        end
    end    
    return T
end


function computeHessian!(H::FloatMatrix,
                         r::FloatVector,
                         SM::FloatMatrix)
    # compute the hessian matrix of the network
    N = size(SM,1)
    @inbounds for j in 1:N
        sxx = 0.0
        syy = 0.0 
        sxy = 0.0
        j2 = 2j
        xj = r[j2-1]
        yj = r[j2]
        @fastmath for i in 1:N
            i2 = 2i
            # fill H with zeros
            H[i2-1,j2-1] = 0.0
            H[i2,j2] = 0.0
            H[i2,j2-1] = 0.0
            H[i2-1,j2] = 0.0
            k = SM[i,j]
            if k > 0 && i != j
                dx = r[i2-1] - xj
                dy = r[i2] - yj
                c = -k / (dx^2 + dy^2)
                Exx = c * dx^2 
                Eyy = c * dy^2 
                Exy = c * dx*dy 
                sxx -= Exx
                syy -= Eyy
                sxy -= Exy
                H[i2-1,j2-1] = Exx
                H[i2,j2] = Eyy
                H[i2,j2-1] = Exy
                H[i2-1,j2] = Exy
            end
        end
        H[j2-1,j2-1] = sxx
        H[j2,j2] = syy
        H[j2,j2-1] = sxy
        H[j2-1,j2] = sxy
    end
    return nothing
end

computeHessian!(Q::Network) = computeHessian!(Q.H,Q.r,Q.SM)

function computeHessian(r::FloatVector, SM::FloatMatrix)
    N = size(SM,1)
    H = Matrix{Float64}(undef, 2N, 2N)
    computeHessian!(H,r,SM)
    return H
end

function updateHessian!(Q::Network, nodeIndex::Int)
    # update hessian for all entries connected to node=node index.
    # SM has to have the new spring constant already.
    r = Q.r
    SM = Q.SM
    H = Q.H
    j = nodeIndex
    N = size(Q.SM,1)

    sxx = 0.0
    syy = 0.0 
    sxy = 0.0
    j2 = 2j
    xj = r[j2-1]
    yj = r[j2]
    @inbounds @fastmath for i in 1:N
        k = SM[i,j]
        if k > 0 && i != j
            i2 = 2i
            dx = r[i2-1] - xj
            dy = r[i2] - yj
            c = -k / (dx^2 + dy^2)
            Exx = c * dx^2 
            Eyy = c * dy^2 
            Exy = c * dx*dy 
            sxx -= Exx
            syy -= Eyy
            sxy -= Exy
            
            #save current entries from i,j
            hxx = H[i2-1,j2-1]
            hyy = H[i2,j2]
            hxy = H[i2,j2-1]

            #update entries in i,j
            H[i2-1,j2-1] = Exx
            H[i2,j2] = Eyy
            H[i2,j2-1] = Exy
            H[i2-1,j2] = Exy
            
            #update entries in j,i
            H[j2-1,i2-1] = Exx
            H[j2,i2] = Eyy
            H[j2,i2-1] = Exy
            H[j2-1,i2] = Exy

            #update entries i,i
            H[i2-1,i2-1] += (hxx - Exx)
            H[i2,i2] += (hyy - Eyy)
            H[i2,i2-1] += (hxy - Exy)
            H[i2-1,i2] += (hxy - Exy)
        end
    end
    H[j2-1,j2-1] = sxx
    H[j2,j2] = syy
    H[j2,j2-1] = sxy
    H[j2-1,j2] = sxy
    return nothing
end



function computeResponse!(B::FloatMatrix,
                          H::FloatMatrix,
                          site::IntVector,
                          strain::FloatVector)
    # compute the response via Riccardo's methods.
    # dr0 is the imposed displacment.
    # H is the hessian matrix.

    numDims = size(H,1)
 #   actSiteDims = sort([actSite * 2; actSite * 2 .- 1])
    @assert length(site) == length(strain)

    # build constraints vector
    b = zeros(numDims) ###########################
    b[site] = strain


    # build matrix A
    B .= -1.0 .* H ##############################
    B[:, site] = I(numDims)[:, site]
        
    # solve Linear Response problem
    y = B \ (H*b) #############################

    # construct full diplacement
    dr = copy(y) #############################3
    dr[site] = strain

    # construct full force
    f = copy(b)
    f[site] = y[site]

    # calculate energy
    E = f ⋅ dr
    
    return E, dr, f
end

function computeResponse(H::FloatMatrix,
                         site::IntVector,
                         strain::FloatVector)
    B = similar(H)
    E, dr, f = computeResponse!(B,H,site,strain)
end

function computeResponse_LM( H::FloatMatrix, site::IntVector, strain::FloatVector)
    # compute the response via Lagrange Multiplier.
    # dr0 is the imposed displacment.
    # H is the hessian matrix.

    numDims = size(H,1)
    @assert length(site) == length(strain)
    n = length(site)

    # build constraints vector
    b = zeros(numDims+n) # 
    b[numDims+1:end] = strain

    B = zeros(numDims+n,numDims+n)
    B[1:numDims,1:numDims] .= H
    B[site+1:end, site] = I(n)
    B[site, site+1:end] = I(n);

    # solve Linear Response problem
    y = B \ b

    # construct full diplacement
    dr = y[1:numDims]

    # compute force
    f = H * dr

    # calculate energy
    E = f ⋅ dr

    return E, dr, f
end




function computeEnergies(Q::Network, pokes::Pokes)
    # compute the energies of applied strains
    n = length(pokes.sites)
    B = similar(Q.H)
    energies = zeros(n)
    for i in 1:n
        site = pokes.sites[i]
        strain = pokes.strains[i]
        E, dr, f = computeResponse!(B, Q.H, site, strain)
        energies[i] = E
    end
    return energies
end


#function computeEnergiesFast(Q::Network, pokes::Pokes)
#    # compute the energies of applied strains
#    
#    # make sure that all the sites are the same. Make more general when needed.
#    @assert all(Ref(pokes.sites[1]) .== pokes.sites)
#    n = length(pokes.sites)
#    @assert n>=1
#    numDims = size(Q.H,1)
#
#    # initial and build arrays
#    H, B = Q.H, Q.B 
#    fillBMatrix!(B,H, pokes.sites[1]) 
#    c = zeros(numDims) 
#    b = similar(c) 
#    energies = zeros(n)
#
#    # first iteracion
#    c .= 0
#    c[pokes.sites[1]] = pokes.strains[1]
#    mul!(b,H,c)
#    prob = LinearProblem(B,b)
#    linsolve = init(prob, alias_A=true)
#    sol = solve(linsolve)
#    energies[1] = sol ⋅ c 
#    
#    
#    for i in 2:n
#        site = pokes.sites[i]
#        strain = pokes.strains[i]
#        
#        # prepare b in Linear problem Ax=b.
#        c .= 0
#        c[site] = strain
#        mul!(b,H,c)
#        linsolve = LinearSolve.set_b(sol.cache,b)
#        sol = solve(linsolve)
#        energies[i] = sol ⋅ c 
#    end
#    return energies
#end



function computeEnergiesFast(Q::Network, pokes::Pokes)
    # compute the energies of applied strains
    
    # make sure that all the sites are the same. Make more general when needed.
    @assert all(Ref(pokes.sites[1]) .== pokes.sites)
    n = length(pokes.sites)
    @assert n>=1
    numDims = size(Q.H,1)

    # initial and build arrays
    H, B, b, c = Q.H, Q.B, Q.b, Q.c 
    fillBMatrix!(B,H, pokes.sites[1]) 
    energies = Vector{Float64}(undef,n)

    # first iteracion
    c .= 0
    c[pokes.sites[1]] = pokes.strains[1]
    mul!(b,H,c)
    prob = LinearProblem(B,b)
    linsolve = init(prob, alias_A=true, alias_b=true)
    sol = solve(linsolve)
    energies[1] = sol ⋅ c 
    
    
    for i in 2:n
        site = pokes.sites[i]
        strain = pokes.strains[i]
        
        # prepare b in Linear problem Ax=b.
        c .= 0
        c[site] = strain
        mul!(b,H,c)
        linsolve = LinearSolve.set_b(sol.cache,b)
        sol = solve(linsolve)
        energies[i] = sol ⋅ c 
    end
    return energies
end



function fillBMatrix!(B::FloatMatrix, H::FloatMatrix, site::IntVector)
    B .= -1 .* H
    for j in site, i in axes(B,1)
        B[i,j] = ifelse(i==j, 1.0, 0.0)
    end
    return nothing
end


function getTranslationRotation(r::FloatVector)
    
    # r is a vector of positions r = (x1,y1,x2,y2,...)
    N = length(r)
    @assert iseven(N)
    n = Int(N/2)
    
    # move origin to center of mass
    xy = r2xy(r)
    x_mean, y_mean = mean(xy, dims=1)
    xy[:,1] .-= x_mean
    xy[:,2] .-= y_mean
    r = xy2r(xy)

    # compute x and y translations vector
    vx = normalize([ones(n) zeros(n)]'[:])
    vy = normalize([zeros(n) ones(n)]'[:])
    
    # compute rotation vector
    rotMat = [0 -1; 1 0]
    R = zeros(N, N)
    for i in 1:n
        R[2i-1:2i,2i-1:2i] = rotMat
    end
    vr = normalize(R * r)
    return [vx vy vr]
end


function removeTranslationRotation!(dr::FloatVector, r::FloatVector)
    # dr0 : applied displacement
    # r : postion vector

    V = getTranslationRotation(r)
    dr .-= V * V' * dr
    normalize!(dr)
    return nothing
end


function createNormalStrains(numStrains::Int, r::FloatVector, actSiteNodes::IntVector)
    
    sites = sort([actSiteNodes * 2; actSiteNodes * 2 .- 1])
    numZeroModes = 3
    numDims = length(sites)
    if numDims < (numZeroModes + numStrains)
        error("Error: Not enough degrees of freedom in the strain to make $numStrains strains")
    end
    
    r_act = r[sites]
    V = getTranslationRotation(r_act)
    strains = Matrix{Float64}(undef,numDims, numStrains)

    for i in 1:numStrains
        strain = randn(numDims)
        strain .-= V * V' * strain
        normalize!(strain)
        strains[:,i] .= strain
        V = [V strain]
    end
    return [strains[:,i] for i in axes(strains,2)]
end





