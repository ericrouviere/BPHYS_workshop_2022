module Elastic_N6

push!(LOAD_PATH,"./")
using LinearAlgebra, Statistics, PyPlot, Elastic_N3, Elastic_N4, Elastic_N5

export computeResponse, removeTranslationRotation!

function computeResponse(r, H, strain, site)
    # Solve the linear response problem.

    # build vector (Δr_a, 0)
    c = zeros(length(r))
    c[site] = strain

    # multiply 
    b = H * c

    # build B matrix
    B = -H
    B[:,site] = I(size(B,1))[:,site]

    # solve linear problem
    x = B \ b 

    # get Force vector
    F = copy(c)
    F[site] = x[site]

    # get displacement vector
    Δr = copy(x)
    Δr[site] = c[site]

    # compute Energy
    E = dot(F, Δr)

    return E, F, Δr 
end

function getTranslationRotation(r)
    
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

function removeTranslationRotation!(dr, r)
    # dr : applied displacement
    # r : postion vector

    V = getTranslationRotation(r)
    dr .-= V * V' * dr
    #normalize!(dr)
    return nothing
end


end