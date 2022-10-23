module Elastic_N5

push!(LOAD_PATH,"./")
using LinearAlgebra, PyPlot, Elastic_N3, Elastic_N4

export Network, computeHessian!, computeHessian, plotDisplacment!, buildNetworkAndTable


###############################################
struct Network
    seq::Vector{Int}
    r::Vector{Float64}
    S::Matrix{Float64}
    H::Matrix{Float64}
end

##############################################

function computeHessian!(H,r, S)
    # compute the hessian matrix of the network
    N = size(S,1)
    H .= 0
    for j in 1:N # index over nodes, not dimentions
        sxx, syy, sxy = 0.0, 0.0, 0.0
        j2 = 2j
        xj = r[j2-1]
        yj = r[j2]
        for i in 1:N # index over nodes, not dimentions
            i2 = 2i
            k = S[i,j]
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


function computeHessian(r, S)
    N = size(S,1)
    H = Matrix{Float64}(undef, 2N, 2N)
    computeHessian!(H,r,S)
    return H
end

function plotDisplacment!(ax, r, dr; zorder=1,
    width=0.004, color="r", angles="xy", scale_units="xy", scale=0.5)
    xy = r2xy(r)
    x,y = xy[:,1], xy[:,2]
    uv = r2xy(dr)
    u,v = uv[:,1], uv[:,2]
    ax.quiver(x,y,u,v, zorder=zorder, width = width, color=color,
    scale=scale, angles=angles, scale_units=scale_units)
    return nothing
end

function buildNetworkAndTable(W, L, q, disorder, k_min, k_max)
    r, A = buildStructure(W, L, disorder)
    seq = rand(1:q, W*L)
    K = buildSpringTable(q, k_min, k_max)
    S = seq2Springs(seq, K, A);
    H = computeHessian(r, S);
    net = Network(seq, r, S, H)
    return net, K
end

end