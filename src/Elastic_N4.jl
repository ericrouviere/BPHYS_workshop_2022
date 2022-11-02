module Elastic_N4

push!(LOAD_PATH,"./")
using LinearAlgebra, PyPlot, Elastic_N3

export buildSpringTable, seq2Springs, plotBondStiffness!

function buildSpringTable(q, k_min, k_max)
    log_diff = log(10, k_max) - log(10, k_min)
    K =  10 .^ (rand(q, q) .* log_diff .+ log(10,k_min))
    K = Float64.(Symmetric(K, :L))
end

function seq2Springs(seq, K, A)
    # convert a sequence into a matrix that stores
    # the spring constant between node i and node j.
    S = zeros(size(A))
    for j in 1:size(A,2), i in 1:size(A,1)
        S[i,j] = A[i,j] * K[ seq[i], seq[j]]
    end
    return S
end

function plotBondStiffness!(ax,r, A, S;
                            bondWidth=3,
                            nodeSize=30,
                            cmap="Purples",
                            k_min=1e-2,
                            k_max=1e1)
    # Plots network coloring springs by their stiffness
    # uses PyPlot, LinearAlgebra
    cmap = get_cmap(cmap)
    k_min = log(10,k_min)
    k_max = log(10,k_max)
    xy = r2xy(r)
    x = xy[:,1]
    y = xy[:,2]
    N = size(A,1)
    for i =2:N, j in 1:(i-1)
        if A[i,j] != 0
            bondColor = cmap( (log(10,S[i,j]) - k_min) / (k_max - k_min))
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    ax.axis("off")
    return nothing
end



end
