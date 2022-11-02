module Elastic_N3

export drawLattice, buildAdjacency, xy2r,  r2xy, buildStructure, plotNetwork!

function drawLattice(W,L)
    u1 = [cos(pi/6), sin(pi/6)]
    u2 = [0, 1] 
    xy = zeros(W*L, 2)
    k,c = 0,0
    for i in 1:L
        for j in 1:W
            k += 1
            b = j + c
            xy[k,:] = i * u1 .+ b * u2
        end
        if i % 2 ==0
            c -= 1
        end 
    end
    return xy
end


function buildAdjacency(xy)
    numNodes = size(xy, 1)
    A = zeros(numNodes, numNodes) # initialize a matrix
    for i in 1:numNodes,j in 1:numNodes
        # compute the distance between nodes i and j
        if i==j; continue; end
        dx = xy[i,1] - xy[j,1]
        dy = xy[i,2] - xy[j,2]
        dist = sqrt(dx^2 + dy^2)
        if dist < 1.1 # lattice spacing is 1.
            A[i,j] = 1.0
        end
    end
    return A
end


xy2r(xy) = xy'[:]
r2xy(r) = Array(reshape(r, 2, :)')

function buildStructure(W,L,disorder)
    # build the structure and adjacency matrix of  
    xy = drawLattice(W,L)
    A = buildAdjacency(xy)
    xy .+= disorder .* randn(size(xy))
    r = xy2r(xy)
    return r, A  
end

function plotNetwork!(ax, # PyPlot axes object
                      r, # structure vector
                      A; # Adjaceny matrix or Spring matrix
                      bondColor=[0.7,0.7,0.7],
                      bondWidth=3,
                      nodeSize=30)
                      
    # Plots network, modifies figure axis object in place
    xy = r2xy(r)
    x = xy[:,1]
    y = xy[:,2]
    N = size(A,1)
    for i =2:N, j in 1:(i-1) 
        if A[i,j] != 0
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    ax.axis("off")
    return nothing
end


end
