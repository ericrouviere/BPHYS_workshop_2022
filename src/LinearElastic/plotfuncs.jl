function plotNetwork!(ax, r::Vector, adjacencyMatrix::Matrix;
                      bondColor::Vector=[0.7,0.7,0.7],
                      bondWidth=3,
                      nodeSize=30)
    # Plots network, modifies figure axis object in place
    # uses PyPlot, LinearAlgebra
    xy = r2xy(r)
    x = xy[:,1]
    y = xy[:,2]
    N = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)
    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            if A[i,j] != 0
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
            end
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    return nothing
end


function plotNetwork(r::FloatVector, A::FloatMatrix)
    fig, ax = subplots()
    plotNetwork!(ax, r, A)
    return fig
end

plotNetwork!(ax, Q::Network) = plotNetwork!(ax, Q.r, Q.SM)

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


function plotBondStiffness!(ax, r::Vector,
                            adjacencyMatrix::NumMatrix,
                            SM::Matrix;
                            bondWidth=3,
                            nodeSize=30,
                            cmap="Purples",
                            k_min=1e-2,
                            k_max=1e1)
    # Plots network, modifies figure axis object in place
    # uses PyPlot, LinearAlgebra
    cmap = get_cmap(cmap)
    k_min = log(10,k_min)
    k_max = log(10,k_max)

    xy = r2xy(r)
    x = xy[:,1]
    y = xy[:,2]
    N = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)
    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            if A[i,j] != 0
                bondColor = cmap( (log(10,SM[i,j]) - k_min) / (k_max - k_min))
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
            end
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    return nothing
end

function plotBondStiffness!(ax,
                            Q::Network;
                            bondWidth=3,
                            nodeSize=30,
                            cmap="Purples",
                            k_min=1e-2,
                            k_max=1e1)
    plotBondStiffness!(ax, Q.r, Q.SM, Q.SM)
    return nothing
end





