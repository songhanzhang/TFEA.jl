function gen_Nodes_grid(nnx, nny, Δx, Δy)

    n_nodes = nnx * nny
    Nodes = zeros(n_nodes,4)

    for i_y = 1:nny
        for i_x = 1:nnx
            counter = (i_y-1)*nnx + i_x
            Nodes[counter,1] = counter
            Nodes[counter,2] = (i_x-1)*Δx
            Nodes[counter,3] = (i_y-1)*Δy
        end
    end

    return Nodes

end
