function plot_model_elements(Nodes, Elements)
    for i_e = 1:n_elements
        if Elements[i_e,2] == "2D_QuadraticTriangle"
            node_1 = Elements[i_e,5][1]
            node_2 = Elements[i_e,5][2]
            node_3 = Elements[i_e,5][3]
            x1 = Nodes[node_1,2]
            y1 = Nodes[node_1,3]
            x2 = Nodes[node_2,2]
            y2 = Nodes[node_2,3]
            x3 = Nodes[node_3,2]
            y3 = Nodes[node_3,3]
            plot!(Shape([x1,x2,x3], [y1,y2,y3]), color = :pink, linewidth = 1)
        end
    end
end
