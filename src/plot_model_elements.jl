function plot_model_elements(Nodes, Elements)
    for i_e = 1:n_elements
        if Elements[i_e,2] == "2D_QuadTriangle"
            node_1 = Int(Elements[i_e,5][1])
            node_2 = Int(Elements[i_e,5][2])
            node_3 = Int(Elements[i_e,5][3])
            node_4 = Int(Elements[i_e,5][4])
            node_5 = Int(Elements[i_e,5][5])
            node_6 = Int(Elements[i_e,5][6])
            x = zeros(6)
            y = zeros(6)
            x[1] = Nodes[node_1,2]
            y[1] = Nodes[node_1,3]
            x[2] = Nodes[node_2,2]
            y[2] = Nodes[node_2,3]
            x[3] = Nodes[node_3,2]
            y[3] = Nodes[node_3,3]
            x[4] = Nodes[node_4,2]
            y[4] = Nodes[node_4,3]
            x[5] = Nodes[node_5,2]
            y[5] = Nodes[node_5,3]
            x[6] = Nodes[node_6,2]
            y[6] = Nodes[node_6,3]
            coor_interp = [ 0.0  0.0
                            0.1  0.0
                            0.2  0.0
                            0.3  0.0
                            0.4  0.0
                            0.5  0.0
                            0.6  0.0
                            0.7  0.0
                            0.8  0.0
                            0.9  0.0
                            1.0  0.0
                            0.9  0.1
                            0.8  0.2
                            0.7  0.3
                            0.6  0.4
                            0.5  0.5
                            0.4  0.6
                            0.3  0.7
                            0.2  0.8
                            0.1  0.9
                            0.0  1.0
                            0.0  0.9
                            0.0  0.8
                            0.0  0.7
                            0.0  0.6
                            0.0  0.5
                            0.0  0.4
                            0.0  0.3
                            0.0  0.2
                            0.0  0.1 ]
            n_interp = size(coor_interp,1)
            x_profile = zeros(n_interp)
            y_profile = zeros(n_interp)
            for i_interp = 1:n_interp
                ξ = coor_interp[i_interp,1]
                η = coor_interp[i_interp,2]
                Nb = zeros(6)
                Nb[1] = (1-ξ-η) * (1-2*ξ-2*η)
                Nb[2] = ξ * (2*ξ-1)
                Nb[3] = η * (2*η-1)
                Nb[4] = 4 * ξ * (1-ξ-η)
                Nb[5] = 4 * ξ * η
                Nb[6] = 4 * η * (1-ξ-η)
                x_profile[i_interp] = transpose(Nb)*x
                y_profile[i_interp] = transpose(Nb)*y
            end
            plot!(Shape(x_profile,y_profile), color = :gray90, linewidth = 1)
        end
    end
end
