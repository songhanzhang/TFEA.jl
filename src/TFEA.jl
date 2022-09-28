module TFEA

include("gen_Nodes_grid.jl")

function cal_list_DOF(n_nodes::Int64, n_dir)
    n_DOF = length(n_dir) * n_nodes
    list_DOF = zeros(n_DOF,2)
    counter = 0
    for i_node = 1:n_nodes
        for i_dir = 1:length(n_dir)
            counter += 1
            list_DOF[counter,1] = counter
            list_DOF[counter,2] = i_node + 0.1*n_dir[i_dir]
        end
    end
    return n_DOF, list_DOF
end

function cal_Kg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])

    println("\n")
    println("*** Kg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Kg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/10))) == 0
            println("    Progress ", i_e/n_elements*100, " %")
        end

        if Elements[i_e,2] == "2D_Bar"

            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            A = Reals[i_real,2][1]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            Le = sqrt((xj-xi)^2+(yj-yi)^2)
            Ke_bar = E*A/Le*[  1   0  -1   0
                               0   0   0   0
                              -1   0   1   0
                               0   0   0   0]
            ss = (yj-yi)/Le
            cs = (xj-xi)/Le
            Te = [ cs   ss    0    0
                  -ss   cs    0    0
                    0    0   cs   ss
                    0    0  -ss   cs ]
            Ke = transpose(Te) * Ke_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Kg[DOFs,DOFs] += Ke
        end
    end
    return Kg
end

function cal_Mg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])

    println("\n")
    println("*** Mg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Mg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/10))) == 0
            println("    Progress ", i_e/n_elements*100, " %")
        end

        if Elements[i_e,2] == "2D_Bar"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            A = Reals[i_real,2][1]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            Le = sqrt((xj-xi)^2+(yj-yi)^2)
            Me_bar = ρ*A/Le*[ 1/3   0  1/6   0
                                0   0    0   0
                              1/6   0  1/3   0
                                0   0    0   0 ]
            ss = (yj-yi)/Le
            cs = (xj-xi)/Le
            Te = [ cs   ss    0    0
                  -ss   cs    0    0
                    0    0   cs   ss
                    0    0  -ss   cs ]
            Me = transpose(Te) * Me_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Mg[DOFs,DOFs] += Me
        end
    end
    return Mg
end

end
