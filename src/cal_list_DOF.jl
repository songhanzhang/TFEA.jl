function cal_list_DOF(n_nodes::Int64, n_dir_list)

    n_DOF = length(n_dir_list) * n_nodes
    list_DOF = zeros(n_DOF,2)
    counter = 0
    for i_node = 1:n_nodes
        for i_dir = 1:length(n_dir_list)
            counter += 1
            list_DOF[counter,1] = counter
            list_DOF[counter,2] = i_node + 0.1*n_dir_list[i_dir]
        end
    end

    return n_DOF, list_DOF

end
