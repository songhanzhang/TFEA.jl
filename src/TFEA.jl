module TFEA

# Write your package code here.
function TFEA_test()
    println("Test successful!")
end

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

end
