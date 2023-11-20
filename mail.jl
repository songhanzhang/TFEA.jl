# push!(LOAD_PATH,"/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src")
# import TFEA

include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")

using MAT

file = matopen("/Users/songhan.zhang/Documents/MATLAB/2023-QuadraticTriangle/model.mat")
Nodes = read(file, "Nodes")
n_nodes = size(Nodes,1)
Elements_node = read(file, "Elements")
n_elements = size(Elements_node,1)

Nodes = [ 1:1:n_nodes  Nodes ]
Elements = Array{Any,2}(undef,n_elements,5)
for i_e = 1:n_elements
    Elements[i_e,1] = i_e
    Elements[i_e,2] = "2D_QuadraticTriangle"
    Elements[i_e,3] = 1
    Elements[i_e,4] = 1
    Elements[i_e,5] = Elements_node[i_e,:]
end
