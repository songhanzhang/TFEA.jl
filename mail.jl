# push!(LOAD_PATH,"/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src")
# import TFEA

include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")

using MAT
using LinearAlgebra
using SparseArrays

file = matopen("/Users/songhan.zhang/Documents/MATLAB/2023-QuadraticTriangle/model.mat")
Nodes_import = read(file, "Nodes")
n_nodes = size(Nodes,1)
Elements_import = read(file, "Elements")
n_elements = size(Elements_import,1)

Nodes = [ 1:1:n_nodes  Nodes_import ]
Elements = Array{Any,2}(undef,n_elements,5)
for i_e = 1:n_elements
    Elements[i_e,1] = i_e
    Elements[i_e,2] = "2D_QuadraticTriangle"
    Elements[i_e,3] = 1
    Elements[i_e,4] = 1
    Elements[i_e,5] = Elements_import[i_e,:]
end
Materials = [ 1  (2e11, 7850, 0.33) ]
Reals = [ 1  (1) ]

(n_DOF, list_DOF) = cal_list_DOF(n_nodes, [1,2])

(Kg, Mg) = cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])
