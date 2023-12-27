include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")

using MAT
using LinearAlgebra
using SparseArrays
using Plots
using Printf
using Measures

Nodes = [
    1  -1   0   0
    2   0   0   0
    3   1   0   0
    4   0  -1   0
    5   0   1   0
]

Elements = [
    1  "2D_Euler_Beam"  1  1  (1,2)
    2  "2D_Euler_Beam"  1  1  (2,3)
    3  "2D_Euler_Beam"  1  1  (4,2)
    4  "2D_Euler_Beam"  1  1  (2,5)
]

Materials = [ 1  (2e11, 7850, 0.3) ]

Reals = [ 1  (1) ]
