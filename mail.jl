# push!(LOAD_PATH,"/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src")
# import TFEA

include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")

using MAT
using LinearAlgebra
using SparseArrays
using Plots
using Printf

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
Materials = [ 1  (2e11, 7850, 0.3) ]
Reals = [ 1  (1) ]

(n_DOF, list_DOF) = cal_list_DOF(n_nodes, [1,2])

(Kg, Mg) = cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])
Cg = spzeros(n_DOF,n_DOF)

plot_mode = plot(
    size = (560,400),
    dpi = 600,
    grid = false,
    legend = false,
    frame_style = :box,
    aspect_ratio = :equal,
    tickfontsize = 10
)
for i_e = 1:n_elements
    node_1 = Int(Elements[i_e,5][1])
    node_2 = Int(Elements[i_e,5][2])
    node_3 = Int(Elements[i_e,5][3])
    x = zeros(4)
    y = zeros(4)
    x[1] = Nodes[node_1,2]
    y[1] = Nodes[node_1,3]
    x[2] = Nodes[node_2,2]
    y[2] = Nodes[node_2,3]
    x[3] = Nodes[node_3,2]
    y[3] = Nodes[node_3,3]
    x[4] = Nodes[node_1,2]
    y[4] = Nodes[node_1,3]
    plot!(x,y,label = "",
          w = 1.0, color = :dodgerblue, linestyle = :solid)
end
plot!()

# %% Excitation source
# node 137 (0.0, 0.3)
Time_label = 0:1e-7:1.0e-3
Ns = length(Time_label)
ΔT = Time_label[2] - Time_label[1]
ctf = 20e3
n_peaks = 5
n_steps = Int(floor(n_peaks/ctf/ΔT))
tar_DOF = Int(findall(isequal(137.1),list_DOF[:,2])[1])
theta = 0.0
Fg = zeros(n_DOF,length(Time_label))
T_0 = 11
for ii = T_0 : T_0 + n_steps
    Fg[tar_DOF,ii] = (100*(1-cos(2*pi*ctf*(ii-T_0)*ΔT/n_peaks))*sin(2*pi*ctf*(ii-T_0)*ΔT))*cos(theta)
end

fig_F = plot(size = (560,200),
              legend = false,
              grid = false,
              frame_style = :box,
              tickfontsize = 10)
plot!(fig_F, Time_label*1e6, Fg[tar_DOF_y,:], label = "",
      color = :black, w = 1, linestyle = :solid)
xlabel!(fig_F, "Time (μs)", guidefontsize = 10)
ylabel!(fig_F, "Force (N)")

# %% Solve - direct time integral
Ug = zeros(n_DOF,length(Time_label))
Ug_0 = zeros(n_DOF,1)
Ug_t_0 = zeros(n_DOF,1)
Ug_tt_0 = zeros(n_DOF,1)
for i_t = 3:length(Time_label)
    if mod(i_t,10) == 0
        @printf "%.3f%% completed ... \n" i_t/Ns*100
    end
        K_eq = Mg/ΔT^2 + Cg/(2*ΔT)
        F_eq = Fg[:,i_t-1] + (2*Mg/ΔT^2-Kg)*Ug[:,i_t-1] - (Mg/ΔT^2-Cg/(2*ΔT))*Ug[:,i_t-2]
        Ug[:,i_t] = K_eq\F_eq
end
