include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")
work_path = "/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/"
using MAT
using LinearAlgebra
using SparseArrays
using Plots
using Printf
using Measures

Nodes = zeros(49,3)
counter = 0
locs = [ 0, 1.5, 3, 5, 7, 8.5, 10 ]*1e-3
for ii = 1:7
    for jj = 1:7
        counter += 1
        Nodes[counter,1] = counter
        Nodes[counter,2] = locs[ii]
        Nodes[counter,3] = locs[jj]
    end
end
n_nodes = size(Nodes,1)

Elements = [
    1   "2D_QuadTriangle"  1  1  (1,3,15,2,9,8)
    2   "2D_QuadTriangle"  1  1  (1+2,3+2,15+2,2+2,9+2,8+2)
    3   "2D_QuadTriangle"  1  1  (1+4,3+4,15+4,2+4,9+4,8+4)
    4   "2D_QuadTriangle"  1  1  (1+14,3+14,15+14,2+14,9+14,8+14)
    5   "2D_QuadTriangle"  2  1  (1+16,3+16,15+16,2+16,9+16,8+16)
    6   "2D_QuadTriangle"  1  1  (1+18,3+18,15+18,2+18,9+18,8+18)
    7   "2D_QuadTriangle"  1  1  (1+28,3+28,15+28,2+28,9+28,8+28)
    8   "2D_QuadTriangle"  1  1  (1+30,3+30,15+30,2+30,9+30,8+30)
    9   "2D_QuadTriangle"  1  1  (1+32,3+32,15+32,2+32,9+32,8+32)
    10  "2D_QuadTriangle"  1  1  (3,17,15,10,16,9)
    11  "2D_QuadTriangle"  1  1  (3+2,17+2,15+2,10+2,16+2,9+2)
    12  "2D_QuadTriangle"  1  1  (3+4,17+4,15+4,10+4,16+4,9+4)
    13  "2D_QuadTriangle"  1  1  (3+14,17+14,15+14,10+14,16+14,9+14)
    14  "2D_QuadTriangle"  2  1  (3+16,17+16,15+16,10+16,16+16,9+16)
    15  "2D_QuadTriangle"  1  1  (3+18,17+18,15+18,10+18,16+18,9+18)
    16  "2D_QuadTriangle"  1  1  (3+28,17+28,15+28,10+28,16+28,9+28)
    17  "2D_QuadTriangle"  1  1  (3+30,17+30,15+30,10+30,16+30,9+30)
    18  "2D_QuadTriangle"  1  1  (3+32,17+32,15+32,10+32,16+32,9+32)
]
n_elements = size(Elements,1)

Materials = [ 1  (2e9,   1000, 0.3)
              2  (200e9, 8000, 0.3) ]

Reals = [ 1  (1) ]

fig_model = plot(
    size = (600,600),
    dpi = 900,
    legend = false,
    grid = false,
    tickfontsize = 10,
    frame_style = :box
)
plot_model_elements(Nodes,Elements)
plot_model_nodes(Nodes)
for i_n = 1:n_nodes
    if size(Nodes,2) == 3
        xi = Nodes[i_n,2]
        yi = Nodes[i_n,3]
        annotate!([xi], [yi],
                  Plots.text(string(i_n), 10, :dodgerblue, :bottom))
    else
        xi = Nodes[i_n,2]
        yi = Nodes[i_n,3]
        zi = Nodes[i_n,4]
        annotate!([xi], [yi], [zi],
                  Plots.text(string(i_n), 10, :dodgerblue, :bottom))
    end
end
for i_e = 1:n_elements
    node_1 = Elements[i_e,5][1]
    node_2 = Elements[i_e,5][2]
    node_3 = Elements[i_e,5][3]
    x1 = Nodes[node_1,2]
    y1 = Nodes[node_1,3]
    x2 = Nodes[node_2,2]
    y2 = Nodes[node_2,3]
    x3 = Nodes[node_3,2]
    y3 = Nodes[node_3,3]
    annotate!([(x1+x2+x3)/3], [(y1+y2+y3)/3],
              Plots.text(string("  ",i_e), 10, :gray70, :bottom))
end
savefig(string(work_path,"fig_model_meta.pdf"))

(n_DOF, list_DOF) = cal_list_DOF(n_nodes, [1,2])

(Kg,Mg,Cg) = cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF)

k_ax = zeros(3000,2)
for ii = 1:1000
    k_ax[ii,1] = ii/1000*pi/0.01
    k_ax[ii,2] = 0
end
for ii = 1001:2000
    k_ax[ii,1] = pi/0.01
    k_ax[ii,2] = (ii-1000)/1000*pi/0.01
end
for ii = 2001:3000
    k_ax[ii,1] = (3001-ii)/1000*pi/0.01
    k_ax[ii,2] = (3001-ii)/1000*pi/0.01
end

ω_save = zeros(n_DOF-13*2,3000)*0im
for ii = 1:3000

    kx = k_ax[ii,1]
    ky = k_ax[ii,2]

    P = zeros(n_DOF,n_DOF)*0im
    for ii = 1:n_DOF
        P[ii,ii] = 1
    end
    DOF_master = [3,4,5,6,7,8,9,10,11,12,15,16,29,30,43,44,57,58,71,72,1,2,1,2,1,2]
    DOF_slave  = [87,88,89,90,91,92,93,94,95,96,27,28,41,42,55,56,69,70,83,84,85,86,97,98,13,14]
    relation = zeros(26)*0im
    relation[1:10] = ones(10)*exp(1im*kx*0.01)
    relation[11:20] = ones(10)*exp(1im*ky*0.01)
    relation[21] = exp(1im*kx*0.01)
    relation[22] = exp(1im*kx*0.01)
    relation[23] = exp(1im*(kx+ky)*0.01)
    relation[24] = exp(1im*(kx+ky)*0.01)
    relation[25] = exp(1im*ky*0.01)
    relation[26] = exp(1im*ky*0.01)

    for jj = 1:26
        P[DOF_slave[jj],DOF_master[jj]] = relation[jj]
    end

    P = P[:,setdiff(1:n_DOF,DOF_slave)]

    Kg_bc = P'*Kg*P
    Mg_bc = P'*Mg*P

    (λ,Φ) = eigen(Kg_bc,Mg_bc)
    ω_save[:,ii] = sqrt.(Complex.(λ))
end
ω_save = real(ω_save)
for ii = 1:3000
    ω_save[:,ii] = sort(ω_save[:,ii])
end

fig_kw = plot(
    size = (600,500),
    dpi = 900,
    legend = false,
    grid = false,
    frame_style = :box,
    tickfontsize = 10
)
color_ax = [:dodgerblue,:green3,:firebrick,:blue,:pink,:gray,:hotpink1,:orange]
for ii = 1:8
    scatter!(
        ω_save[ii,:]/2/pi/1e3, label = "",
        color = color_ax[ii], markerstrokecolor = color_ax[ii], markersize = 1
)
end
plot!(xticks = ([0,1000,2000,3000],["O","Γ","M","K"]))
plot!(Shape([0,3000,3000,0],
      [64,64,77,77]), label = "",
      color = :gray80, linecolor = :gray80, opacity=.9, linewidth = 0)
xlabel!("Wavevector")
ylabel!("Frequency (kHz)")
plot!(ylims = (0,140))
