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
n_nodes = size(Nodes,1)

Elements = [
    1  "2D_Euler_Beam"  1  1  (1,2)
    2  "2D_Euler_Beam"  1  1  (2,3)
    3  "2D_Euler_Beam"  1  1  (4,2)
    4  "2D_Euler_Beam"  1  1  (2,5)
]
n_elements = size(Elements,1)

Materials = [ 1  (2e11, 7850, 0.3) ]

Reals = [ 1  (1, 1, 1, 1/64*pi*0.01^4) ]

(n_DOF, list_DOF) = cal_list_DOF(n_nodes, [1,2,6])

(Kg,Mg,Cg) = cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF)

kx = 2
ky = 0
P = zeros(15,15)*0im
for ii = 1:15
    P[ii,ii] = 1
end
DOF_master = 1
DOF_slave  = 7
relation = exp(1im*kx*2)
P = P[:,setdiff(1:n_DOF,DOF_slave)]
P[DOF_slave,DOF_master] = relation

Kg_bc = P'*Kg*P
Mg_bc = P'*Mg*P

(λ,Φ) = eigen(Kg_bc,Mg_bc)

k_ax = zeros(3000,2)
for ii = 1:1000
    k_ax[ii,1] = ii/1000*2*pi
    k_ax[ii,2] = 0
end
for ii = 1001:2000
    k_ax[ii,1] = 2*pi
    k_ax[ii,2] = (ii-1000)/1000*2*pi
end
for ii = 3000:-1:2001
    k_ax[ii,1] = (ii-2000)/1000*2*pi
    k_ax[ii,2] = (ii-2000)/1000*2*pi
end

ω_save = zeros(9,3000)*0im
for ii = 1:3000

    kx = k_ax[ii,1]
    ky = k_ax[ii,2]

    P = zeros(15,15)*0im
    for ii = 1:15
        P[ii,ii] = 1
    end
    DOF_master = [1,2,3,10,11,12]
    DOF_slave  = [7,8,9,13,14,15]
    relation = [exp(1im*kx*2),exp(1im*kx*2),exp(1im*kx*2),exp(1im*ky*2),exp(1im*ky*2),exp(1im*ky*2)]
    for jj = 1:6
        P[DOF_slave[jj],DOF_master[jj]] = relation[jj]
    end

    P = P[:,setdiff(1:n_DOF,DOF_slave)]

    Kg_bc = P'*Kg*P
    Mg_bc = P'*Mg*P

    (λ,Φ) = eigen(Kg_bc,Mg_bc)
    ω_save[:,ii] = sqrt.(Complex.(λ))
end
ω_save = real(ω_save)

fig_kw = plot(
    size = (600,1000),
    dpi = 900,
    legend = false,
    grid = false,
    frame_style = :box,
    tickfontsize = 10
)
for ii = 1:9
    scatter!(
        ω_save[ii,:]/1e3, label = "",
        color = :dodgerblue, markerstrokecolor = :dodgerblue, markersize = 1
)
end
plot!(left_margin = 6mm)
plot!(xticks = ([0,1000,2000,3000],["O","Γ","M","K"]))
plot!(Shape([0,3000,3000,0],
      [6.5,6.5,8.7,8.7]), label = "",
      color = :gray80, linecolor = :gray80, opacity=.9, linewidth = 0)
xlabel!("Wavevector")
ylabel!("ω (kHz)")
P = zeros(15,15)*0im
for ii = 1:15
    P[ii,ii] = 1
end
DOF_master = [1,2,3,10,11,12]
DOF_slave  = [7,8,9,13,14,15]
relation = [exp(1im*kx*2),exp(1im*kx*2),exp(1im*kx*2),exp(1im*ky*2),exp(1im*ky*2),exp(1im*ky*2)]
P = P[:,setdiff(1:n_DOF,DOF_slave)]
for ii = 1:6
    P[DOF_slave[ii],DOF_master[ii]] = relation[ii]
end
P = P[:,setdiff(1:n_DOF,DOF_slave)]
