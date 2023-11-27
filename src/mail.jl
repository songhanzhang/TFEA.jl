# push!(LOAD_PATH,"/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src")
# import TFEA

include("/Users/songhan.zhang/Documents/GitHub/TFEA.jl/src/TFEA.jl")

using MAT
using LinearAlgebra
using SparseArrays
using Plots
using Printf
using Measures

file = matopen("/Users/songhan.zhang/Documents/MATLAB/2023-QuadraticTriangle/model.mat")
Nodes_import = read(file, "Nodes")
n_nodes = size(Nodes_import,1)
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

plot_model = plot(
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
savefig("/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/model.pdf")

# %% Excitation source
# node 137 (0.0, 0.3)
Time_label = 0:1e-7:1.0e-3
Ns = length(Time_label)
ΔT = Time_label[2] - Time_label[1]
ctf = 100e3
n_peaks = 5
n_steps = Int(floor(n_peaks/ctf/ΔT))
tar_DOF = Int(findall(isequal(137.1),list_DOF[:,2])[1])
theta = 0.0
Fg = zeros(n_DOF,length(Time_label))
T_0 = 11
for ii = T_0 : T_0 + n_steps
    Fg[tar_DOF,ii] = (100*(1-cos(2*pi*ctf*(ii-T_0)*ΔT/n_peaks))*sin(2*pi*ctf*(ii-T_0)*ΔT))*cos(theta)
end

fig_F = plot(size = (560,300),
              legend = false,
              grid = false,
              frame_style = :box,
              tickfontsize = 10)
plot!(fig_F, Time_label*1e6, Fg[tar_DOF,:], label = "",
      color = :black, w = 1, linestyle = :solid)
xlabel!(fig_F, "Time (μs)", guidefontsize = 10)
ylabel!(fig_F, "Force (N)")
savefig("/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/excitation.pdf")

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

plot(Time_label, Ug[(7792-1)*2+1,:])
plot!(Time_label, Ug[(137-1)*2+1,:])
u_sel_1 = Ug[(137-1)*2+1,:]*1e9
u_sel_2 = Ug[(7792-1)*2+1,:]*1e9
t_ax = Time_label[:,:]

fig_u_sel = plot(
    size = (600,300),
    dpi = 300,
    legend = false,
    grid = false,
    frame_style = :box,
    tickfontsize = 10
)
plot!(t_ax*1e6, u_sel_1, w = 1.5, color = :pink1)
plot!(t_ax*1e6, u_sel_2, w = 1.5, color = :dodgerblue)
xlabel!("Time (μs)")
ylabel!("Displacement (μm)")
plot!(bottom_margin = 3mm)
savefig("/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/ut.pdf")

# %% Plot wave field
ani = @animate for i_t = 20:20:2000
# for i_t = 50:50:10000
    println("i_t = ", i_t)
    u_xy = transpose(reshape(Ug[:,i_t],2,n_nodes))

    x_ax = 0:0.005:1
    y_ax = 0:0.005:0.6

    ux_mat = zeros(length(x_ax),length(y_ax))

    for (i_x,x) in enumerate(x_ax)

        if mod(i_x, Int(floor(length(x_ax)/10))) == 0
            println("Progress: ", i_x/length(x_ax)*100, " %")
        end

        for (i_y,y) in enumerate(y_ax)
            dist = zeros(4,5)
            for ii = 1:4
                for jj = 1:5
                    x_central = 0.2 + (ii-1)*0.1
                    y_central = 0.1 + (jj-1)*0.1
                    dist[ii,jj] = norm([x-x_central, y-y_central])
                end
            end
            if minimum(dist) < 0.02
                ux_mat[i_x,i_y] = NaN
                continue
            end
            for i_e = 1:n_elements
                node_1 = Int(Elements[i_e,5][1])
                node_2 = Int(Elements[i_e,5][2])
                node_3 = Int(Elements[i_e,5][3])
                x1 = Nodes[node_1,2]
                y1 = Nodes[node_1,3]
                x2 = Nodes[node_2,2]
                y2 = Nodes[node_2,3]
                x3 = Nodes[node_3,2]
                y3 = Nodes[node_3,3]
                if judge_point_inside_triangle(x1,y1,x2,y2,x3,y3,x,y)
                    St = ( x2*y3 + x3*y1 + x1*y2 - x2*y1 - x1*y3 - x3*y2 )/2
                    S1 = ( x3*yp + xp*y2 + x2*y3 - x3*y2 - x2*yp - xp*y3 )/2
                    S2 = ( x1*yp + xp*y3 + x3*y1 - x1*y3 - x3*yp - xp*y1 )/2
                    ξ = S2/St
                    η = S3/St
                    Nb = zeros(6)
                    Nb[1] = (1-ξ-η)*(1-2*ξ-2*η)
                    Nb[2] = ξ*(2*ξ-1)
                    Nb[3] = η*(2*η-1)
                    Nb[4] = 4*ξ*(1-ξ-η)
                    Nb[5] = 4*ξ*η
                    Nb[6] = 4*η*(1-ξ-η)
                    DOF_1  = Int(findall(isequal(node_1+0.1),list_DOF[:,2])[1])
                    DOF_3  = Int(findall(isequal(node_2+0.1),list_DOF[:,2])[1])
                    DOF_5  = Int(findall(isequal(node_3+0.1),list_DOF[:,2])[1])
                    DOF_7  = Int(findall(isequal(node_4+0.1),list_DOF[:,2])[1])
                    DOF_9  = Int(findall(isequal(node_5+0.1),list_DOF[:,2])[1])
                    DOF_11 = Int(findall(isequal(node_6+0.1),list_DOF[:,2])[1])
                    DOFs = [DOF_1;DOF_3;DOF_5;DOF_7;DOF_9;DOF_11]
                    Ue = Ug[DOFs,i_t]
                    ux_mat[i_x,i_y] = transpose(Nb)*Ue
                    continue
                end
            end
        end
    end
    fig_ux = plot(size = (600,350),
                  dpi = 300,
                  legend = false,
                  grid = false,
                  frame_style = :box,
                  tickfontsize = 10,
                  aspect_ratio = :equal)
    heatmap!(transpose(ux_mat), c = :jet, aspect_ratio = :equal, clim = (-5e-3,5e-3))
    title_content = @sprintf "t = %0.2f μs" Time_label[i_t]*1e6
    title!(title_content, titlefont = 10)
end
gif(
    ani,
    "/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/ani.gif",
    fps=20
)


i_t = 2300
u_xy = transpose(reshape(Ug[:,i_t],2,n_nodes))

x_ax = 0:0.002:1
y_ax = 0:0.002:0.6

t_ax = Time_label[10:10:10000]
ux_mat = zeros(length(x_ax),length(y_ax),length(t_ax))

for (i_x,x) in enumerate(x_ax)

    if mod(i_x, Int(floor(length(x_ax)/100))) == 0
        println("Progress: ", i_x/length(x_ax)*100, " %")
    end

    for (i_y,y) in enumerate(y_ax)
        dist = zeros(4,5)
        for ii = 1:4
            for jj = 1:5
                x_central = 0.2 + (ii-1)*0.1
                y_central = 0.1 + (jj-1)*0.1
                dist[ii,jj] = norm([x-x_central, y-y_central])
            end
        end
        if minimum(dist) < 0.02
            ux_mat[i_x,i_y,:] *= NaN
            continue
        end
        for i_e = 1:n_elements
            node_1 = Int(Elements[i_e,5][1])
            node_2 = Int(Elements[i_e,5][2])
            node_3 = Int(Elements[i_e,5][3])
            node_4 = Int(Elements[i_e,5][4])
            node_5 = Int(Elements[i_e,5][5])
            node_6 = Int(Elements[i_e,5][6])
            x1 = Nodes[node_1,2]
            y1 = Nodes[node_1,3]
            x2 = Nodes[node_2,2]
            y2 = Nodes[node_2,3]
            x3 = Nodes[node_3,2]
            y3 = Nodes[node_3,3]
            if judge_point_inside_triangle(x1,y1,x2,y2,x3,y3,x,y)
                St = ( x2*y3 + x3*y1 + x1*y2 - x2*y1 - x1*y3 - x3*y2 )/2
                # S1 = ( x3*y  +  x*y2 + x2*y3 - x3*y2 - x2*y  -  x*y3 )/2
                S2 = ( x1*y  +  x*y3 + x3*y1 - x1*y3 - x3*y  -  x*y1 )/2
                S3 = ( x2*y  +  x*y1 + x1*y2 - x2*y1 - x1*y  -  x*y2 )/2
                ξ = S2/St
                η = S3/St
                Nb = zeros(6)
                Nb[1] = (1-ξ-η)*(1-2*ξ-2*η)
                Nb[2] = ξ*(2*ξ-1)
                Nb[3] = η*(2*η-1)
                Nb[4] = 4*ξ*(1-ξ-η)
                Nb[5] = 4*ξ*η
                Nb[6] = 4*η*(1-ξ-η)
                DOF_1 = (node_1-1)*2 + 1
                DOF_2 = (node_2-1)*2 + 1
                DOF_3 = (node_3-1)*2 + 1
                DOF_4 = (node_4-1)*2 + 1
                DOF_5 = (node_5-1)*2 + 1
                DOF_6 = (node_6-1)*2 + 1
                DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6]
                Ue = Ug[DOFs,10:10:10000]
                ux_mat[i_x,i_y,:] = transpose(Nb)*Ue
                continue
            end
        end
    end
end

ani = @animate for i_t = 1:1:100
    fig_ux = plot(size = (600,350),
                  dpi = 300,
                  legend = false,
                  grid = false,
                  frame_style = :box,
                  tickfontsize = 10,
                  aspect_ratio = :equal)
    heatmap!(transpose(ux_mat[:,:,i_t]),
             c = :balance, aspect_ratio = :equal, clim = (-2e-10,2e-10))
    title_content = @sprintf "t = %0.2f μs" t_ax[i_t]*1e6
    title!(title_content, titlefont = 10)
end
gif(
    ani,
    "/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/ani.gif",
    fps=20
)

fig_ux = plot(size = (600,380),
              dpi = 600,
              legend = false,
              grid = false,
              frame_style = :box,
              tickfontsize = 10,
              aspect_ratio = :equal)
heatmap!(fig_ux, x_ax, y_ax, transpose(ux_mat),
         c = :balance, aspect_ratio = :equal,
         clim = (-2e-10,2e-10),
         xlims = (0,1), ylims = (0,0.6))
title!("t = 230 μs", titlefontsize = 10)
xlabel!("x (m)")
ylabel!("y (m)")
savefig("/Users/songhan.zhang/Documents/Julia/2023-TFEA-v1120-AcMetaMat/wave_field_230.png")

# contour(transpose(ux_mat),
#         level = 1600, linewidth = 0, fillrange = true, c = :jet, aspect_ratio = :equal, clim = (-1e-3,1e-3))
