function cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF;
                  Nodes_a = [], pml_interface = [], model_boundary = [], eta_max = 0)

    println("\n")
    println("*** Kg and Mg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Kg = spzeros(n_DOF,n_DOF)
    Mg = spzeros(n_DOF,n_DOF)
    Cg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/1))) == 0
            println("Progress ", i_e/n_elements*100, " %")
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
            Le = sqrt((xj-xi)^2 + (yj-yi)^2)
            Ke_bar = cal_Ke_bar_2DBar(Le,E,A)
            Te = cal_Te_2DBar(xi,yi,xj,yj)
            Ke = transpose(Te) * Ke_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Kg[DOFs,DOFs] += Ke
        elseif Elements[i_e,2] == "2D_Euler_Beam"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            A = Reals[i_real,2][1]
            Izz = Reals[i_real,2][4]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            Le = sqrt((xj-xi)^2 + (yj-yi)^2)
            (Ke_bar,Me_bar) = cal_KeMe_bar_2DEulerBeam(Le,E,ρ,A,Izz)
            Te = cal_Te_2DBeam(xi,yi,xj,yj)
            Ke = transpose(Te) * Ke_bar * Te
            Me = transpose(Te) * Me_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(i_node + 0.6),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_5 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOF_6 = findall(isequal(j_node + 0.6),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6]
            Kg[DOFs,DOFs] += Ke
            Mg[DOFs,DOFs] += Me
        elseif Elements[i_e,2] == "3D_Euler_Beam"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            ν = Materials[i_mat,2][3]
            G = E/(2*(1+ν))
            A = Reals[i_real,2][1]
            Ixx = Reals[i_real,2][2]
            Iyy = Reals[i_real,2][3]
            Izz = Reals[i_real,2][4]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            k_node = Elements[i_e,5][3]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            zi = Nodes[i_node,4]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            zj = Nodes[j_node,4]
            xk = Nodes[k_node,2]
            yk = Nodes[k_node,3]
            zk = Nodes[k_node,4]
            Le = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            Ke_bar = cal_Ke_bar_3DEulerBeam(Le,E,G,A,Ixx,Iyy,Izz)
            Te = cal_Te_3DBeam(xi,yi,zi,xj,yj,zj,xk,yk,zk)
            Ke = transpose(Te) * Ke_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(i_node + 0.6),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_5 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOF_6 = findall(isequal(j_node + 0.6),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6]
            Kg[DOFs,DOFs] += Ke
        elseif Elements[i_e,2] == "3D_Tmsk_Beam"

        elseif Elements[i_e,2] == "2D_LGL_36n"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            v = Materials[i_mat,2][3]
            t = Reals[i_real,2][1]
            Nodes_xy = zeros(36,2)
            for ii = 1:36
                ii_node = Elements[i_e,5][ii]
                Nodes_xy[ii,1:2] = [Nodes[ii_node,2] Nodes[ii_node,3]]
            end
            Ke = cal_Ke_2D_LGL_36(E,ν,t,Nodes_xy)
            DOFs = Array{Int64,1}(undef,72)
            for i_node = 1:36
                for i_dir = 1:2
                    n_node = Elements[i_e,5][i_node]
                    DOFs[(i_node-1)*2+i_dir] = Int(findall(isequal(n_node+0.1*i_dir),list_DOF[:,2])[1])
                end
            end
            Kg[DOFs,DOFs] += Ke
        elseif Elements[i_e,2] == "2D_QuadTriangle"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            ν = Materials[i_mat,2][3]
            t = Reals[i_real,2][1]
            node_1 = Int(Elements[i_e,5][1])
            node_2 = Int(Elements[i_e,5][2])
            node_3 = Int(Elements[i_e,5][3])
            node_4 = Int(Elements[i_e,5][4])
            node_5 = Int(Elements[i_e,5][5])
            node_6 = Int(Elements[i_e,5][6])
            x = zeros(6)
            y = zeros(6)
            x[1] = Nodes[node_1,2]
            y[1] = Nodes[node_1,3]
            x[2] = Nodes[node_2,2]
            y[2] = Nodes[node_2,3]
            x[3] = Nodes[node_3,2]
            y[3] = Nodes[node_3,3]
            x[4] = Nodes[node_4,2]
            y[4] = Nodes[node_4,3]
            x[5] = Nodes[node_5,2]
            y[5] = Nodes[node_5,3]
            x[6] = Nodes[node_6,2]
            y[6] = Nodes[node_6,3]
            (Ke,Me,Ce) = cal_KeMe_QuadTriangle(x,y,E,ν,ρ,"PlaneStrain",pml_interface,model_boundary,eta_max)
            DOF_1  = Int(findall(isequal(node_1+0.1),list_DOF[:,2])[1])
            DOF_2  = Int(findall(isequal(node_1+0.2),list_DOF[:,2])[1])
            DOF_3  = Int(findall(isequal(node_2+0.1),list_DOF[:,2])[1])
            DOF_4  = Int(findall(isequal(node_2+0.2),list_DOF[:,2])[1])
            DOF_5  = Int(findall(isequal(node_3+0.1),list_DOF[:,2])[1])
            DOF_6  = Int(findall(isequal(node_3+0.2),list_DOF[:,2])[1])
            DOF_7  = Int(findall(isequal(node_4+0.1),list_DOF[:,2])[1])
            DOF_8  = Int(findall(isequal(node_4+0.2),list_DOF[:,2])[1])
            DOF_9  = Int(findall(isequal(node_5+0.1),list_DOF[:,2])[1])
            DOF_10 = Int(findall(isequal(node_5+0.2),list_DOF[:,2])[1])
            DOF_11 = Int(findall(isequal(node_6+0.1),list_DOF[:,2])[1])
            DOF_12 = Int(findall(isequal(node_6+0.2),list_DOF[:,2])[1])
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6;DOF_7;DOF_8;DOF_9;DOF_10;DOF_11;DOF_12]
            Kg[DOFs,DOFs] += Ke
            Mg[DOFs,DOFs] += Me
            Cg[DOFs,DOFs] += Ce
        elseif Elements[i_e,2] == "2D_Acoustics_4N"
            node_1 = Elements[i_e,5][1]
            node_2 = Elements[i_e,5][2]
            node_3 = Elements[i_e,5][3]
            node_4 = Elements[i_e,5][4]
            x = zeros(4)
            y = zeros(4)
            x[1] = Nodes[node_1,2]
            y[1] = Nodes[node_1,3]
            x[2] = Nodes[node_2,2]
            y[2] = Nodes[node_2,3]
            x[3] = Nodes[node_3,2]
            y[3] = Nodes[node_3,3]
            x[4] = Nodes[node_4,2]
            y[4] = Nodes[node_4,3]
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            c_air = Materials[i_mat,2][1]
            t = Reals[i_real,2][1]
            (Ke,Me) = cal_KeMe_2D_Acoustics_4N(x,y,c_air,t)
            DOF_1 = findall(isequal(node_1 + 0.7),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(node_2 + 0.7),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(node_3 + 0.7),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(node_4 + 0.7),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Kg[DOFs,DOFs] += Ke
            Mg[DOFs,DOFs] += Me
        elseif Elements[i_e,2] == "2D_QuadTriangle_ac"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            c_air = Materials[i_mat,2][1]
            t = Reals[i_real,2][1]
            node_1 = Int(Elements[i_e,5][1])
            node_2 = Int(Elements[i_e,5][2])
            node_3 = Int(Elements[i_e,5][3])
            node_4 = Int(Elements[i_e,5][4])
            node_5 = Int(Elements[i_e,5][5])
            node_6 = Int(Elements[i_e,5][6])
            x = zeros(6)
            y = zeros(6)
            x[1] = Nodes[node_1,2]
            y[1] = Nodes[node_1,3]
            x[2] = Nodes[node_2,2]
            y[2] = Nodes[node_2,3]
            x[3] = Nodes[node_3,2]
            y[3] = Nodes[node_3,3]
            x[4] = Nodes[node_4,2]
            y[4] = Nodes[node_4,3]
            x[5] = Nodes[node_5,2]
            y[5] = Nodes[node_5,3]
            x[6] = Nodes[node_6,2]
            y[6] = Nodes[node_6,3]
            (Ke,Me,Ce) = cal_KeMe_QuadTriangle_ac(x,y,c_air,t,pml_interface,model_boundary,eta_max)
            DOF_1  = Int(findall(isequal(node_1+0.1),list_DOF[:,2])[1])
            DOF_2  = Int(findall(isequal(node_1+0.2),list_DOF[:,2])[1])
            DOF_3  = Int(findall(isequal(node_2+0.1),list_DOF[:,2])[1])
            DOF_4  = Int(findall(isequal(node_2+0.2),list_DOF[:,2])[1])
            DOF_5  = Int(findall(isequal(node_3+0.1),list_DOF[:,2])[1])
            DOF_6  = Int(findall(isequal(node_3+0.2),list_DOF[:,2])[1])
            DOF_7  = Int(findall(isequal(node_4+0.1),list_DOF[:,2])[1])
            DOF_8  = Int(findall(isequal(node_4+0.2),list_DOF[:,2])[1])
            DOF_9  = Int(findall(isequal(node_5+0.1),list_DOF[:,2])[1])
            DOF_10 = Int(findall(isequal(node_5+0.2),list_DOF[:,2])[1])
            DOF_11 = Int(findall(isequal(node_6+0.1),list_DOF[:,2])[1])
            DOF_12 = Int(findall(isequal(node_6+0.2),list_DOF[:,2])[1])
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6;DOF_7;DOF_8;DOF_9;DOF_10;DOF_11;DOF_12]
            Kg[DOFs,DOFs] += Ke
            Mg[DOFs,DOFs] += Me
            Cg[DOFs,DOFs] += Ce
        elseif Elements[i_e,2] == "3D_Hexahedral"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            ν = Materials[i_mat,2][3]
            Nodes_xy = zeros(8,3)
            for ii = 1:8
                ii_node = Elements[i_e,5][ii]
                Node_xy[ii,1:3] = [ Nodes[ii_node,2]  Nodes[ii_node,3]  Nodes[ii_node,4] ] 
            end
            tp = 1/sqrt(3)
            Gauss = [ -tp -tp -tp
                      -tp -tp  tp
                      -tp  tp -tp
                      -tp  tp  tp
                       tp -tp -tp
                       tp -tp  tp
                       tp  tp -tp
                       tp  tp  tp ]
            Me = zeros(24,24)
            Ke = zeros(24,24)
            for i_Gauss = 1:size(Gauss,1)
                ξ   = Gauss(i_Gauss,1)
                η  = Gauss(i_Gauss,2)
                ζ = Gauss(i_Gauss,3)
                N = zeros(1,8)
                N[1,1] = 1/8*(1-ξ)*(1-η)*(1+ζ)
                N[1,2] = 1/8*(1+ξ)*(1-η)*(1+ζ)
                N[1,3] = 1/8*(1+ξ)*(1-η)*(1-ζ)
                N[1,4] = 1/8*(1-ξ)*(1-η)*(1-ζ)
                N[1,5] = 1/8*(1-ξ)*(1+η)*(1+ζ)
                N[1,6] = 1/8*(1+ξ)*(1+η)*(1+ζ)
                N[1,7] = 1/8*(1+ξ)*(1+η)*(1-ζ)
                N[1,8] = 1/8*(1-ξ)*(1+η)*(1-ζ)
                dN_dξ = zeros(1,8)
                dN_dξ[1,1] = -1/8*(1-η)*(1+ζ)
                dN_dξ[1,2] =  1/8*(1-η)*(1+ζ)
                dN_dξ[1,3] =  1/8*(1-η)*(1-ζ)
                dN_dξ[1,4] = -1/8*(1-η)*(1-ζ)
                dN_dξ[1,5] = -1/8*(1+η)*(1+ζ)
                dN_dξ[1,6] =  1/8*(1+η)*(1+ζ)
                dN_dξ[1,7] =  1/8*(1+η)*(1-ζ)
                dN_dξ[1,8] = -1/8*(1+η)*(1-ζ)
                dN_dη = zeros(1,8)
                dN_dη[1,1] = -1/8*(1-ξ)*(1+ζ)
                dN_dη[1,2] = -1/8*(1+ξ)*(1+ζ)
                dN_dη[1,3] = -1/8*(1+ξ)*(1-ζ)
                dN_dη[1,4] = -1/8*(1-ξ)*(1-ζ)
                dN_dη[1,5] =  1/8*(1-ξ)*(1+ζ)
                dN_dη[1,6] =  1/8*(1+ξ)*(1+ζ)
                dN_dη[1,7] =  1/8*(1+ξ)*(1-ζ)
                dN_dη[1,8] =  1/8*(1-ξ)*(1-ζ)
                dN_dζ = zeros(1,8)
                dN_dζ[1,1] =  1/8*(1-ξ)*(1-η)
                dN_dζ[1,2] =  1/8*(1+ξ)*(1-η)
                dN_dζ[1,3] = -1/8*(1+ξ)*(1-η)
                dN_dζ[1,4] = -1/8*(1-ξ)*(1-η)
                dN_dζ[1,5] =  1/8*(1-ξ)*(1+η)
                dN_dζ[1,6] =  1/8*(1+ξ)*(1+η)
                dN_dζ[1,7] = -1/8*(1+ξ)*(1+η)
                dN_dζ[1,8] = -1/8*(1-ξ)*(1+η)
                dx_dξ = dN_dξ * Nodes_xy[:,1]
                dx_dη = dN_dη * Nodes_xy[:,1]
                dx_dζ = dN_dζ * Nodes_xy[:,1]
                dy_dξ = dN_dξ * Nodes_xy[:,2]
                dy_dη = dN_dη * Nodes_xy[:,2]
                dy_dζ = dN_dζ * Nodes_xy[:,2]
                dz_dξ = dN_dξ * Nodes_xy[:,3]
                dz_dη = dN_dη * Nodes_xy[:,3]
                dz_dζ = dN_dζ * Nodes_xy[:,3]
                J = [ dx_dξ  dy_dξ  dz_dξ
                      dx_dη  dy_dη  dz_dη
                      dx_dζ  dy_dζ  dz_dζ ]
                B = zeros(6,24)
                for ii = 1:8
                    Temp = J\[dN_dξ[1,ii]; dN_dη[1,ii]; dN_dζ[1,ii]]
                    dNi_dx = Temp[1,1]
                    dNi_dy = Temp[2,1]
                    dNi_dz = Temp[3,1]
                    B[1,(ii-1)*3+1] = dNi_dx
                    B[2,(ii-1)*3+2] = dNi_dy
                    B[3,(ii-1)*3+3] = dNi_dz
                    B[4,(ii-1)*3+1] = dNi_dy
                    B[4,(ii-1)*3+2] = dNi_dx
                    B[5,(ii-1)*3+1] = dNi_dz
                    B[5,(ii-1)*3+3] = dNi_dx                    
                    B[6,(ii-1)*3+2] = dNi_dz
                    B[6,(ii-1)*3+3] = dNi_dy
                end
                λ = E*ν/((1+ν)*(1-2*ν))
                μ = E/(2*(1+ν))
                D = [ λ+2*μ λ     λ     0 0 0
                      λ     λ+2*μ λ     0 0 0
                      λ     λ     λ+2*μ 0 0 0
                      0     0     0     μ 0 0
                      0     0     0     0 μ 0
                      0     0     0     0 0 μ ]
                NN = kron(N, [1 0 0; 0 1 0; 0 0 1])
                Me = Me + ρ*transpose(NN)*NN*abs(det(J))*H
                Ke = Ke + transpose(B)*D*B*abs(det(J))*H                
            end
            DOF_1  = Int(findall(isequal(node_1+0.1),list_DOF[:,2])[1])
            DOF_2  = Int(findall(isequal(node_1+0.2),list_DOF[:,2])[1])
            DOF_3  = Int(findall(isequal(node_1+0.3),list_DOF[:,2])[1])
            DOF_4  = Int(findall(isequal(node_2+0.1),list_DOF[:,2])[1])
            DOF_5  = Int(findall(isequal(node_2+0.2),list_DOF[:,2])[1])
            DOF_6  = Int(findall(isequal(node_2+0.3),list_DOF[:,2])[1])
            DOF_7  = Int(findall(isequal(node_3+0.1),list_DOF[:,2])[1])
            DOF_8  = Int(findall(isequal(node_3+0.2),list_DOF[:,2])[1])
            DOF_9  = Int(findall(isequal(node_3+0.3),list_DOF[:,2])[1])
            DOF_10 = Int(findall(isequal(node_4+0.1),list_DOF[:,2])[1])
            DOF_11 = Int(findall(isequal(node_4+0.2),list_DOF[:,2])[1])
            DOF_12 = Int(findall(isequal(node_4+0.3),list_DOF[:,2])[1])
            DOF_13 = Int(findall(isequal(node_5+0.1),list_DOF[:,2])[1])
            DOF_14 = Int(findall(isequal(node_5+0.2),list_DOF[:,2])[1])
            DOF_15 = Int(findall(isequal(node_5+0.3),list_DOF[:,2])[1])
            DOF_16 = Int(findall(isequal(node_6+0.1),list_DOF[:,2])[1])
            DOF_17 = Int(findall(isequal(node_6+0.2),list_DOF[:,2])[1])
            DOF_18 = Int(findall(isequal(node_6+0.3),list_DOF[:,2])[1])
            DOF_19 = Int(findall(isequal(node_7+0.1),list_DOF[:,2])[1])
            DOF_20 = Int(findall(isequal(node_7+0.2),list_DOF[:,2])[1])
            DOF_21 = Int(findall(isequal(node_7+0.3),list_DOF[:,2])[1])
            DOF_22 = Int(findall(isequal(node_8+0.1),list_DOF[:,2])[1])
            DOF_23 = Int(findall(isequal(node_8+0.2),list_DOF[:,2])[1])
            DOF_24 = Int(findall(isequal(node_8+0.3),list_DOF[:,2])[1])
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6;DOF_7;DOF_8;DOF_9;DOF_10;DOF_11;DOF_12;DOF_13;DOF_14;DOF_15;DOF_16;DOF_17;DOF_18;DOF_19;DOF_20;DOF_21;DOF_22;DOF_23;DOF_24]
            Kg[DOFs,DOFs] += Ke
            Mg[DOFs,DOFs] += Me
        end
    end

    return Kg, Mg, Cg

end
