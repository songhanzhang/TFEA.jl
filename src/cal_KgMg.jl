function cal_KgMg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])

    println("\n")
    println("*** Kg and Mg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Kg = spzeros(n_DOF,n_DOF)
    Mg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/10))) == 0
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
            A = Reals[i_real,2][1]
            Izz = Real[i_real,2][4]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            Le = sqrt((xj-xi)^2 + (yj-yi)^2)
            Ke_bar = cal_Ke_bar_2DEulerBeam(Le,E,A,Izz)
            Te = cal_Te_2DBeam(xi,yi,xj,yj)
            Ke = transpose(Te) * Ke_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(i_node + 0.6),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_5 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOF_6 = findall(isequal(j_node + 0.6),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4;DOF_5;DOF_6]
            Kg[DOFs,DOFs] += Ke
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
        elseif Elements[i_e,2] == "2D_QuadraticTriangle"
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
            PlaneType = "PlaneStrain"
            (Ke,Me) = cal_KeMe_QuadraticTriangular(x,y,E,ν,ρ,PlaneType)
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
        elseif Elements[i_e,2] == "2D_Acoustics_4N"
            node_1 = Elements[i_e,5][1]
            node_2 = Elements[i_e,5][2]
            node_3 = Elements[i_e,5][3]
            node_4 = Elements[i_e,5][4]
            x1 = Nodes[node_1,2]
            y1 = Nodes[node_1,3]
            x2 = Nodes[node_2,2]
            y2 = Nodes[node_2,3]
            x3 = Nodes[node_3,2]
            y3 = Nodes[node_3,3]
            x4 = Nodes[node_4,2]
            y4 = Nodes[node_4,3]
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            c_air = Materials[i_mat,2][1]
            t = Reals[i_real,2][1]
            

            Gauss = [ -1/(sqrt(3))  -1/(sqrt(3))  1
                       1/(sqrt(3))  -1/(sqrt(3))  1
                       1/(sqrt(3))   1/(sqrt(3))  1
                      -1/(sqrt(3))   1/(sqrt(3))  1 ]
            Ke = zeros(4,4)
            for i_Gauss = 1:size(Gauss,1)
                ξ = Gauss[i_Gauss,1]
                η = Gauss[i_Gauss,2]
                H = Gauss[i_Gauss,3]
                N1 = 1/4*(1-ξ)*(1-η)
                N2 = 1/4*(1+ξ)*(1-η)
                N3 = 1/4*(1+ξ)*(1+η)
                N4 = 1/4*(1-ξ)*(1+η)
                dN1_dξ = -1/4*(1-η)
                dN2_dξ =  1/4*(1-η)
                dN3_dξ =  1/4*(1+η)
                dN4_dξ = -1/4*(1+η)
                dN1_dη = -1/4*(1-ξ)
                dN2_dη = -1/4*(1+ξ)
                dN3_dη =  1/4*(1+ξ)
                dN4_dη =  1/4*(1-ξ)
                dx_dξ = dN1_dξ*x1 + dN2_dξ*x2 + dN3_dξ*x3 + dN4_dξ*x4
                dy_dξ = dN1_dξ*y1 + dN2_dξ*y2 + dN3_dξ*y3 + dN4_dξ*y4
                dx_dη = dN1_dη*x1 + dN2_dη*x2 + dN3_dη*x3 + dN4_dη*x4
                dy_dη = dN1_dη*y1 + dN2_dη*y2 + dN3_dη*y3 + dN4_dη*y4
                J = [ dx_dξ  dy_dξ
                      dx_dη  dy_dη ]
                Temp = J\[ dN1_dξ; dN1_dη ]
                dN1_dx = Temp[1]
                dN1_dy = Temp[2]
                Temp = J\[ dN2_dξ; dN2_dη ]
                dN2_dx = Temp[1]
                dN2_dy = Temp[2]
                Temp = J\[ dN3_dξ; dN3_dη ]
                dN3_dx = Temp[1]
                dN3_dy = Temp[2]
                Temp = J\[ dN4_dξ; dN4_dη ]
                dN4_dx = Temp[1]
                dN4_dy = Temp[2]
                dN_dx = [ dN1_dx  dN2_dx  dN3_dx  dN4_dx ]
                dN_dy = [ dN1_dy  dN2_dy  dN3_dy  dN4_dy ]
                Ke += c_air^2*t*(transpose(dN_dx)*dN_dx + transpose(dN_dy)*dN_dy)*abs(det(J))*H
            end
            DOF_1 = findall(isequal(node_1 + 0.7),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(node_2 + 0.7),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(node_3 + 0.7),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(node_4 + 0.7),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Kg[DOFs,DOFs] += Ke
        end
    end

    return Kg, Mg

end
