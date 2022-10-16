function cal_Kg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])

    println("\n")
    println("*** Kg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Kg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/10))) == 0
            println("    Progress ", i_e/n_elements*100, " %")
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
        elseif Element[i_e,2] == "3D_Euler_Beam"
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
        elseif Element[i_e,2] == "3D_Tmsk_Beam"
            
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
        end
    end

    return Kg

end
