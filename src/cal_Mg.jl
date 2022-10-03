function cal_Mg(Nodes, Elements, Materials, Reals, list_DOF; Nodes_a = [])

    println("\n")
    println("*** Mg evaluation started ...\n")

    n_DOF = size(list_DOF,1)
    n_elements = size(Elements,1)
    Mg = spzeros(n_DOF,n_DOF)

    for i_e = 1:n_elements

        if mod(i_e, Int(floor(n_elements/10))) == 0
            println("    Progress ", i_e/n_elements*100, " %")
        end

        if Elements[i_e,2] == "2D_Bar"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            E = Materials[i_mat,2][1]
            ρ = Materials[i_mat,2][2]
            A = Reals[i_real,2][1]
            i_node = Elements[i_e,5][1]
            j_node = Elements[i_e,5][2]
            xi = Nodes[i_node,2]
            yi = Nodes[i_node,3]
            xj = Nodes[j_node,2]
            yj = Nodes[j_node,3]
            e_nodes = [ xi yi
                        xj yj ]
            Le = sqrt((xj-xi)^2 + (yj-yi)^2)
            Me_bar = cal_Me_bar_2DBar(Le,ρ,A)
            Te = cal_Te_2DBar(e_nodes)
            Me = transpose(Te) * Me_bar * Te
            DOF_1 = findall(isequal(i_node + 0.1),list_DOF[:,2])[1]
            DOF_2 = findall(isequal(i_node + 0.2),list_DOF[:,2])[1]
            DOF_3 = findall(isequal(j_node + 0.1),list_DOF[:,2])[1]
            DOF_4 = findall(isequal(j_node + 0.2),list_DOF[:,2])[1]
            DOFs = [DOF_1;DOF_2;DOF_3;DOF_4]
            Mg[DOFs,DOFs] += Me
        elseif Elements[i_e,2] == "2D_LGL_36n"
            i_mat = Elements[i_e,3]
            i_real = Elements[i_e,4]
            ρ = Materials[i_mat,2][2]
            t = Reals[i_real,2][1]
            Nodes_xy = zeros(36,2)
            for ii = 1:36
                ii_node = Elements[i_e,5][ii]
                Nodes_xy[ii,1:2] = [Nodes[ii_node,2] Nodes[ii_node,3]]
            end
            ξ_label = zeros(6,1)
            ξ_label[1] = -1
            ξ_label[2] = -sqrt(1/3+2/(3*sqrt(7)))
            ξ_label[3] = -sqrt(1/3-2/(3*sqrt(7)))
            ξ_label[4] =  sqrt(1/3-2/(3*sqrt(7)))
            ξ_label[5] =  sqrt(1/3+2/(3*sqrt(7)))
            ξ_label[6] =  1
            η_label = ξ_label[:]
            Gauss = zeros(36,3)
            for (i_ξ, ξ) in enumerate(ξ_label)
                for (i_η, η) in enumerate(η_label)
                    Gauss[(i_ξ-1)*6+i_η,1] = ξ
                    Gauss[(i_ξ-1)*6+i_η,2] = η
                    Gauss[(i_ξ-1)*6+i_η,3] = 1/15/(1/8*(63*ξ^5-70*ξ^3+15*ξ))^2 * 1/15/(1/8*(63*η^5-70*η^3+15*η))^2
                end
            end
            Me = zeros(72,72)
            for i_Gauss = 1:size(Gauss,1)
                ξ  = Gauss[i_Gauss,1]
                η  = Gauss[i_Gauss,2]
                H  = abs(Gauss[i_Gauss,3])
                N = ones(1,36)
                for aa = 1:6
                    for bb = 1:6
                        for ii = 1:6
                            if ii != aa
                                N[1, (aa-1)*6+bb] *= (ξ-ξ_label[ii])/(ξ_label[aa]-ξ_label[ii])
                            end
                        end
                        for jj = 1:6
                            if jj != bb
                                N[1, (aa-1)*6+bb] *= (η-η_label[jj])/(η_label[bb]-η_label[jj])
                            end
                        end
                    end
                end
                N_xi = zeros(1,36)
                for aa = 1:6
                    for bb = 1:6
                        for kk = 1:6
                            if kk != aa
                                Temp = 1
                                for ii = 1:6
                                    if ii != aa && ii != kk
                                        Temp *= (ξ-ξ_label[ii])/(ξ_label[aa]-ξ_label[ii])
                                    elseif ii == kk
                                        Temp *= 1/(ξ_label[aa] - ξ_label[kk])
                                    end
                                end
                                N_xi[1, (aa-1)*6+bb] += Temp
                            end
                        end
                        for jj = 1:6
                            if jj != bb
                                N_xi[1, (aa-1)*6+bb] *= (η-η_label[jj])/(η_label[bb]-η_label[jj])
                            end
                        end
                    end
                end
                N_eta = zeros(1,36)
                for aa = 1:6
                    for bb = 1:6
                        for kk = 1:6
                            if kk != bb
                                Temp = 1
                                for ii = 1:6
                                    if ii != bb && ii != kk
                                        Temp *= (η-η_label[ii])/(η_label[bb]-η_label[ii])
                                    elseif ii == kk
                                        Temp *= 1/(η_label[bb] - η_label[kk])
                                    end
                                end
                                N_eta[1, (aa-1)*6+bb] += Temp
                            end
                        end
                        for jj = 1:6
                            if jj != aa
                                N_eta[1, (aa-1)*6+bb] *= (ξ-ξ_label[jj])/(ξ_label[aa]-ξ_label[jj])
                            end
                        end
                    end
                end
                x_xi  = N_xi  * Nodes_xy[:,1]
                x_eta = N_eta * Nodes_xy[:,1]
                y_xi  = N_xi  * Nodes_xy[:,2]
                y_eta = N_eta * Nodes_xy[:,2]
                J = [ x_xi   y_xi
                      x_eta  y_eta ]
                NN = kron(N, [1 0; 0 1])
                Me = Me + ρ*t*transpose(NN)*NN*abs(det(J))*H
            end
            DOFs = Array{Int64,1}(undef,72)
            for i_node = 1:36
                for i_dir = 1:2
                    n_node = Elements[i_e,5][i_node]
                    DOFs[(i_node-1)*2+i_dir] = Int(findall(isequal(n_node+0.1*i_dir),list_DOF[:,2])[1])
                end
            end
            for ii = 1:72
                for jj = 1:72
                    if ii != jj
                        Me[ii,jj] = 0
                    end
                end
            end
            Mg[DOFs,DOFs] = Mg[DOFs,DOFs] + Me
        end
    end

    return Mg

end
