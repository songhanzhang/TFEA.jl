function cal_Ke_2D_LGL_36(E,ν,t,Nodes_xy)

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

    Ke = zeros(72,72)
    for i_Gauss = 1:size(Gauss,1)
      ξ  = Gauss[i_Gauss,1]
      η  = Gauss[i_Gauss,2]
      H  = abs(Gauss[i_Gauss,3])
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
      D = E/(1-v^2) * [ 1  v  0
                        v  1  0
                        0  0  (1-v)/2 ]
      B = zeros(3,72)
      for ii = 1:36
          Ni_xi  = N_xi[ii]
          Ni_eta = N_eta[ii]
          Ni_xy = J\[Ni_xi; Ni_eta]
          B[1,(ii-1)*2+1] = Ni_xy[1];
          B[2,(ii-1)*2+2] = Ni_xy[2];
          B[3,(ii-1)*2+1] = Ni_xy[2];
          B[3,(ii-1)*2+2] = Ni_xy[1];
      end
      Ke = Ke + t*transpose(B)*D*B*abs(det(J))*H
    end

    return Ke

end
