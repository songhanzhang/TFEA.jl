function cal_KeMe_Hexa(x,y,z,E,ρ,ν,pml_interface,model_boundary,eta_max)
    
    Nodes_xy = [ x y z ]

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
    Ce = zeros(24,24)
    Ke = zeros(24,24)
    for i_Gauss = 1:size(Gauss,1)
        ξ = Gauss[i_Gauss,1]
        η = Gauss[i_Gauss,2]
        ζ = Gauss[i_Gauss,3]
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
        Me = Me + ρ*transpose(NN)*NN*abs(det(J))
        Ke = Ke + transpose(B)*D*B*abs(det(J))
        if !isempty(pml_interface)
            x_gauss = transpose(Nb)*x
            y_gauss = transpose(Nb)*y
            # eta_pml = cal_pml_eta(pml_interface, model_boundary, x_gauss, y_gauss, eta_max)
            x_Lb = pml_interface[1]
            x_Rb = pml_interface[2]
            y_Bb = pml_interface[3]
            y_Tb = pml_interface[4]

            x_min = model_boundary[1]
            x_max = model_boundary[2]
            y_min = model_boundary[3]
            y_max = model_boundary[4]

            if x_gauss < x_Lb
                eta_pml = eta_max*((x_Lb-x_gauss)/(x_Lb-x_min))^3
            elseif x_gauss > x_Rb
                eta_pml = eta_max*((x_gauss-x_Rb)/(x_max-x_Rb))^3
            elseif y_gauss < y_Bb
                eta_pml = eta_max*((y_Bb-y_gauss)/(y_Bb-y_min))^3
            elseif y_gauss > y_Tb
                eta_pml = eta_max*((y_gauss-y_Tb)/(y_max-y_Tb))^3
            else
                eta_pml = 0.0
            end
            Ce += Me * eta_pml
        end
    end

    return Ke,Me,Ce

end