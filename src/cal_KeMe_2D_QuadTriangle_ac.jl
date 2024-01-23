function cal_KeMe_QuadTriangle_ac(x,y,c_air,t;pml_interface = [],model_boundary = [], eta_max = 0)

    Gauss = [ 0.0915762135  0.8168475730  0.1099517437
              0.0915762135  0.0915762135  0.1099517437
              0.8168475730  0.0915762135  0.1099517437
              0.4459484909  0.1081030182  0.2233815897
              0.4459484909  0.4459484909  0.2233815897
              0.1081030182  0.4459484909  0.2233815897 ];

    Ke = zeros(12,12)
    Me = zeros(12,12)
    Ce = zeros(12,12)

    for i_Gauss = 1:6
        ξ = Gauss[i_Gauss,1]
        η = Gauss[i_Gauss,2]
        H = Gauss[i_Gauss,3]
        Nb = zeros(6)
        Nb[1] = (1-ξ-η) * (1-2*ξ-2*η)
        Nb[2] = ξ * (2*ξ-1)
        Nb[3] = η * (2*η-1)
        Nb[4] = 4 * ξ * (1-ξ-η)
        Nb[5] = 4 * ξ * η
        Nb[6] = 4 * η * (1-ξ-η)
        dN_dξ = zeros(6)
        dN_dξ[1] = 4*η + 4*ξ - 3
        dN_dξ[2] = 4*ξ - 1
        dN_dξ[3] = 0
        dN_dξ[4] = 4 - 8*ξ - 4*η
        dN_dξ[5] = 4*η
        dN_dξ[6] = -4*η
        dN_dη = zeros(6)
        dN_dη[1] = 4*η + 4*ξ - 3
        dN_dη[2] = 0
        dN_dη[3] = 4*η - 1
        dN_dη[4] = -4*ξ
        dN_dη[5] = 4*ξ
        dN_dη[6] = 4 - 4*ξ - 8*η
        dx_dξ = dN_dξ' * x
        dy_dξ = dN_dξ' * y
        dx_dη = dN_dη' * x
        dy_dη = dN_dη' * y
        J = [ dx_dξ  dy_dξ
              dx_dη  dy_dη ]
        dN_dx = zeros(6)
        dN_dy = zeros(6)
        for ii = 1:6
            Temp = J \ [ dN_dξ[ii]; dN_dη[ii] ]
            dN_dx[ii] = Temp[1,1]
            dN_dy[ii] = Temp[2,1]
        end
        N = kron(transpose(Nb), [ 1 0; 0 1 ])
        Ke += c_air^2*t*(transpose(dN_dx)*dN_dx + transpose(dN_dy)*dN_dy)*abs(det(J))*H
        Me += t*transpose(N)*N*abs(det(J))*H
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
            Ce += t*transpose(N)*N*abs(det(J))*H * eta_pml
        end
        if det(J) < 0
            println("\nAttension: The determinate of the Jacobian is negative!")
        end
    end

    return Ke, Me, Ce

end
