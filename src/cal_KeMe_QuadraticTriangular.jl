function cal_KeMe_QuadraticTriangular(x,y,E,ν,ρ,type)

    Gauss = [ 0.0915762135  0.8168475730  0.1099517437
              0.0915762135  0.0915762135  0.1099517437
              0.8168475730  0.0915762135  0.1099517437
              0.4459484909  0.1081030182  0.2233815897
              0.4459484909  0.4459484909  0.2233815897
              0.1081030182  0.4459484909  0.2233815897 ];
    λ = E*ν/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    if type == "PlaneStress"
        D = E/(1-ν^2)*[ 1 ν 0
                        ν 1 0
                        0 0 (1-ν)/2 ]
    elseif type == "PlaneStrain"
        D = [ λ+2*μ  λ      0
              λ      λ+2*μ  0
              0      0      μ ]
    else
        println("\nPlease specify the plane problem type, plane stress or plane strain?")
    end
    Ke = zeros(12,12)
    Me = zeros(12,12)
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
        B = [ dN_dx[1] 0 dN_dx[2] 0 dN_dx[3] 0 dN_dx[4] 0 dN_dx[5] 0 dN_dx[6] 0
              0 dN_dy[1] 0 dN_dy[2] 0 dN_dy[3] 0 dN_dy[4] 0 dN_dy[5] 0 dN_dy[6]
              dN_dy[1] dN_dx[1] dN_dy[2] dN_dx[2] dN_dy[3] dN_dx[3] dN_dy[4] dN_dx[4] dN_dy[5] dN_dx[5] dN_dy[6] dN_dx[6] ]
        N = kron(transpose(Nb), [ 1 0; 0 1 ])
        Ke += transpose(B)*D*B*det(J)*H
        Me += ρ*transpose(N)*N*det(J)*H
        if det(J) < 0
            println("\nAttension: The determinate of the Jacobian is negative!")
        end
    end

    return Ke, Me

end
