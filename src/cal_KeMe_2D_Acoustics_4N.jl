function cal_KeMe_2D_Acoustics_4N(x,y,c_air,t)

    Gauss = [ -1/(sqrt(3))  -1/(sqrt(3))  1
               1/(sqrt(3))  -1/(sqrt(3))  1
               1/(sqrt(3))   1/(sqrt(3))  1
              -1/(sqrt(3))   1/(sqrt(3))  1 ]
    Ke = zeros(4,4)
    Me = zeros(4,4)
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
        dx_dξ = dN1_dξ*x[1] + dN2_dξ*x[2] + dN3_dξ*x[3] + dN4_dξ*x[4]
        dy_dξ = dN1_dξ*y[1] + dN2_dξ*y[2] + dN3_dξ*y[3] + dN4_dξ*y[4]
        dx_dη = dN1_dη*x[1] + dN2_dη*x[2] + dN3_dη*x[3] + dN4_dη*x[4]
        dy_dη = dN1_dη*y[1] + dN2_dη*y[2] + dN3_dη*y[3] + dN4_dη*y[4]
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
        N = [ N1  N2  N3  N4 ]
        Ke += c_air^2*t*(transpose(dN_dx)*dN_dx + transpose(dN_dy)*dN_dy)*abs(det(J))*H
        Me += t*transpose(N)*N*abs(det(J))*H
    end

    return Ke, Me

end
