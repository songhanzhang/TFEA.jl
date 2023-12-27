function cal_Ke_bar_2DEulerBeam(Le,E,A,Izz)

    Ke_bar = [ E*A/Le               0              0  -E*A/Le               0              0
                    0   12*E*Izz/Le^3   6*E*Izz/Le^2        0  -12*E*Izz/Le^3   6*E*Izz/Le^2
                    0    6*E*Izz/Le^2     4*E*Izz/Le        0   -6*E*Izz/Le^2     2*E*Izz/Le
              -E*A/Le               0              0   E*A/Le               0              0
                    0  -12*E*Izz/Le^3  -6*E*Izz/Le^2        0   12*E*Izz/Le^3  -6*E*Izz/Le^2
                    0    6*E*Izz/Le^2     2*E*Izz/Le        0   -6*E*Izz/Le^2     4*E*Izz/Le ]

    return Ke_bar

end
