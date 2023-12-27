function cal_KeMe_bar_2DEulerBeam(Le,E,ρ,A,Izz)

    Ke_bar = [ E*A/Le               0              0  -E*A/Le               0              0
                    0   12*E*Izz/Le^3   6*E*Izz/Le^2        0  -12*E*Izz/Le^3   6*E*Izz/Le^2
                    0    6*E*Izz/Le^2     4*E*Izz/Le        0   -6*E*Izz/Le^2     2*E*Izz/Le
              -E*A/Le               0              0   E*A/Le               0              0
                    0  -12*E*Izz/Le^3  -6*E*Izz/Le^2        0   12*E*Izz/Le^3  -6*E*Izz/Le^2
                    0    6*E*Izz/Le^2     2*E*Izz/Le        0   -6*E*Izz/Le^2     4*E*Izz/Le ]

    Me_bar = ρ*A*Le*[ 1/3             0             0  1/6             0             0
                        0         13/35   (11*Le)/210    0          9/70  -(13*Le)/420
                        0   (11*Le)/210      Le^2/105    0   (13*Le)/420     -Le^2/140
                      1/6             0             0  1/3             0             0
                        0          9/70   (13*Le)/420    0         13/35  -(11*Le)/210
                        0  -(13*Le)/420     -Le^2/140    0  -(11*Le)/210      Le^2/105 ]

    return Ke_bar, Me_bar

end
