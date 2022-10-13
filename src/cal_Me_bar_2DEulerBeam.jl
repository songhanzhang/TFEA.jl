function cal_Me_bar_2DEulerBeam(Le,ρ,A)

    Me_bar = ρ*A*Le * [ 1/3           0          0  1/6           0           0
                          0       13/35  11*Le/210    0        9/70  -13*Le/420
                          0   11*Le/210   Le^2/105    0   13*Le/420   -Le^2/140
                        1/6           0          0  1/3           0           0
                          0        9/70  13*Le/420    0       13/35  -11*Le/210
                          0  -13*Le/420  -Le^2/140    0  -11*Le/210    Le^2/105 ]

    return Me_bar

end
