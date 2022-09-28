function cal_Ke_bar_2DBar(Le,E,A)

    Ke_bar = E*A/Le*[  1   0  -1   0
                       0   0   0   0
                      -1   0   1   0
                       0   0   0   0]

    return Ke_bar

end
