function cal_Me_bar_2DBar(Le,ρ,A)

    Me_bar = ρ*A/Le*[ 1/3   0  1/6   0
                        0   0    0   0
                      1/6   0  1/3   0
                        0   0    0   0 ]

    return Me_bar

end
