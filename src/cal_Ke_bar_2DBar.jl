function cal_Ke_bar_2DBar(e_nodes,E,A)

    xi = e_nodes[1,1]
    yi = e_nodes[1,2]
    xj = e_nodes[2,1]
    yj = e_nodes[2,2]

    Le = sqrt((xj-xi)^2+(yj-yi)^2)

    Ke_bar = E*A/Le*[  1   0  -1   0
                       0   0   0   0
                      -1   0   1   0
                       0   0   0   0]

    return Ke_bar

end
