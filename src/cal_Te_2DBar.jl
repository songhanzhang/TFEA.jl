function cal_Te_2DBar(e_nodes)

    xi = e_nodes[1,1]
    yi = e_nodes[1,2]
    xj = e_nodes[2,1]
    yj = e_nodes[2,2]

    ss = (yj-yi)/Le
    cs = (xj-xi)/Le
    
    Te = [ cs   ss    0    0
          -ss   cs    0    0
            0    0   cs   ss
            0    0  -ss   cs ]

    return Te

end
