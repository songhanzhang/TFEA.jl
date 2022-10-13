function cal_Te_2DBar(xi,yi,xj,yj)

    Le = sqrt((xj-xi)^2 + (yj-yi)^2)

    ss = (yj-yi)/Le
    cs = (xj-xi)/Le

    Te = [ cs   ss    0    0
          -ss   cs    0    0
            0    0   cs   ss
            0    0  -ss   cs ]

    return Te

end
