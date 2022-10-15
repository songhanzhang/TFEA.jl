function cal_Ke_bar_3DEulerBeam(Le,E,G,A,Ixx,Iyy,Izz)

    Ke_bar = zeros(12,12)

    Ke_axial = E*A/Le*[  1  -1
                        -1   1 ]
    sel_DOF = [ 1, 7 ]
    Ke_bar[sel_DOF,sel_DOF] += Ke_axial

    Ke_torsion = G*Ixx/Le*[  1  -1
                            -1   1 ]
    sel_DOF = [ 4, 10 ]
    Ke_bar[sel_DOF,sel_DOF] += Ke_torsion

    Ke_b_in = [ 12*E*Izz/Le^3   6*E*Izz/Le^2  -12*E*Izz/Le^3   6*E*Izz/Le^2
                 6*E*Izz/Le^2   4*E*Izz/Le     -6*E*Izz/Le^2   2*E*Izz/Le
               -12*E*Izz/Le^3  -6*E*Izz/Le^2   12*E*Izz/Le^3  -6*E*Izz/Le^2
                 6*E*Izz/Le^2   2*E*Izz/Le     -6*E*Izz/Le^2   4*E*Izz/Le ]
    sel_DOF = [ 2, 6, 8, 12 ]
    Ke_bar[sel_DOF,sel_DOF] += Ke_b_in

    Ke_b_out = [ 12*E*Iyy/Le^3   -6*E*Iyy/Le^2  -12*E*Iyy/Le^3   -6*E*Iyy/Le^2
                 -6*E*Iyy/Le^2    4*E*Iyy/Le      6*E*Iyy/Le^2    2*E*Iyy/Le
                -12*E*Iyy/Le^3    6*E*Iyy/Le^2   12*E*Iyy/Le^3    6*E*Iyy/Le^2
                 -6*E*Iyy/Le^2    2*E*Iyy/Le      6*E*Iyy/Le^2    4*E*Iyy/Le ]
    sel_DOF = [ 3, 5, 9, 11 ]
    Ke_bar[sel_DOF,sel_DOF] += Ke_b_out

    return

end
