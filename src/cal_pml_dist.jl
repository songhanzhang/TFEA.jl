function cal_pml_eta(pml_interface, model_boundary, x, y, damp_max)

    x_Lb = pml_interface[1]
    x_Rb = pml_interface[2]
    y_Bb = pml_interface[3]
    y_Tb = pml_interface[4]

    x_min = model_boundary[1]
    x_max = model_boundary[2]
    y_min = model_boundary[3]
    y_max = model_boundary[4]

    if x < x_Lb
        eta_pml = damp_max*((x_Lb-x)/(x_Lb-x_min))^3
    elseif x > x_Rb
        eta_pml = damp_max*((x-x_Rb)/(x_max-x_Rb))^3
    elseif y < y_Bb
        eta_pml = damp_max*((y_Bb-y)/(y_Bb-y_min))^3
    elseif y > y_Tb
        eta_pml = damp_max*((y-y_Tb)/(y_max-y_Tb))^3

    return eta_pml

end
