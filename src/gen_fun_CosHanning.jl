function gen_fun_CosHanning(t_ax,fc,n_peaks,T_0)

    Ns = length(t_ax)
    ΔT = t_ax[2] - t_ax[1]
    n_steps = Int(floor(n_peaks/fc/ΔT))

    y = zeros(Ns)
    for ii = T_0 : T_0 + n_steps
        y[ii] = (1-cos(2*pi*fc*(ii-T_0)*ΔT/n_peaks))*sin(2*pi*fc*(ii-T_0)*ΔT)
    end

    return y

end
