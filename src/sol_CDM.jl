function sol_CDM(Mg,Cg,Kg,Fg,t_ax)

    ΔT = t_ax[2] - t_ax[1]
    n_DOF = size(Mg,1)
    Ns = length(t_ax)
    Ug = zeros(n_DOF,length(t_ax))

    @printf "\nDirect time integral (central difference method) started ... \n"
    for i_t = 3:length(t_ax)
        if mod(i_t,10) == 0
            @printf "    %.3f%% completed ... \n" i_t/Ns*100
        end
            K_eq = Mg/ΔT^2 + Cg/(2*ΔT)
            F_eq = Fg[:,i_t-1] + (2*Mg/ΔT^2-Kg)*Ug[:,i_t-1] - (Mg/ΔT^2-Cg/(2*ΔT))*Ug[:,i_t-2]
            Ug[:,i_t] = K_eq\F_eq
    end
    @printf "\nAnalysis completed!\n"

    return Ug

end
