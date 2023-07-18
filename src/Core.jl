mutable struct Core

    amo :: AdvancedMatrixOperators
    psi_solver :: StreamfunctionSolver
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)
        psi_solver = StreamfunctionSolver(amo=amo, pp=ev.pp)

        M_Ψ = - 1 / (ev.pp.J * ev.pp.k^4) * sparse( [ amo.V_mask_V * amo.V_interp_T * amo.T_f_T  amo.V_∂y_T ] )

        ops = Dict(
            :ydiff    => amo.T_DIVy_V * spdiagm(0=>ev.pp.L) * amo.V_∂y_T,
            :M_Ψ      => M_Ψ,
        )

        return new(
            amo,
            psi_solver,
            ops,
        )    
    end
end
