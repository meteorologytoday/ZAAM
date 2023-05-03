mutable struct Core

    amo :: AdvancedMatrixOperators
    psi_solver :: StreamfunctionSolver
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)
        psi_solver = StreamfunctionSolver(amo=amo, pp=ev.pp)

        ops = Dict(
            :ydiff    => amo.T_DIVy_V * ev.pp.L * amo.V_∂y_T,
            :V_solve_T => psi_solver.V_solve_T,
        )

        return new(
            amo,
            psi_solver,
            ops,
        )    
    end
end
