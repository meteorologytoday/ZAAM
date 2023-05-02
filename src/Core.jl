mutable struct Core

    amo :: AdvancedMatrixOperators
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)
        ops = Dict(
            :ydiff => amo.T_DIVy_V * ev.pp.L * amo.V_∂y_T,
        )
        return new(
            amo,
            ops,
        )    
    end
end
