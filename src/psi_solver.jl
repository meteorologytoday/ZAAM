mutable struct StreamfunctionSolver

    amo          :: AdvancedMatrixOperators

    eW_send_W    :: AbstractArray{Float64, 2}
    W_send_eW    :: AbstractArray{Float64, 2}

    V_solveΨ_V    :: AbstractArray{Float64, 2}
    V_solveΨ_T    :: AbstractArray{Float64, 2}
 
    V_solveΓ_V    :: AbstractArray{Float64, 2}
    V_solveΓ_T    :: AbstractArray{Float64, 2}
    
    function StreamfunctionSolver(;
        pp :: PhyParams,
        amo :: AdvancedMatrixOperators,
    )
    
        bmo = amo.bmo
        Ny = bmo.Ny

        # Create eV
        # need a mask excluding bnd points
        mask_eff_V = reshape(amo.V_mask_V * ones(Float64, bmo.V_pts), bmo.V_dim...)

        # Create coversion matrix and its inverse
        V_num = reshape(collect(1:length(mask_eff_V)), size(mask_eff_V)...)
        active_num_eff_V = V_num[ mask_eff_V .==1 ]

        eV_send_V = bmo.V_I_V[active_num_eff_V, :]; dropzeros!(eV_send_V)
        V_send_eV = sparse(eV_send_V')
 
        # construct operators on the left-hand-side
        eV_LAPy_eV = sparse(eV_send_V   * amo.V_LAPy_V   * V_send_eV)

        eV_L_eV = sparse(eV_send_V * ( 
            pp.J * pp.L * pp.k^2 * amo.V_LAPy_V - ( amo.V_f_V * amo.V_f_V .+ pp.J^2 * pp.k^4 )
        ) * V_send_eV)

        eV_invL_eV = inv(Matrix(eV_L_eV))
        
        V_solveΨ_V = - (pp.J * pp.k^4)^(-1) * (bmo.V_I_V + amo.V_f_V * V_send_eV * eV_invL_eV * eV_send_V * amo.V_f_V )
        V_solveΨ_T = V_solveΨ_V * amo.V_∂y_T
        

        V_solveΓ_V = V_send_eV * eV_invL_eV * eV_send_V * amo.V_f_V
        V_solveΓ_T = V_solveΓ_V * amo.V_∂y_T



        return new(
            amo,

            eV_send_V,
            V_send_eV,

            V_solveΨ_V,  # starting from dBdy, for testing
            V_solveΨ_T,

            V_solveΓ_V,  # starting from dBdy, for testing
            V_solveΓ_T,

        ) 

    end
end

