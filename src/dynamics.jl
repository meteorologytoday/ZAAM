function step!(
    m :: Model,
    dt :: Float64,
)

    st = m.st

    B = st.B
    
    A = getA(m)
    F = getF(m)

    B[:] = ( I - dt * A ) \ ( B + F )
    
end


function getA(
    m :: Model,
)

    return m.co.ops[:ydiff]


end

function getF(
    m :: Model,
)

    return m.st.B * 0


end




function computeT_TOA(
    m :: Model,
    B :: AbstractArray,
)

    pp = m.ev.pp
    gd = m.ev.gd

    return θ0 * (1 .+ (B .+ (pp.N * gd.H) / 2) / g0)

end


function getF_TOAMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd

    Npts = gd.Ny

    # F_TOA

    T_TOA = computeT_TOA(m, m.st.B)
    tmp = σ_boltz * T_TOA.^3
    A  = spdiagm(0 => (- tmp * θ0 / g0))
    F = tmp * (θ0 + pp.N * gd.H / (2 * g0))

    jacobian_A = A

    if return_jacobian
        
        return A, F, jacobian_A
        
    else
        return A, F
    end
end

function getF_AOMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd

    Npts = gd.Ny

    # F_ao

    A  = - pp.C / pp.τ * θ0 / g0 * sparse(I, Npts, Npts)
    F = zeros(Float64, Npts) .+ ( T_ocn - θ0 * ( 1 - N * H / (2 * g0)) ) * pp.C / pp.τ
   
    jacobian_A = A

    if return_jacobian
        
        return A, F, jacobian_A
        
    else
        return A, F
    end
end
