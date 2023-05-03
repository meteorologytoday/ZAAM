mutable struct PhyParams

    J :: Float64   # momentum vertical diffusivity 
    L :: Float64   # momentum horizontal diffusivity
    K :: Float64   # buoyancy horizontal diffusivity
    N :: Float64   # Stability of the atmosphere 
    C :: Float64   
    k :: Float64   

end


mutable struct Env


    pp :: PhyParams
    gd :: Grid

    function Env(;
        ϕs :: Float64,
        ϕn :: Float64,
        npts :: Integer,
        J :: Float64,
        L :: Float64,
        K :: Float64,
        H :: Float64,
        N :: Float64,
        C :: Float64,
    )

        pp = PhyParams(J, L, K, N, C, π/H)

        gd = Grid(
            Ny=npts,
            Ω=Ω,
            ϕn = ϕn,
            ϕs = ϕs,
            H = H,
            R = R,
        )

        return new(
            pp, gd,
        )
        
    end
end
