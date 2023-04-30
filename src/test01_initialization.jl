include("ZAAM.jl")




m = ZAAM.Model(ZAAM.Env(
    ϕn =  deg2rad(90),
    ϕs =  deg2rad(-90),
    npts = 20,
    J = 1e3,
    L = 1e5,
    K = 1e5,
    H = 10e3, 
    N = 1e-6,
    τ = 1e6,
    C = 1e-2,
))



