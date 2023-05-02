include("ZAAM.jl")




m = ZAAM.Model(ZAAM.Env(
    ϕn =  deg2rad(90),
    ϕs =  deg2rad(-90),
    npts = 100,
    J = 1e3,
    L = 1e5,
    K = 1e5,
    H = 10e3, 
    N = 1e-6,
    τ = 1e6,
    C = 1e-2,
))


B = m.st.B
B[50] = 1.0

dt = 86400.0
steps = 100

save_B = zeros(Float64, length(B), steps)


for t=1:steps
   
    save_B[:, t] = B
    
    println("Step Model: ", t) 

    ZAAM.step!(m, dt)


end


println("Try to output data...")

println("Loading libraries...")
using YAXArrays, NetCDF
using Dates
println("Done")

axlist = [
    RangeAxis("y", collect(1:length(m.ev.gd.ϕ_T))),
    RangeAxis("time", collect(1:steps)),
]

ds = YAXArray(axlist, save_B)

savecube(ds, "test.nc", driver=:netcdf, overwrite=true)
