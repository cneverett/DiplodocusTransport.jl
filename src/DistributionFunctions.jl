function MaxwellJuttner(meanp::Vector{Float32},T::Float32,mass::Float32)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    mc2 = 8.199e-14
    c = 3e8
    kb = 1.38e-23
    mc2divkb = mc2/kb

    T64 = Float64(T)
    mass64 = Float64(mass)

    meanp64 = Float64.(meanp)

    MJPoints = zeros(Float64,size(meanp))
    # besselk doesn't do well with small T besselkx = besselk * e^x does better.
    @. MJPoints = (c * meanp64^2)/(kb*T*besselkx(2,mass64*mc2divkb/T64)) * exp((1-sqrt(1+meanp64^2))*mass64*mc2divkb/T64)

    #return  @. MJPoints * (MJPoints >= 1e-30) 
    return Float32.(MJPoints)

end