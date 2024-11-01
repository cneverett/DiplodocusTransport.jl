function MaxwellJuttner_Distribution(meanp::Vector{Float32},T::Float32,mass::Float32;n::Float32=1f6)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    T64 = Float64(T)
    m64 = Float64(mass)
    n64 = Float64(n)

    meanp64 = Float64.(meanp)

    MJPoints = zeros(Float64,size(meanp))
    # besselk doesn't do well with small T besselkx = besselk * e^x does better.
    theta = mEle*c^2/(kb*T64)
    @. MJPoints = (meanp64^2)*(n64*theta/(2*m64^2))* (1/besselkx(2,m64*theta)) * exp((m64-sqrt(m64^2+meanp64^2))*theta)

    #@. MJPoints = (c * meanp64^2)/(kb*T*besselkx(2,1/theta)) * exp((1-sqrt(1+meanp64^2))/theta)

    #return  @. MJPoints * (MJPoints >= 1e-30) 

    return Float32.(MJPoints)

end