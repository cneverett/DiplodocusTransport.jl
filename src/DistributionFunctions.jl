function MaxwellJuttner_Distribution(meanp::Vector{Float64},T::Float64,mass::Float64;n::Float64=1e6)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    MJPoints = zeros(Float64,size(meanp))
    # besselk doesn't do well with small T besselkx = besselk * e^x does better.
    theta = mEle*c^2/(kb*T)
    @. MJPoints = (meanp^2)*(n*theta/(2*m^2))* (1/besselkx(2,m*theta)) * exp((m-sqrt(m^2+meanp^2))*theta)

    return MJPoints

end