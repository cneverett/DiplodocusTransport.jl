"""
    numberDensity(f0,dp,dμ)

Returns the average number density of a distribution function `f`. Number density is defined as the zeroth moment of the distribution function. i.e. 
    ```math
    n = \\int \\mathrm{d}p\\mathrm{d}\\mu f(p,μ) = \\sum_{i,j} f_{ij} \\Delta p_i \\Delta μ_j
    ``` 
where `dp` = ``\\Delta p_i`` is a vector of momentum intervals and `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals.
"""
function NumberDensity(f::Array{Float32},dp::Vector{Float32},dμ::Vector{Float32})

    n = dp' * f * dμ
    
    return n
    
end

"""
    momentum(f0,dp,meanp,dμ)

Returns the average momentum of a distribution function `f`. average momentum is defined as the first moment of the distribution function. i.e. 
    ```math
    \\braket{p} = \\frac{\\int \\mathrm{d}p\\mathrm{d}\\mu pf(p,μ)}{n} = \\sum_{i,j} f_{ij} \\Delta p_i\\braket{p}_i \\Delta μ_j / n
    ``` 
where `dp` = ``\\Delta p_i`` is a vector of momentum intervals, `meanp` = ``\\braket{p}_i`` is a vector of the average momentum value per bin, `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals, and `numberDensity` = ``n`` is the average number density calculated using the function [`numberDensity`](@ref).
"""
function Momentum(f::Array{Float32},dp::Vector{Float32},meanp::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32)

    dpmeanp = dp .* meanp
    momentum = dpmeanp' * f * dμ
    momentum /= numberDensity

    return momentum
end


"""
    Energy(f0,ΔE,dμ,numberDensity)

Returns the average TOTAL energy of a distribution function `f`. average energy is defined as the first moment of the distribution function. i.e. 
    ```math
    \\braket{p} = \\frac{\\int \\mathrm{d}p\\mathrm{d}\\mu p^0f(p,μ)}{n} = \\sum_{i,j} f_{ij} \\Delta E_i \\Delta μ_j / n
    ``` 
where `ΔE` = ``\\Delta E_i`` is a vector of the average "energy" value per bin (has dimensions of momentum squared), `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals, and `numberDensity` = ``n`` is the average number density calculated using the function [`numberDensity`](@ref).
"""
function Energy(f::Array{Float32},ΔE::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32)

    energy = ΔE' * f * dμ
    energy /= numberDensity

    return energy
    
end

"""
    Temperature(f,ΔEkin,dμ,numberDensity)

Returns the Temperature of a distribution function `f`. Average, kinetic energy is defined from the first moment of the distribution function minus the zeroth moment times the mass, i.e. ``p^0-m``: 
    ```math
    \\braket{Ekin} = \\frac{\\int \\mathrm{d}p\\mathrm{d}\\mu (p^0-m)f(p,μ)}{n} = \\sum_{i,j} f_{ij} \\Delta Ekin_i \\Delta μ_j / n
    ``` 
where `ΔEkin` = ``\\Delta Ekin_i`` is a vector of the average "kinetric energy" value per bin (has dimensions of momentum squared), `dp` = ``\\Delta p_i``. `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals, and `numberDensity` = ``n`` is the average number density calculated using the function [`numberDensity`](@ref).
"""
function Temperature(f::Array{Float32},ΔEkin::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32)

    Ekin = transpose(ΔEkin) * f * dμ
    Ekin /= numberDensity

    #  Ekin = 3/2 kB * T 
    # Ekin = E - mass 
    #println(Ekin)
    # units of Kelvin
    Temperature = (Ekin) * (9.11e-31 * 3f8^2 / 1.38f-23) * (2/3)  # E = 3 kb T 

    return Temperature
    
end

