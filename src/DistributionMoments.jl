"""
    numberDensity(u,dp,dμ;mode="AXI")

Returns the average number density of a distribution function `u` output from solver. Number density is defined as the zeroth moment of the distribution function. i.e. 
    ```math
    n = \\int \\mathrm{d}p\\mathrm{d}\\mu f(p,μ) = \\sum_{i,j} f_{ij} \\Delta p_i \\Delta μ_j
    ``` 
where `dp` = ``\\Delta p_i`` is a vector of momentum intervals and `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals.
"""
function NumberDensity(u::Vector{Float32},nump,numt,dp::Vector{Float32},dμ::Vector{Float32};mode="AXI")

    if mode=="AXI"
        f = copy(reshape(u,(nump,numt)))
        #for i in axes(f,1), j in axes(f,2)
        #    f[i,j] /= dp[i] * dμ[j]
        #end
        n = dp' * f * dμ
    elseif mode=="ISO"
        f = u
        n = dp' * f * 2 # 2 is total range of dμ
    end  
    
    return n
    
end

"""
    momentum(f0,dp,meanp,dμ;mode="AXI")

Returns the average momentum of a distribution function `u` output from solver. average momentum is defined as the first moment of the distribution function. i.e. 
    ```math
    \\braket{p} = \\frac{\\int \\mathrm{d}p\\mathrm{d}\\mu pf(p,μ)}{n} = \\sum_{i,j} f_{ij} \\Delta p_i\\braket{p}_i \\Delta μ_j / n
    ``` 
where `dp` = ``\\Delta p_i`` is a vector of momentum intervals, `meanp` = ``\\braket{p}_i`` is a vector of the average momentum value per bin, `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals, and `numberDensity` = ``n`` is the average number density calculated using the function [`numberDensity`](@ref).
"""
function Momentum(u::Vector{Float32},nump,numt,dp::Vector{Float32},meanp::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32;mode="AXI")

    
    dpmeanp = dp .* meanp
    
    if mode=="AXI"
        f = copy(reshape(u,(nump,numt)))
        #for i in axes(f,1), j in axes(f,2)
        #    f[i,j] /= dp[i] * dμ[j]
        #end
        momentum = dpmeanp' * f * dμ
    elseif mode=="ISO"
        f = u
        momentum = dpmeanp' * f * 2 # 2 is total range of dμ
    end
   
    momentum /= numberDensity

    return momentum
end


"""
    Energy(f0,ΔE,dμ,numberDensity)

Returns the average TOTAL energy of a distribution function `u` output from solver. average energy is defined as the first moment of the distribution function. i.e. 
    ```math
    \\braket{p} = \\frac{\\int \\mathrm{d}p\\mathrm{d}\\mu p^0f(p,μ)}{n} = \\sum_{i,j} f_{ij} \\Delta E_i \\Delta μ_j / n
    ``` 
where `ΔE` = ``\\Delta E_i`` is a vector of the average "energy" value per bin (has dimensions of momentum squared), `dμ` = ``\\Delta μ_j`` is a vector of cosine (momentum space) angle intervals, and `numberDensity` = ``n`` is the average number density calculated using the function [`numberDensity`](@ref).
"""
function Energy(u::Vector{Float32},nump,numt,ΔE::Vector{Float32},dp::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32;mode="AXI")


    if mode=="AXI"
        f = copy(reshape(u,(nump,numt)))
        #for i in axes(f,1), j in axes(f,2)
        #    f[i,j] /= dp[i] * dμ[j]
        #end
        energy = ΔE' * f * dμ
    elseif mode=="ISO"
        f = u
        energy = ΔE' * f * 2 # 2 is total range of dμ
    end

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
function Temperature(u::Vector{Float32},nump,numt,ΔEkin::Vector{Float32},dp::Vector{Float32},dμ::Vector{Float32},numberDensity::Float32;mode="AXI")

    if mode=="AXI"
        f = copy(reshape(u,(nump,numt)))
        #for i in axes(f,1), j in axes(f,2)
        #    f[i,j] /= dp[i] * dμ[j]
        #end
        Ekin = transpose(ΔEkin) * f * dμ
    elseif mode=="ISO"
        f = u
        Ekin = transpose(ΔEkin) * f * 2 # 2 is total range of dμ
    end

    Ekin /= numberDensity

    #  Ekin = 3/2 kB * T 
    # Ekin = E - mass 
    #println(Ekin)
    # units of Kelvin
    Temperature = (Ekin) * (9.11e-31 * 3f8^2 / 1.38f-23) * (2/3)  # E = 3 kb T 

    return Temperature
    
end

"""
    FourFlow(u,dp,du,dE,du2)

Returns the four flow vector Ua from the flattened axysmmetric distribution function f1D.
"""
function FourFlow(f1D::Vector{Float32},nump,numt,dp::Vector{Float64},du::Vector{Float64},dE::Vector{Float64},du2::Vector{Float64})

    f2D = zeros(Float64,nump,numt)
    @. f2D = Float64(reshape(f1D,(nump,numt)))

    Na = zeros(Float64,4)
    Na[1] = dp' * f2 * du
    Na[4] = dE' * f2 * du2

    return Na

end

function HydroFourVelocity(Na::Vector{Float64})

    Ua = zeros(Float64,4)

    metric = zeros(Float64,4,4)
    metric[1,1] = -1
    metric[2,2] = 1
    metric[3,3] = 1
    metric[4,4] = 1

    sqrtNa2 = sqrt(Na'*metric*Na)
    Ua .= Na/sqrtNa2

    return Ua

end

function ProjectionTensor(Ua::Vector{Float64})

    Δab = zeros(Float64,4,4)

    metric = zeros(Float64,4,4)
    metric[1,1] = -1
    metric[2,2] = 1
    metric[3,3] = 1
    metric[4,4] = 1

    for i in 1:4, j in 1:4
        Δab[i,j] = metric[i,j] - Ua[i]*Ua[j]
    end

    return Δab

end

function StressEnergyTensor()

    Tab = zeros(Float64,4,4)

    return Tab
end
