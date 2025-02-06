"""
    getGridValues(Lists)

Returns a tuple of `Vector{Vector{Float32}}` of particle masses (normalised), and momentum magnitude and angle cosine grid bounds, differences, and means for each particle species.
"""
function getGridValues(Lists)

    (name_list,p_up_list,p_low_list,p_grid_list,p_num_list,u_grid_list,u_num_list,interaction_list_Binary,interaction_list_Sync) = Lists

    num_species = length(name_list);
    pr_list = Vector{Vector{Float32}}(undef,num_species);
    ur_list = Vector{Vector{Float32}}(undef,num_species);
    dp_list = Vector{Vector{Float32}}(undef,num_species);
    du_list = Vector{Vector{Float32}}(undef,num_species);
    #ΔE_list = Vector{Vector{Float32}}(undef,num_species);
    #ΔEkin_list = Vector{Vector{Float32}}(undef,num_species);
    meanp_list = Vector{Vector{Float32}}(undef,num_species);
    meanu_list = Vector{Vector{Float32}}(undef,num_species);
    mass_list = Vector{Float32}(undef,num_species);
    for i in eachindex(name_list)
        mass_list[i] = getfield(BCI,Symbol("mu"*name_list[i]));
        pr_list[i] = BCI.bounds(p_low_list[i],p_up_list[i],p_num_list[i],p_grid_list[i]);
        ur_list[i] = BCI.bounds(BCI.u_low,BCI.u_up,u_num_list[i],u_grid_list[i]);
        dp_list[i] = BCI.deltaVector(pr_list[i]);
        du_list[i] = BCI.deltaVector(ur_list[i]);    
        meanp_list[i] = BCI.meanVector(pr_list[i]);
        meanu_list[i] = BCI.meanVector(ur_list[i]);
    end

    return (mass_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list)

end

"""
    getInitialScaling(sol,Lists)

Returns a tuple of `Vector{Float32}` of initial number density, energy density, and temperature values for each particle species, given the output of a simulation `sol` and the input parameters `Lists`.    
"""
function getInitialScaling(sol,Lists,pr_list,ur_list,mass_list)

    (name_list,p_up_list,p_low_list,p_grid_list,p_num_list,u_grid_list,u_num_list,interaction_list_Binary,interaction_list_Sync) = Lists

    num_species = length(name_list);
    numInit_list = zeros(Float32,num_species);
    engInit_list = zeros(Float32,num_species);
    tempInit_list = zeros(Float32,num_species);
    for i in eachindex(name_list)
        Na = FourFlow(sol.u[1].x[i],p_num_list[i],u_num_list[i],pr_list[i],ur_list[i],mass_list[i])
        #println(Na)
        Ua = HydroFourVelocity(Na)
        #println(Ua)
        Δab = ProjectionTensor(Ua)
        #println(Δab)
        Tab = StressEnergyTensor(sol.u[1].x[i],p_num_list[i],u_num_list[i],pr_list[i],ur_list[i],mass_list[i])
        #println(Tab)

        n = ScalarNumberDensity(Na,Ua)
        #println(n)
        e = ScalarEnergyDensity(Tab,Ua,n)
        #println(e)
        p = ScalarPressure(Tab,Δab)
        #println(p)
        T = ScalarTemperature(p,n)
        #println(T)

        numInit_list[i] = n
        engInit_list[i] = e
        tempInit_list[i] = T
    end

    return (numInit_list,engInit_list,tempInit_list)

end


"""
    deltaVector(valr)

Inputs a `num+1` long `Vector{Float}` quantitiy values (domain bounds) and returns a `num` long `Vector{Float}` of differeces (domain widths).

# Examples
```julia-repl
julia> deltaVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0])
3-element Vector{Float64}:
 9.0
 90.0
 900.0
```
"""
function deltaVector(valr::Vector{T}) where T <: Union{Float32,Float64}
    num = size(valr)[1]-1  # number of grid cells
    Δ = zeros(T,num)
    
    for ii in 1:num 
        Δ[ii] += abs(valr[ii+1] - valr[ii]) # abs accounts for Δμ = Δcosθ = cospi(t i) - cospi(t 1+1)
    end

    return Δ
end

# ================================================================ #

# ============== Mean Values in Grid Cells ======================= #
"""
    meanVector(valr)

Inputs a `num+1` long `Vector{Float}` of domain bounds and returns a `num` long `Vector{Float}` of mean value in domain range.

# Examples
```julia-repl
julia> meanVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0])
3-element Vector{Float64}:
 5.5
 55.0
 550.0
```
"""
function meanVector(valr::Vector{T}) where T <: Union{Float32,Float64}
    num = size(valr)[1]-1  # number of grid cells
    mean = zeros(T,num)
    
    for ii in 1:num 
        mean[ii] = (valr[ii+1] + valr[ii])/2
    end

    return mean
end

# =============== Values on Grid Boundaries ====================== #
"""
    prange(pL,pH,num_p)

Returns a `num_p+1` long `Vector{Float}` of p-space grid bounds NOT in Log10 space, between the low and high p bounds `pL` and `pH`.

# Examples
```julia-repl
julia> prange(-5e0,4e0,9)
10-element Vector{Float64}:
 1.0e-5
 1.0e-4
 1.0e-3
 0.01
 0.1
 1.0
 10.0
 100.0
 1000.0
 10000.0
```
"""
function p_range(pL::T,pH::T,num_p::Int64) where T <: Union{Float32,Float64}
    # returns a vector{T} of p grid bounds NOT in Log10 space
    return 10 .^[range(pL,pH,num_p+1);]
end

"""
    u_range(uL,uH,num_u)

Returns a `num_u+1` long `Vector{Float}` of u-space grid bounds NOT in Log10 space, between the low and high u bounds `uL` and `uH`.

# Examples
```julia-repl
julia> u_range(-1.0,1.0,8)
9-element Vector{Float64}:
 -1.0
 -0.75
 -0.5
 -0.25
  0.0
  0.25
  0.5
  0.75
  1.0
```
"""
function u_range(uL::T,uH::T,num_u::Int64) where T <: Union{Float32,Float64}
    return [range(uL,uH,num_u+1);]
end

# ================================================================ #

# =============== Momentum ^ power values on the grid ============ #
function p_power_range(L::T,H::T,num::Int64,power::Int64) where T <: Union{Float32,Float64}

    return 10 .^[range(L,H,num+1);].^power

end


function t_power_range(L::T,H::T,num::Int64,power::Int64) where T <: Union{Float32,Float64}

    return [range(L,H,num+1);].^power
    
end










# ================================================================ #

# ========================= "Delta Energy" ======================= #
"""
    deltaEVector(pr,mu)

Inputs a `num+1` long `Vector{Float}` of p grid boundries and the particle `mu` value (normalised mass) and returns a `num` long `Vector{Float}` of average energy values per grid cell.

# Examples
```julia-repl
julia> deltaEVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0], 1.0e0)
3-element Vector{Float64}:
 50.600693
 4951.15
 495001.16
```
"""
function deltaEVector(pr::Vector{T},mu::T) where T <: Union{Float32,Float64}

    num = size(pr)[1]-1  # number of grid cells
    E = zeros(T,num+1)
    ΔE = zeros(T,num)

    if mu == zero(T)
        for ii in 1:num 
            ΔE[ii] = pr[ii+1]^2-pr[ii]^2
            ΔE[ii] /= 2
        end
    else 
        for ii in 1:num+1 
            E[ii] = mu + pr[ii]^2/(sqrt(mu^2+pr[ii]^2)+mu)
        end 
        for ii in 1:num 
            ΔE[ii] = (pr[ii+1]-pr[ii])*mu
            ΔE[ii] += pr[ii+1]^3/(E[ii+1]+mu) - pr[ii]^3/(E[ii]+mu) 
            ΔE[ii] += mu^2*(asinh(pr[ii+1]/mu)-asinh(pr[ii]/mu))
            ΔE[ii] /= 2
        end
    end 

    return ΔE
end

# ================================================================ #

# ====================== "Delta Kinetic Energy" ================== #
"""
    deltaEkinVector(pr,mu)

Inputs a `num+1` long `Vector{Float}` of p grid boundries and the particle `mu` value (normalised mass) and returns a `num` long `Vector{Float}` of average kinetic energy values per grid cell.

# Examples
```julia-repl
julia> deltaEkinVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0], 1.0e0)
3-element Vector{Float64}:
     46.10069600605712
   4906.1506753523645
 494551.15128635924
```
"""
function deltaEkinVector(pr::Vector{T},mu::T) where T <: Union{Float32,Float64}
    # inputs a (num+1) vector{Float} of p grid boundries and the particle mu value and return a (num) vector{Float} of average energy values per grid cell

    num = size(pr)[1]-1  # number of grid cells
    E = zeros(T,num+1)
    ΔEkin = zeros(T,num)

    if mu == zero(T) # same as ΔE
        for ii in 1:num 
            ΔEkin[ii] = pr[ii+1]^2-pr[ii]^2
            ΔEkin[ii] /= 2
        end
    else 
        for ii in 1:num+1 
            E[ii] = mu + pr[ii]^2/(sqrt(mu^2+pr[ii]^2)+mu)
        end 
        for ii in 1:num 
            # ΔEkin = ΔE - Δp*m
            ΔEkin[ii] = pr[ii+1]^3/(E[ii+1]+mu) - pr[ii]^3/(E[ii]+mu) 
            ΔEkin[ii] += mu^2*(asinh(pr[ii+1]/mu)-asinh(pr[ii]/mu))
            ΔEkin[ii] /= 2
        end
    end 

    return ΔEkin
end

# =============================================================== #