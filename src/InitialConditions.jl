function InitialConditions(Lists::ListStruct)

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    num_species = length(name_list)
    f0_list2D = Vector{Array{Float32,2}}(undef,num_species)
    for i in 1:num_species
        f0_list2D[i] = fill(Float32(0),p_num_list[i],u_num_list[i])
    end

    for i in eachindex(name_list)
        f0_list2D[i][40:41,:] .= 1f3 
    end

    f0_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f0_list[i] = reshape(f0_list2D[i],p_num_list[i]*u_num_list[i])
    end

    u0 = zeros(Float32,sum(p_num_list.*u_num_list))

    f_list_to_u!(u0,f0_list)

    #test = reshape(u0,(p_num_list[1],u_num_list[1]))
    
    #test[40:41,:] .= 1f3

    #u0 .= reshape(test,prod(size(test)))

    return u0

end

"""
    Initial_PowerLaw(Lists,species,pmin,pmax,index,num_Init)

A power-law distribution is typically defined by N(E) ∝ E^(-index). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index) * p/E = pE^(-index-1). Averaging this over a cell gives f(p)_avg = [E^(1-index)/(1-index)]/[p] where [] denote evalution at the cell bounds.
"""
function Initial_PowerLaw(Lists::ListStruct,species::String,pmin::T,pmax::T,umin::T,umax::T,index::Float32,num_Init::Float32) where T <: Union{Float32,Int64}

    
    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list

    species_index = findfirst(==(species),name_list)
    u0_2D_species = zeros(Float32,p_num_list[species_index],u_num_list[species_index])

    # Set initial conditions goes here
    pu = p_up_list[species_index]
    pl = p_low_list[species_index]
    p_num = p_num_list[species_index]
    p_grid = p_grid_list[species_index]
    u_num = u_num_list[species_index]
    u_grid = u_grid_list[species_index]
    pr = BCI.bounds(pl,pu,p_num,p_grid)
    dp = BCI.deltaVector(pr)
    du = BCI.deltaVector(BCI.bounds(BCI.u_low,BCI.u_up,u_num,u_grid))
    mass = getfield(BCI,Symbol("mu"*name_list[species_index]))

    type = zero(T)
    if typeof(type)==Float32
        pmin_index = location(pl,pu,p_num,pmin,p_grid)
        pmax_index = location(pl,pu,p_num,pmax,p_grid)
        umin_index = location(BCI.u_low,BCI.u_up,u_num,umin,u_grid)
        umax_index = location(BCI.u_low,BCI.u_up,u_num,umax,u_grid)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
    end
    

    # power law averaged over cell width.
    for i in pmin_index:pmax_index, j in umin_index:umax_index
        u0_2D_species[i,j] = sqrt(mass^2+pr[i+1]^2)^(1-index) - sqrt(mass^2+pr[i]^2)^(1-index)
        u0_2D_species[i,j] /= 1-index
        u0_2D_species[i,j] /= (pr[i+1]-pr[i])
    end
    # set values and normlaise to initial number density (in m^{-3})
    num = dp' * u0_2D_species * du
    u0_2D_species *= num_Init/num

    # scale by dp*du 
    for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
        u0_2D_species[i,j] *= dp[i] * du[j]
    end

    u0_2D_species = Float32.(u0_2D_species)

    u0_species = reshape(u0_2D_species,p_num*u_num)

    return u0_species
end

function Initial_Constant(Lists::ListStruct,species::String,pmin::T,pmax::T,umin::T,umax::T,num_Init::Float32;mode="AXI") where T <: Union{Float32,Int64}

    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list

    species_index = findfirst(==(species),name_list)
    u0_2D_species = zeros(Float32,p_num_list[species_index],u_num_list[species_index])

    pu = p_up_list[species_index]
    pl = p_low_list[species_index]
    p_grid = p_grid_list[species_index]
    p_num = p_num_list[species_index]
    u_grid = u_grid_list[species_index]
    u_num = u_num_list[species_index]

    type = zero(T)
    if typeof(type)==Float32
        pmin_index = BCI.location(pl,pu,p_num,pmin,p_grid)
        pmax_index = BCI.location(pl,pu,p_num,pmax,p_grid)
        umin_index = BCI.location(BCI.u_low,BCI.u_up,u_num,umin,u_grid)
        umax_index = BCI.location(BCI.u_low,BCI.u_up,u_num,umax,u_grid)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
    end

    dp = BCI.deltaVector(BCI.bounds(pl,pu,p_num,p_grid))
    du = BCI.deltaVector(BCI.bounds(BCI.u_low,BCI.u_up,u_num,u_grid))

    # set values and normlaise to initial number density (in m^{-3})
    for i in pmin_index:pmax_index, j in umin_index:umax_index
        u0_2D_species[i,j] = 1e0
        #u0_2D_species[i,j] = dp[i] * du[j]
    end
    num = dp' * u0_2D_species * du
    #num = sum(u0_2D_species)
    u0_2D_species *= num_Init/num
    u0_2D_species = Float32.(u0_2D_species)

    # scale by dp*du 
    for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
        u0_2D_species[i,j] *= dp[i] * du[j]
    end

    if mode=="AXI"
        u0_species = reshape(u0_2D_species,p_num*u_num)
    elseif mode=="ISO"
        # f(p) = 2*f(p,μ)
        u0_species = dropdims(sum(u0_2D_species,dims=2),dims=2) * 2 / u_num_list[species_index]
    end

    return u0_species
end

function Initial_Temperature(Lists::ListStruct,species::String,T::Float32,num_Init::Float32;mode="AXI")
    
    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list

    species_index = findfirst(==(species),name_list)
    pu = p_up_list[species_index]
    pl = p_low_list[species_index]
    p_grid = p_grid_list[species_index]
    p_num = p_num_list[species_index]
    u_num = u_num_list[species_index]
    u_grid = u_grid_list[species_index]

    u0_2D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index])

    meanp = BCI.meanVector(BCI.bounds(pl,pu,p_num,p_grid))
    dp = BCI.deltaVector(BCI.bounds(pl,pu,p_num,p_grid))
    du = BCI.deltaVector(BCI.bounds(BCI.u_low,BCI.u_up,u_num,u_grid))
    mass = getfield(BCI,Symbol("mu"*name_list[species_index]))

    u0_2D_species .= MaxwellJuttner_Distribution(Float32.(meanp),T,Float32(mass))
    
    # set values and normalise to initial number density (in m^{-3})
    num = dp' * u0_2D_species * du
    u0_2D_species *= num_Init/num

    # scale by dp*du 
    for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
        u0_2D_species[i,j] *= dp[i] * du[j]
    end

    u0_2D_species = Float32.(u0_2D_species)

    if mode=="AXI"
        u0_species = reshape(u0_2D_species,p_num_list[species_index]*u_num_list[species_index])
    elseif mode=="ISO"
        # f(p) = 2*f(p,μ)
        u0_species = dropdims(sum(u0_2D_species,dims=2),dims=2) * 2 / u_num_list[species_index]
    end

    return u0_species

end

