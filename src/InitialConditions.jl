function InitialConditions(Lists)

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    num_species = length(name_list)
    f0_list2D = Vector{Array{Float32,2}}(undef,num_species)
    for i in 1:num_species
        f0_list2D[i] = fill(Float32(0),nump_list[i],numt_list[i])
    end

    for i in eachindex(name_list)
        f0_list2D[i][40:41,:] .= 1f3 
    end

    f0_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f0_list[i] = reshape(f0_list2D[i],nump_list[i]*numt_list[i])
    end

    u0 = zeros(Float32,sum(nump_list.*numt_list))

    f_list_to_u!(u0,f0_list)

    #test = reshape(u0,(nump_list[1],numt_list[1]))
    
    #test[40:41,:] .= 1f3

    #u0 .= reshape(test,prod(size(test)))

    return u0

end

function Inital_PowerLaw(Lists,species::String,pmin::Float32,pmax::Float32,index::Float32,max_density::Float32)

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    species_index = findfirst(==(species),name_list)
    u0_2D_species = Array{Float32,2}(undef,nump_list[species_index],numt_list[species_index])

    # Set initial conditions goes here
    pu = pu_list[species_index]
    pl = pl_list[species_index]
    nump = nump_list[species_index]

    pmin_index = location_p(pu,pl,nump,pmin)
    pmax_index = location_p(pu,pl,nump,pmax)
    prange = BoltzmannCollisionIntegral.prange(pl,pu,nump)

    # power law p^(-index) averaged over cell width.
    for i in pmin_index:pmax_index
        u0_2D_species[i,:] .= (prange[i+1]^(1-index)-prange[i]^(1-index))/((1-index)*(prange[i+1]-prange[i]))
    end
    # scale correctly
    max = maximum(u0_2D_species)
    u0_2D_species .*= max_density/max

    u0_species = reshape(u0_2D_species,nump_list[species_index],numt_list[species_index])

    return u0_species
end

function Initial_Constant(Lists,species::String,pmin::T,pmax::T,umin::T,umax::T,num_Init::Float32;mode="AXI") where T <: Union{Float32,Int64}

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    species_index = findfirst(==(species),name_list)
    u0_2D_species = zeros(Float32,nump_list[species_index],numt_list[species_index])

    pu = pu_list[species_index]
    pl = pl_list[species_index]
    nump = nump_list[species_index]
    numt = numt_list[species_index]

    type = zero(T)
    if typeof(type)==Float32
        pmin_index = location_p(pu,pl,nump,pmin)
        pmax_index = location_p(pu,pl,nump,pmax)
        umin_index = location_t(numt,umin)
        umax_index = location_t(numt,umax)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
    end

    Δp = BoltzmannCollisionIntegral.deltaVector(BoltzmannCollisionIntegral.prange(pl,pu,nump))
    Δμ = BoltzmannCollisionIntegral.deltaVector(BoltzmannCollisionIntegral.trange(numt))

    # set values and normlaise to initial number density (in m^{-3})
    for i in pmin_index:pmax_index, j in umin_index:umax_index
        u0_2D_species[i,j] = 1e0
        #u0_2D_species[i,j] = Δp[i] * Δμ[j]
    end
    num = Δp' * u0_2D_species * Δμ
    #num = sum(u0_2D_species)
    u0_2D_species *= num_Init/num
    u0_2D_species = Float32.(u0_2D_species)

    if mode=="AXI"
        u0_species = reshape(u0_2D_species,nump_list[species_index]*numt_list[species_index])
    elseif mode=="ISO"
        # f(p) = 2*f(p,μ)
        u0_species = dropdims(sum(u0_2D_species,dims=2),dims=2) * 2 / numt_list[species_index]
    end

    return u0_species
end

function Initial_Temperature(Lists,species::String,T::Float32,num_Init::Float32;mode="AXI")
    
    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    species_index = findfirst(==(species),name_list)
    pu = pu_list[species_index]
    pl = pl_list[species_index]
    nump = nump_list[species_index]
    numt = numt_list[species_index]

    u0_2D_species = zeros(Float32,nump_list[species_index],numt_list[species_index])

    meanp = BoltzmannCollisionIntegral.meanVector(BoltzmannCollisionIntegral.prange(pl,pu,nump))
    Δp = BoltzmannCollisionIntegral.deltaVector(BoltzmannCollisionIntegral.prange(pl,pu,nump))
    Δμ = BoltzmannCollisionIntegral.deltaVector(BoltzmannCollisionIntegral.trange(numt))
    mass = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name_list[species_index]))

    u0_2D_species .= MaxwellJuttner(meanp,T,Float32(mass))

    # set values and normlaise to initial number density (in m^{-3})
    num = Δp' * u0_2D_species * Δμ
    u0_2D_species *= num_Init/num

    #for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
    #    u0_2D_species[i,j] *= Δp[i] * Δμ[j]
    #end

    u0_2D_species = Float32.(u0_2D_species)


    if mode=="AXI"
        u0_species = reshape(u0_2D_species,nump_list[species_index]*numt_list[species_index])
    elseif mode=="ISO"
        # f(p) = 2*f(p,μ)
        u0_species = dropdims(sum(u0_2D_species,dims=2),dims=2) * 2 / numt_list[species_index]
    end

end

function location_p(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform, val in log10
    logp = val
    loc = logp != l ? ceil(Int64,Float32(num)*(logp-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc : loc>num ? num+1 : 1 # assignes 1 for under, num+1 for over and loc for in range
end

function location_t(numt::Int64,val::Float64)
    # function for generating poisition in array. Bins MUST be uniform
    return val != tl ? ceil(Int64,Float64(numt)*(val-tl)/(tu-tl)) : Int64(1) 
end
