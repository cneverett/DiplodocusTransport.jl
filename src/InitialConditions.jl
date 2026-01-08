"""
    Initialise_Initial_Condition(PhaseSpace;Precision=Float32)

Returns a zero vector with elements for the initial conditions of all particles at all positions in phase space. 
"""
function Initialise_Initial_Condition(PhaseSpace::PhaseSpaceStruct;Precision::DataType=Float32)
    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num
    
    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*x_num*y_num*z_num

    return zeros(Precision,m)

end


"""
    Initial_PowerLaw!(Initial,PhaseSpace,species,pmin,pmax,umin,uman,hmin,hmax,index,num_Init)

Modifies the initial state vector `Initial` with a power law distribution with `index` for `species`. A power-law distribution is typically defined by N(E) ∝ E^(-index). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index) * p/E = pE^(-index-1). Averaging this over a cell gives f(p)_avg = [E^(1-index)/(1-index)]/[p] where [] denote evaluation at the cell bounds.
"""
function Initial_PowerLaw!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::P,pmax::P,umin::S=-1.0,umax::S=1.0,hmin::S=0.0,hmax::S=2.0,index::AbstractFloat,num_Init::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where P <: Union{Float32,Float64,Int64} where S <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    Grids = PhaseSpace.Grids
    pr_list = Grids.pxr_list
    mp_list = Grids.mpx_list
    mu_list = Grids.mpy_list
    dp_list = Grids.dpx_list
    du_list = Grids.dpy_list
    mass_list = Grids.mass_list

    species_index = findfirst(==(species),name_list)
    pr = pr_list[species_index]
    dp = dp_list[species_index]
    du = du_list[species_index]
    dh = du_list[species_index]

    f0_3D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    # Set initial conditions goes here
    pu = Float64(p_up_list[species_index])
    pl = Float64(p_low_list[species_index])
    p_num = p_num_list[species_index]
    p_grid = p_grid_list[species_index]
    u_num = u_num_list[species_index]
    u_grid = u_grid_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]
    mass = mass_list[species_index]

    typeP = zero(P)
    if (typeof(typeP) == Float32) || (typeof(typeP) == Float64) 
        pmin_index = location(pl,pu,p_num,Float64(pmin),p_grid)
        pmax_index = location(pl,pu,p_num,Float64(pmax),p_grid) 
    elseif typeof(typeP)==Int64
        pmin_index = pmin
        pmax_index = pmax
    end
    typeS = zero(S)
    if (typeof(typeS) == Float32) || (typeof(typeS) == Float64) 
        umin_index = location(CONST_u0,CONST_u1,u_num,Float64(umin),u_grid)
        umax_index = location(CONST_u0,CONST_u1,u_num,Float64(umax),u_grid)
        hmin_index = location(CONST_h0,CONST_h1,h_num,Float64(hmin),h_grid)
        hmax_index = location(CONST_h0,CONST_h1,h_num,Float64(hmax),h_grid)
    elseif typeof(typeS)==Int64
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end
    
    #println(pmin_index,pmax_index,umin_index,umax_index)

    # power law averaged over cell width.
    for px in pmin_index:pmax_index, py in umin_index:umax_index, pz in hmin_index:hmax_index 
        if index != 1
            f0_3D_species[px,py,pz] = sqrt(mass^2+pr[px+1]^2)^(1-index) - sqrt(mass^2+pr[px]^2)^(1-index)
            f0_3D_species[px,py,pz] /= 1-index
            f0_3D_species[px,py,pz] *= du[py] * dh[pz] / (4pi)
        else 
            f0_3D_species[px,py,pz] = log(sqrt(mass^2+pr[px+1]^2)) - log(sqrt(mass^2+pr[px]^2))
            f0_3D_species[px,py,pz] /= 2.0
            f0_3D_species[px,py,pz] *= du[py] * dh[pz] / (4pi)
        end
    end
    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    # scale by dp*du 
    #=for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
        u0_2D_species[i,j] *= dp[i] * du[j]
    end=#

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    Initial_BoostedPowerLaw!(Initial,PhaseSpace,species,pmin,pmax,Gamma,index,num_Init)

Takes an isotropic power-law distribution, with minimum momentum `pmin`, maximum momentum `pmax` and `index` in some frame propagating with Lorentz factor `Gamma` in the local z-direction and modifies the initial state vector (distribution), with a number density of `num_Init`.
"""
function Initial_BoostedPowerLaw!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::S,pmax::S,Gamma::S,index::AbstractFloat,num_Init::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where S <: Union{Float32,Float64} where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    Grids = PhaseSpace.Grids
    species_index = findfirst(==(species),name_list)
    mass = Grids.mass_list[species_index]
    pu = Float64(p_up_list[species_index])
    pl = Float64(p_low_list[species_index])
    p_num = p_num_list[species_index]
    p_grid = p_grid_list[species_index]
    u_num = u_num_list[species_index]
    u_grid = u_grid_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    h_r = Grids.pzr_list[species_index]

    f0_3D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    # Set Rapidity 
    w::Float64 = acosh(Gamma)
    cw = cosh(w)
    sw = sinh(w)
    # de-boost pmin and pmax 
    Emin::Float64 = sqrt(mass^2+pmin^2)
    Emax::Float64 = sqrt(mass^2+pmax^2)
    p1 = sqrt((cw*Emin-sw*pmin)^2-mass^2)
    p2 = sqrt((cw*Emin+sw*pmin)^2-mass^2)
    p3 = sqrt((cw*Emax-sw*pmax)^2-mass^2)
    p4 = sqrt((cw*Emax+sw*pmax)^2-mass^2)
    pmin_UB::Float64 = min(p1,p2,p3,p4)
    pmax_UB::Float64 = max(p1,p2,p3,p4)
    Emin_UB::Float64 = sqrt(mass^2+pmin_UB^2)
    Emax_UB::Float64 = sqrt(mass^2+pmax_UB^2)

    println("Boost Gamma: ",Gamma," w:",w)
    println("p1,p2,p3,p4: ",p1,", ",p2,", ",p3,", ",p4)
    println("Emin: ",Emin," Emax: ",Emax)
    println("Unboosted pmin: ",pmin_UB," Unboosted pmax: ",pmax_UB)
    println("Unboosted Emin: ",Emin_UB," Unboosted Emax: ",Emax_UB)

    pmin_index = location(pl,pu,p_num,pmin_UB,p_grid)
    pmax_index = location(pl,pu,p_num,pmax_UB,p_grid)

    if pmax_index > p_num
        error("Unboosted pmax index exceeds momentum grid size. p_max: ",pmax_UB)
    end

    DistributionToDIP_TrapeziumIntegration!(f0_3D_species,p_r,u_r,h_r,Distribution_Boosted,Gamma,mass,Distribution_PowerLaw,index,pmin,pmax,mass)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    Initial_Constant!(Initial,PhaseSpace,species,pmin,pmax,umin,umax,hmin,hmax,num_Init)

Divides the initial number density `num_Init` equally among momentum-space bins in the range of `pmin` to `pmax`, `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values. This is then applies to the initial state vector `Initial`.
"""
function Initial_Constant!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::T,pmax::T,umin::T,umax::T,hmin::T,hmax::T,num_Init::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where T <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    #Grids = PhaseSpace.Grids
    #dp_list = Grids.dpx_list
    #du_list = Grids.dpy_list

    species_index = findfirst(==(species),name_list)
    #dp = dp_list[species_index]
    #du = du_list[species_index]
    f0_3D_species = zeros(Float32,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    pu = p_up_list[species_index]
    pl = p_low_list[species_index]
    p_grid = p_grid_list[species_index]
    p_num = p_num_list[species_index]
    u_grid = u_grid_list[species_index]
    u_num = u_num_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]

    type = zero(T)
    if (typeof(type) == Float32) || (typeof(type) == Float64) 
        pmin_index = location(pl,pu,p_num,Float64(pmin),p_grid)
        pmax_index = location(pl,pu,p_num,Float64(pmax),p_grid)
        umin_index = location(CONST_u0,CONST_u1,u_num,Float64(umin),u_grid)
        umax_index = location(CONST_u0,CONST_u1,u_num,Float64(umax),u_grid)
        hmin_index = location(CONST_h0,CONST_h1,h_num,Float64(hmin),h_grid)
        hmax_index = location(CONST_h0,CONST_h1,h_num,Float64(hmax),h_grid)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end

    # set values and normalise to initial number density (in m^{-3})
    for px in pmin_index:pmax_index, py in umin_index:umax_index, pz in hmin_index:hmax_index 
        f0_3D_species[px,py,pz] = 1e0
    end

    num = sum(f0_3D_species)
    f0_3D_species *= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    Initial_MaxwellJuttner(PhaseSpace,species,T,umin,umax,hmin,hmax,num_Init)

Modeifies the initial state vector `Initial` with a Maxwell-Juttner distribution for `species` of temperature `T` in Kelvin with a number density of `num_Init` and angular range `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function Initial_MaxwellJuttner!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;T::Float64,umin::S,umax::S,hmin::S,hmax::S,num_Init::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1)  where S <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    species_index = findfirst(==(species),name_list)

    u_grid = u_grid_list[species_index]
    u_num = u_num_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]

    type = zero(S)
    if (typeof(type) == Float32) || (typeof(type) == Float64) 
        umin_index = location(CONST_u0,CONST_u1,u_num,Float64(umin),u_grid)
        umax_index = location(CONST_u0,CONST_u1,u_num,Float64(umax),u_grid)
        hmin_index = location(CONST_h0,CONST_h1,h_num,Float64(hmin),h_grid)
        hmax_index = location(CONST_h0,CONST_h1,h_num,Float64(hmax),h_grid)
    elseif typeof(type)==Int64
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end

    f0_3D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    for py in umin_index:umax_index, pz in hmin_index:hmax_index 
        @view(f0_3D_species[:,py,pz]) .= MaxwellJuttner_Distribution(PhaseSpace,species,T)
    end
    
    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num_list[species_index]*u_num_list[species_index]*h_num_list[species_index])

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing

end

"""
    Initial_BlackBody(PhaseSpace,species,T,umin,umax,hmin,hmax,num_Init)

Modeifies the initial state vector `Initial` with a Black-Body distribution for `species` of temperature `T` in Kelvin with a number density of `num_Init` and angular range `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function Initial_BlackBody!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;T::Float64,umin::S,umax::S,hmin::S,hmax::S,num_Init::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1)  where S <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    species_index = findfirst(==(species),name_list)

    u_grid = u_grid_list[species_index]
    u_num = u_num_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]

    type = zero(S)
    if (typeof(type) == Float32) || (typeof(type) == Float64) 
        umin_index = location(CONST_u0,CONST_u1,u_num,Float64(umin),u_grid)
        umax_index = location(CONST_u0,CONST_u1,u_num,Float64(umax),u_grid)
        hmin_index = location(CONST_h0,CONST_h1,h_num,Float64(hmin),h_grid)
        hmax_index = location(CONST_h0,CONST_h1,h_num,Float64(hmax),h_grid)
    elseif typeof(type)==Int64
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end

    f0_3D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    for py in umin_index:umax_index, pz in hmin_index:hmax_index 
        @view(f0_3D_species[:,py,pz]) .= BlackBody_Distribution(PhaseSpace,species,T)
    end
    
    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num_list[species_index]*u_num_list[species_index]*h_num_list[species_index])

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing

end


