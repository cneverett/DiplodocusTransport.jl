"""
    Initial_PowerLaw(Lists,species,pmin,pmax,index,num_Init)

A power-law distribution is typically defined by N(E) ∝ E^(-index). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index) * p/E = pE^(-index-1). Averaging this over a cell gives f(p)_avg = [E^(1-index)/(1-index)]/[p] where [] denote evaluation at the cell bounds.
"""
function Initial_PowerLaw(PhaseSpace::PhaseSpaceStruct,species::String,pmin::S,pmax::S,umin::S,umax::S,hmin::S,hmax::S,index::Float32,num_Init::AbstractFloat) where S <: Union{Float32,Float64,Int64}

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

    species_index = findfirst(==(species),name_list)
    pr = pr_list[species_index]
    dp = dp_list[species_index]
    du = du_list[species_index]

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
    #pr = DC.bounds(pl,pu,p_num,p_grid)
    #dp = DC.deltaVector(pr)
    #du = DC.deltaVector(DC.bounds(DC.u_low,DC.u_up,u_num,u_grid))
    mass = getfield(DC,Symbol("mu"*name_list[species_index]))

    type = zero(S)
    if typeof(type)==Float32
        pmin_index = DC.location(pl,pu,p_num,Float64(pmin),p_grid)
        pmax_index = DC.location(pl,pu,p_num,Float64(pmax),p_grid)
        umin_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umin),u_grid)
        umax_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umax),u_grid)
        hmin_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmin),h_grid)
        hmax_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmax),h_grid)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end
    
    #println(pmin_index,pmax_index,umin_index,umax_index)

    # power law averaged over cell width.
    for px in pmin_index:pmax_index, py in umin_index:umax_index, pz in hmin_index:hmax_index 
        f0_3D_species[px,py,pz] = sqrt(mass^2+pr[px+1]^2)^(1-index) - sqrt(mass^2+pr[px]^2)^(1-index)
        f0_3D_species[px,py,pz] /= 1-index
        #f0_3D_species[px,py,pz] /= (pr[px+1]-pr[px])
    end
    # set values and normlaise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    # scale by dp*du 
    #=for i in axes(u0_2D_species,1), j in axes(u0_2D_species,2)
        u0_2D_species[i,j] *= dp[i] * du[j]
    end=#

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    return Float32.(f0_species)
end

"""
    Initial_Constant(PhaseSpace,species,pmin,pmax,umin,umax,hmin,hmax,num_Init)

Divides the initial number density `num_Init` equally among momentum-space bins in the range of `pmin` to `pmax`, `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function Initial_Constant(PhaseSpace::PhaseSpaceStruct,species::String,pmin::T,pmax::T,umin::T,umax::T,hmin::T,hmax::T,num_Init::AbstractFloat) where T <: Union{Float32,Float64,Int64}

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
        pmin_index = DC.location(pl,pu,p_num,Float64(pmin),p_grid)
        pmax_index = DC.location(pl,pu,p_num,Float64(pmax),p_grid)
        umin_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umin),u_grid)
        umax_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umax),u_grid)
        hmin_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmin),h_grid)
        hmax_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmax),h_grid)
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

    return Float32.(f0_species)
end

"""
    Initial_MaxwellJuttner(PhaseSpace,species,T,umin,umax,hmin,hmax,num_Init)

Returns initial conditions for `species` as a Maxwell-Juttner distribution of temperature `T` in Kelvin with a number density of `num_Init` and angular range `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function Initial_MaxwellJuttner(PhaseSpace::PhaseSpaceStruct,species::String,T::Float64,umin::S,umax::S,hmin::S,hmax::S,num_Init::AbstractFloat)  where S <: Union{Float32,Float64,Int64}

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
        umin_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umin),u_grid)
        umax_index = DC.location(DC.u_low,DC.u_up,u_num,Float64(umax),u_grid)
        hmin_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmin),h_grid)
        hmax_index = DC.location(DC.h_low,DC.h_up,h_num,Float64(hmax),h_grid)
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

    return Float32.(f0_species)

end

