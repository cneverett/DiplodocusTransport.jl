"""
    Initialise_Initial_Condition(PhaseSpace)

Returns a zero vector with elements for the initial conditions of all particles at all positions in phase space. 
"""
function Initialise_Initial_Condition(PhaseSpace::PhaseSpaceStruct)

    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num
    
    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*x_num*y_num*z_num

    return zeros(Float32,m)

end

function Location_Species_To_StateVector(Initial::Vector{Float32},PhaseSpace::PhaseSpaceStruct;x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1,species_index::Int64=0)

    if iszero(species_index)
        error("Species not defined")
    end
    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num

    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list

    px_num = px_num_list[species_index]
    py_num = py_num_list[species_index]
    pz_num = pz_num_list[species_index]

    name_list = PhaseSpace.name_list

    off_space = off_space = (x_idx-1)*y_num*z_num+(y_idx-1)*z_num+z_idx-1
    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end
    off_name = offset[species_index]

    start_idx = 1+off_name+off_space
    end_idx = pz_num*px_num*py_num+off_name+off_space

    return @view(Initial[start_idx:end_idx])

end


"""
    Initial_PowerLaw(Lists,species,pmin,pmax,umin,uman,hmin,hmax,index,num_Init)

A power-law distribution is typically defined by N(E) ∝ E^(-index). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index) * p/E = pE^(-index-1). Averaging this over a cell gives f(p)_avg = [E^(1-index)/(1-index)]/[p] where [] denote evaluation at the cell bounds.
"""
function Initial_PowerLaw(PhaseSpace::PhaseSpaceStruct,species::String,pmin::S,pmax::S,umin::S,umax::S,hmin::S,hmax::S,index::AbstractFloat,num_Init::AbstractFloat) where S <: Union{Float32,Float64,Int64}

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

function Initial_PowerLaw!(Initial::Vector{Float32},PhaseSpace::PhaseSpaceStruct,species::String,pmin::S,pmax::S,umin::S,umax::S,hmin::S,hmax::S,index::AbstractFloat,num_Init::AbstractFloat) where S <: Union{Float32,Float64,Int64}

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

    Initial_local = Location_Species_To_StateVector(Initial,PhaseSpace,species_index=species_index)

    Initial_local .+= f0_species

    return nothing
end

"""
    Initial_UnBoostedPowerLaw(PhaseSpace,species,pmin,pmax,Gamma,index,num_Init)

Takes an isotropic power-law distribution, with minimum momentum `pmin`, maximum momentum `pmax` and `index` in some frame propagating with Lorentz factor `Gamma` in the z-direction and returns the distribution as observed by a static observer, with a number density of `num_Init`.
"""
function Initial_UnBoostedPowerLaw(PhaseSpace::PhaseSpaceStruct,species::String,pmin::S,pmax::S,Gamma::S,index::AbstractFloat,num_Init::AbstractFloat) where S <: Union{Float32,Float64}

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
    mass = getfield(DC,Symbol("mu"*name_list[species_index]))
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
``
    # Set Rapidity 
    w::Float64 = acosh(Gamma)
    cw = cosh(w)
    sw = sinh(w)
    # de-boost pmin and pmax 
    Emin::Float64 = sqrt(mass^2+pmin^2)
    Emax::Float64 = sqrt(mass^2+pmax^2)
    pmin_UB::Float64 = sqrt((cosh(w)*Emin-sinh(w)*pmin)^2-mass^2)
    pmax_UB::Float64 = sqrt((cosh(w)*Emax-sinh(w)*pmax)^2-mass^2)
    Emin_UB::Float64 = sqrt(mass^2+pmin_UB^2)
    Emax_UB::Float64 = sqrt(mass^2+pmax_UB^2)


    pmin_index = DC.location(pl,pu,p_num,pmin_UB,p_grid)
    pmax_index = DC.location(pl,pu,p_num,pmax_UB,p_grid)

    # power law averaged over cell width.
    for px in pmin_index:pmax_index, py in 1:u_num, pz in 1:h_num 
        # f calculated using simple trapezium rule (fabc is f at the a bound of p, b bound of u and c bound of h)
        pp = p_r[px+1]
        pm = p_r[px]
        up = u_r[py+1]
        um = u_r[py]
        hp = h_r[pz+1]
        hm = h_r[pz]
        fmmm = pm^2 * (cw*sqrt(mass^2+pm^2)-sw*pm*um)^(-index) / sqrt(mass^2+pm^2) / sqrt((cw*sqrt(mass^2+pm^2)-sw*pm*um)^2-mass^2)
        fpmm = pp^2 * (cw*sqrt(mass^2+pp^2)-sw*pp*um)^(-index) / sqrt(mass^2+pp^2) / sqrt((cw*sqrt(mass^2+pp^2)-sw*pp*um)^2-mass^2)
        fmpm = pm^2 * (cw*sqrt(mass^2+pm^2)-sw*pm*up)^(-index) / sqrt(mass^2+pm^2) / sqrt((cw*sqrt(mass^2+pm^2)-sw*pm*up)^2-mass^2)
        fmmp = pm^2 * (cw*sqrt(mass^2+pm^2)-sw*pm*um)^(-index) / sqrt(mass^2+pm^2) / sqrt((cw*sqrt(mass^2+pm^2)-sw*pm*um)^2-mass^2)
        fppm = pp^2 * (cw*sqrt(mass^2+pp^2)-sw*pp*up)^(-index) / sqrt(mass^2+pp^2) / sqrt((cw*sqrt(mass^2+pp^2)-sw*pp*up)^2-mass^2)
        fpmp = pp^2 * (cw*sqrt(mass^2+pp^2)-sw*pp*um)^(-index) / sqrt(mass^2+pp^2) / sqrt((cw*sqrt(mass^2+pp^2)-sw*pp*um)^2-mass^2)
        fmpp = pm^2 * (cw*sqrt(mass^2+pm^2)-sw*pm*up)^(-index) / sqrt(mass^2+pm^2) / sqrt((cw*sqrt(mass^2+pm^2)-sw*pm*up)^2-mass^2)
        fppp = pp^2 * (cw*sqrt(mass^2+pp^2)-sw*pp*up)^(-index) / sqrt(mass^2+pp^2) / sqrt((cw*sqrt(mass^2+pp^2)-sw*pp*up)^2-mass^2)
        f0_3D_species[px,py,pz] = (fmmm + fpmm + fmpm + fmmp + fppm + fpmp + fmpp + fppp) / 8
        if px == pmin_index
            f0_3D_species[px,py,pz] *= (p_r[px+1]-pmin_UB) / (p_r[px+1]-p_r[px])
        elseif px == pmax_index
            f0_3D_species[px,py,pz] *= (pmax_UB-p_r[px]) / (p_r[px+1]-p_r[px])
        end
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

