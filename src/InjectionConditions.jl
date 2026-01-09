"""
    Initialise_Injection_Condition(PhaseSpace)

Returns a zero vector with elements for the injection conditions of all particles at all positions in phase space. 
"""
function Initialise_Injection_Condition(PhaseSpace::PhaseSpaceStruct)

    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num
    
    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*x_num*y_num*z_num

    return zeros(Float64,m)

end


"""
    Injection_PowerLaw!(Injection,PhaseSpace,species,pmin,pmax,umin,uman,hmin,hmax,index,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a power law distribution generated using `Initial_PowerLaw!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function Injection_PowerLaw!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::Float64,pmax::Float64,umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,index::AbstractFloat=2.0,num_Inj::AbstractFloat=1.0,rate_Inj::Float64=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where F<:AbstractFloat
    
    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    Initial_PowerLaw!(tmp,PhaseSpace,species;pmin=pmin,pmax=pmax,umin=umin,umax=umax,hmin=hmin,hmax=hmax,index=index,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    tr = PhaseSpace.Grids.tr
    dt0 = tr[2] - tr[1]
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * (tmp/dt0)

    return nothing
end

"""
    Injection_BoostedPowerLaw!(Injection,PhaseSpace,species,pmin,pmax,Gamma,index,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a power law distribution generated using `Initial_BoostedPowerLaw!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function Injection_BoostedPowerLaw!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::Float64,pmax::Float64,Gamma::Float64,index::AbstractFloat,num_Inj::Float64=1.0,rate_Inj::Float64=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where F<:AbstractFloat

     # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    Initial_UnBoostedPowerLaw!(tmp,PhaseSpace,species;pmin=pmin,pmax=pmax,Gamma=Gamma,index=index,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    tr = PhaseSpace.Grids.tr
    dt0 = tr[2] - tr[1]
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * (tmp/dt0)

    return nothing

end

"""
    Injection_Constant!(Injection,PhaseSpace,species,pmin,pmax,umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a constant distribution generated using `Initial_Constant!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function Injection_Constant!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::T,pmax::T,umin::T,umax::T,hmin::T,hmax::T,num_Inj::AbstractFloat=1.0,rate_Inj::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where T <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    Initial_Constant!(tmp,PhaseSpace,species;pmin=pmin,pmax=pmax,umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    tr = PhaseSpace.Grids.tr
    dt0 = tr[2] - tr[1]
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * (tmp/dt0)

    return nothing
end

"""
    Injection_MaxwellJuttner!(Injection,PhaseSpace,species,T,umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a Maxwell-Juttner distribution generated using `Initial_MaxwellJuttner!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function Injection_MaxwellJuttner!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;T::Float64,umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::Float64=1.0,rate_Inj::Float64=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where F<:AbstractFloat

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    Initial_MaxwellJuttner!(tmp,PhaseSpace,species;T=T,umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    tr = PhaseSpace.Grids.tr
    dt0 = tr[2] - tr[1]
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * (tmp/dt0)

    return nothing
end

"""
    Injection_BlackBody!(Injection,PhaseSpace,species,T,umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a Black-Body distribution generated using `Initial_BlackBody!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function Injection_BlackBody!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;T::Float64,umin::S,umax::S,hmin::S,hmax::S,num_Inj::AbstractFloat=1.0,rate_Inj::AbstractFloat=1.0,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1)  where S <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    Initial_BlackBody!(tmp,PhaseSpace,species;T=T,umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)

    tr = PhaseSpace.Grids.tr
    dt0 = tr[2] - tr[1]
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * (tmp/dt0)

    return nothing
end