# Type definitions
abstract type CoordinateType end
abstract type ForceType end
abstract type InteractionType end
abstract type SteppingMethodType <: Function end
abstract type ModeType end
abstract type BoundaryType end

# ModeType definitions
#abstract type IsoType <: ModeType end
struct Iso <: ModeType end
#abstract type AxiType <: ModeType end
struct Axi <: ModeType end
#abstract type AniType <: ModeType end
struct Ani <: ModeType end

# BoundaryType definitions
struct Periodic <: BoundaryType end # particles leaving at one coordinate boundary re-appear on the opposite coordinate boundary
struct Open <: BoundaryType end # particles may leave or enter at this flux boundary. Entering may be specified, could be constrained to only photons with other particles closed 
struct Closed <: BoundaryType end # no particles may leave or enter this boundary (no particle reflection when they hit this boundary)
struct Reflective <: BoundaryType end # particles reflect off this boundary (momentum component normal to boundary reverses sign) TODO: this reversal of sign may not be well defined for specific coordinate systems and momentum grids

# CoordinateType Structs: 
@kwdef struct Cartesian <: CoordinateType # x=x, y=y, z=z
    xp_BC::BoundaryType = Periodic()
    xm_BC::BoundaryType = Periodic()
    yp_BC::BoundaryType = Periodic()
    ym_BC::BoundaryType = Periodic()
    zp_BC::BoundaryType = Periodic()
    zm_BC::BoundaryType = Periodic()
end
@kwdef struct Cylindrical <: CoordinateType # x=ρ, y=ϑ, z=z
    α::Float64 = 0.0# != 0 not implemented
    β::Float64 = 0.0
    γ::Float64 = 0.0# != 0 not implemented 
    xp_BC::BoundaryType = Periodic()
    xm_BC::BoundaryType = Periodic()
    yp_BC::BoundaryType = Closed()
    ym_BC::BoundaryType = Closed()
    zp_BC::BoundaryType = Periodic()
    zm_BC::BoundaryType = Periodic() 
end
#Cylindrical() =  Cylindrical(0.0, 0.0, 0.0)
@kwdef struct Spherical <: CoordinateType # x=r, y=θ, z=ψ or momentum space x=p, y=u, z=ϕ
    xp_BC::BoundaryType = Closed()
    xm_BC::BoundaryType = Closed()
    yp_BC::BoundaryType = Closed()
    ym_BC::BoundaryType = Closed()
    zp_BC::BoundaryType = Periodic()
    zm_BC::BoundaryType = Periodic()
end


struct CylindricalMag <: CoordinateType 
    b1::Float64
    b2::Float64
end

# Force Structs: 
struct CoordinateForce <: ForceType end

mutable struct SyncRadReact <: ForceType
    mode::ModeType
    B::Float64
end

struct ExB <: ForceType
    E0::Float64
    B0::Float64
end

# InteractionType definitions
struct BinaryStruct <: InteractionType
    name1::String
    name2::String
    name3::String
    name4::String
end

"""
    EmiStruct

A struct to contain information pertaining to a particular emissive interaction. 
# Fields
- `name1::String`: name of incoming particle (absorbed)
- `name2::String`: name of outgoing particle 1 (emitted)
- `name3::String`: name of outgoing particle 2 (emitted)
- `EmiName::String`: name of the emissive interaction e.g. "Sync" for synchrotron radiation
- `Ext_sampled::Vector{Float64}`: vector of sampled values of the external parameter for which the emissive interaction matrices are built e.g. sampled magnetic field values for synchrotron radiation
- `mode::ModeType`: mode of the emissive interaction e.g. "Iso", "Axi" or "Ani"
- `Force::Bool`: if name1==name2 the emissive interaction may be treated as reactive force e.g. radiation reaction for synchrotron. In which case this force can be applied directly, within the emission matrix or separately as a regular force.
- `Domain::Union{Vector{Int64},Nothing}`: vector of spatial indices (offset_space) referring to spatial sub-domains in which to apply the emissive interaction. If `nothing` the interaction will be applied to all spatial sub-domains.
"""
@kwdef struct EmiStruct <: InteractionType
    name1::String
    name2::String
    name3::String
    EmiName::String
    Ext_sampled::Vector{Float64}
    mode::ModeType
    Force::Bool = true
    Domain::Union{Vector{Int64},Nothing} = nothing
end

# ElectromagneticField Structs 
abstract type ElectroMagneticFieldType end

@kwdef struct Constant_ElectromagneticField <: ElectroMagneticFieldType
    # Define the constant magnetic field with strength B in the z-direction
    parameters::Vector{Float64} = [1e-5,0.0] # B (Tesla), E (Tesla*c)
    ElectroMagneticFieldFunction::Function = Constant_ElectroMagneticFieldFunction
end