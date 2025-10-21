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
struct Open <: BoundaryType end # particles may leave or enter at this flux boundary. Entering may be specified (to be implemented), could be constrained to only photons with other particles closed (to be implemented)
struct Closed <: BoundaryType end # no particle may leave or enter this boundary

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
@kwdef struct Spherical <: CoordinateType # x=r, y=θ, z=ψ
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

# ForceType Structs: 
struct CoordinateForce <: ForceType end
struct SyncRadReact <: ForceType
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
struct EmiStruct <: InteractionType
    name1::String
    name2::String
    name3::String
    EmiName::String
    Ext::Vector{Float64}
    mode::ModeType
end