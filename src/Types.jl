# Type definitions
abstract type CoordinateType end
abstract type ForceType end
abstract type SteppingMethod <: Function end
abstract type ModeType end

fType = Union{AbstractArray,ArrayPartition}

# ModeType definitions
#abstract type IsoType <: ModeType end
struct Iso <: ModeType end
#abstract type AxiType <: ModeType end
struct Axi <: ModeType end
#abstract type AniType <: ModeType end
struct Ani <: ModeType end

# CoordinateType definitions
abstract type SphericalType <: CoordinateType end
struct Spherical <: SphericalType
end

abstract type CylindricalType <: CoordinateType end
struct Cylindrical <: CylindricalType end

# ForceType definitions
abstract type CoordinateForceType <: ForceType end
struct CoordinateForce <: CoordinateForceType end

abstract type SyncRadReactType <: ForceType end
struct SyncRadReact <: SyncRadReactType
    mode::ModeType
end

# InteractionType definitions
struct BinaryStruct
    name1::String
    name2::String
    name3::String
    name4::String
end
struct EmiStruct
    name1::String
    name2::String
    name3::String
    EmiName::String
    mode::ModeType
end