# Type definitions
abstract type CoordinateType end
abstract type ForceType end
abstract type InteractionType end
abstract type SteppingMethodType <: Function end
abstract type ModeType end

fType = Union{AbstractArray,ArrayPartition}

# ModeType definitions
#abstract type IsoType <: ModeType end
struct Iso <: ModeType end
#abstract type AxiType <: ModeType end
struct Axi <: ModeType end
#abstract type AniType <: ModeType end
struct Ani <: ModeType end

# CoordinateType Structs: 
struct Spherical <: CoordinateType end
struct Cylindrical <: CoordinateType end

# ForceType Structs: 
struct CoordinateForce <: ForceType end
struct SyncRadReact <: ForceType
    mode::ModeType
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
    mode::ModeType
end