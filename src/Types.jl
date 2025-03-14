# Type definitions
abstract type CoordinateType end
abstract type ForceType end
abstract type SteppingMethod <: Function end
abstract type ModeType end

fType = Union{AbstractArray,ArrayPartition}

# ModeType definitions
abstract type IsoType <: ModeType end
struct Iso <: IsoType end
abstract type AxiType <: ModeType end
struct Axi <: AxiType end
abstract type AniType <: ModeType end
struct Ani <: AniType end

# CoordinateType definitions
abstract type SphericalType <: CoordinateType end
struct Spherical <: SphericalType
    mode::ModeType
end

abstract type CylindricalType <: CoordinateType end
struct Cylindrical <: CylindricalType end


# ForceType definitions
abstract type CoordinateForceType <: ForceType end
struct CoordinateForce <: CoordinateForceType end

abstract type SyncRadReactType <: ForceType end
struct SyncRadReact <: SyncRadReactType end