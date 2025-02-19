# Type definitions
abstract type CoordinateType end
abstract type ForceType end
abstract type SteppingMethod <: Function end

fType = Union{AbstractArray,ArrayPartition}

# CoordinateType definitions
abstract type SphericalType <: CoordinateType end
struct Spherical <: SphericalType end
abstract type CylindricalType <: CoordinateType end
struct Cylindrical <: CylindricalType end
# ForceType definitions
abstract type CoordinateForceType <: ForceType end
struct CoordinateForce <: CoordinateForceType end

abstract type SyncRadReactType <: ForceType end
struct SyncRadReact <: SyncRadReactType end