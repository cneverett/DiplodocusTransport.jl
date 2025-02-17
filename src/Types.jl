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
abstract type Cylindrical_RicciType <: ForceType end
struct Cylindrical_Ricci <: Cylindrical_RicciType end
abstract type SynchrotronRadReactType <: ForceType end
struct SynchrotronRadReact <: SynchrotronRadReactType end