"""
    CoordinateTransform(pos,from::AbstractCoordinates,to::AbstractCoordinates)

Returns the coordinates ``pos=[t,x,y,z]`` in the ``to`` coordinate system given the coordinates in the ``from`` coordinate system.
"""
CoordinateTransform(pos,from::AbstractCoordinates,to::AbstractCoordinates) = error("Coordinate transform function not defined for coordinates $(typeof(from)) to coordinates $(typeof(to)).")


function CoordinateTransform(pos::AbstractVector{T},::Paraboloidal,::Cartesian) where T

    t = pos[1]
    ϕ = pos[2]
    u = pos[3]
    v = pos[4]

    x::T = u*v*cos(ϕ)
    y::T = u*v*sin(ϕ)
    z::T = 1/2 * (u^2 - v^2)

    return SVector{4,T}(t,x,y,z)

end

function CoordinateTransform(pos::AbstractVector{T},::Paraboloidal,::Cylindrical) where T

    t = pos[1]
    ϕ = pos[2]
    u = pos[3]
    v = pos[4]

    ρ::T = u*v
    z::T = 1/2 * (u^2 - v^2)

    return SVector{4,T}(t,ρ,ϕ,z)

end