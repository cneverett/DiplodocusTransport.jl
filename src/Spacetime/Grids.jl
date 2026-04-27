abstract type AbstractSpacetimeGrid end

struct UniformGrid <: AbstractSpacetimeGrid
    up::Float64
    low::Float64
    num::Int64
end
struct Log10Grid <: AbstractSpacetimeGrid
    up::Float64
    low::Float64
    num::Int64
end
struct StretchGrid <: AbstractSpacetimeGrid
    low::Float64
    a::Float64
    ratio::Float64
    num::Int64
end

"""
    SpacetimeGrid!(pos,g,::AbstractMetric,::AbstractCoordinates,)

Returns a `num+1` long `Vector{Float}` of grid bounds for a given spacetime coordinate based on the type of `grid` specified.
"""
SpacetimeGrid!(grid::AbstractSpacetimeGrid) = error("Spacetime grid function not defined for this grid $(typeof(grid)).")


function SpacetimeGrid(grid::UniformGrid) 
    #= uniform spacing
        | a | a | a | a | a | 
        low                  up
        a = (up - low)/num
    =#
    return [range(grid.low,grid.up,grid.num+1);]
end

function SpacetimeGrid(grid::Log10Grid) 
    #= log10 spacing
        | 1 | 10 | 100 | 1000 | 10000 | 
        10^low                        10^up
    =#
    return 10 .^[range(grid.low,grid.up,grid.num+1);]
end

function SpacetimeGrid(grid::StretchGrid) 
    #= stretch spacing
        |       |      |           |              |                   |                   | 
       low    (low+a) (low+a)*r   (low+a)*r^2   (low+a)*r^3         (low+a)*r^4      (low+a)*r^5         
    =#
    @assert grid.ratio != 1.0 "Ratio must be different from 1 for stretch grid"
    vec = [(grid.a-grid.low)*grid.ratio^i for i in 0:(grid.num-1)]
    return prepend!(vec, grid.low)
end