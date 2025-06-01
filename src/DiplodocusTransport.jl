module DiplodocusTransport

export Solver, InitialConditions, LoadMatrices, BigMatrices, FluxMatrices
export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct 

using JLD2
import DiplodocusCollisions as DC
using LinearAlgebra
using ProgressMeter
using Bessels
using Statistics
using RecursiveArrayTools
using CairoMakie

include("Types.jl")
include("StructsAndDictionaries.jl")
include("PhaseSpaceFactors.jl")
include("MatrixResizing.jl")
include("DistributionFunctions.jl")
include("DistributionMoments.jl")
include("InitialConditions.jl")
include("ValuesOnTheGrid.jl")

# Collisions
include("Collisions/BuildBigMatrices.jl")
include("Collisions/LoadMatrices.jl")

# Fluxes
include("Fluxes/BuildFluxMatrices.jl")
include("Fluxes/FluxFunctionsCoordinate.jl")
include("Fluxes/FluxFunctionsForces.jl")

# Stepping 
include("SteppingStructs.jl")
include("SteppingMethods.jl")

# Solver
include("Solver.jl")

# post-processing
include("Plotting.jl")

end
