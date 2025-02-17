module BoltzmannEquationSolver

export Solver, InitialConditions, LoadMatrices, BigMatrices, SpaceTimeStruct, FluxMatricesCoordinate, FluxMatricesForce

using JLD2
import BoltzmannCollisionIntegral as BCI
using LinearAlgebra
using ProgressMeter
#using LoopVectorization
#using Tullio
#using TensorOperations
using Bessels
using Statistics
using RecursiveArrayTools
using CairoMakie

include("Types.jl")
include("StructsAndDictionaries.jl")
include("PhaseSpaceFactors.jl")
include("MatrixResizing.jl")
include("AllocateMatrices.jl")
include("LoadMatrices.jl")
include("DistributionFunctions.jl")
include("DistributionMoments.jl")
include("InitialConditions.jl")
include("ValuesOnTheGrid.jl")

# Fluxes
include("Fluxes/FluxStructs.jl")
include("Fluxes/FluxMatrices.jl")
include("Fluxes/FluxFunctionsCoordinate.jl")

# Stepping 
include("SteppingStructs.jl")
include("SteppingMethods.jl")

# Solver
include("Solver.jl")

# post-processing
include("Plotting.jl")

end
