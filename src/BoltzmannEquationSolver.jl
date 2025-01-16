module BoltzmannEquationSolver

export Solver, InitialConditions, LoadMatrices

using JLD2
using BoltzmannCollisionIntegral
using LinearAlgebra
using CairoMakie
using DifferentialEquations
using ProgressLogging
#using LoopVectorization
using Tullio
#using TensorOperations
using Bessels
using Statistics
using RecursiveArrayTools

#include("Setup.jl")
include("StructsAndDictionaries.jl")
include("PhaseSpaceFactors.jl")
include("MatrixResizing.jl")
include("LoadMatrices.jl")
include("DistributionFunctions.jl")
include("DistributionMoments.jl")
include("InitialConditions.jl")
include("ValuesOnTheGrid.jl")
include("Solver.jl")

# post-processing
include("Plotting.jl")

end
