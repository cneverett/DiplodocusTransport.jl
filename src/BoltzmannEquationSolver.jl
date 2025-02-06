module BoltzmannEquationSolver

export Solver, InitialConditions, LoadMatrices, Big_Matrices

using JLD2
import BoltzmannCollisionIntegral as BCI
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
include("AllocateMatrices.jl")
include("LoadMatrices.jl")
include("DistributionFunctions.jl")
include("DistributionMoments.jl")
include("InitialConditions.jl")
include("ValuesOnTheGrid.jl")
include("Solver.jl")

# post-processing
include("Plotting.jl")

end
