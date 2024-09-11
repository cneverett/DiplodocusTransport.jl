module DiscreteBoltzmannEquation

export BoltzmannEquationSolver, InitialConditions, LoadMatricies

using JLD2
using BoltzmannCollisionIntegral
using CairoMakie
using DifferentialEquations
using LoopVectorization
using Bessels
using Statistics

#include("Setup.jl")
include("StructsAndDictionaries.jl")
include("MatrixResizing.jl")
include("LoadMatricies.jl")
include("DistributionFunctions.jl")
include("InitialConditions.jl")
include("Solver.jl")

# post-processing
include("DistributionMoments.jl")
include("Plotting.jl")

end
