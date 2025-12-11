module DiplodocusTransport

    export LoadMatrices, BigMatrices, FluxMatrices
    export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct, OutputStruct, CharacteristicStruct, GridsStruct
    export BinaryStruct, EmiStruct, ForceType
    export Cylindrical, Spherical, Cartesian, Ani, Axi, Iso
    export Periodic, Open, Closed, Reflective
    export CoordinateForce, SyncRadReact
    export BuildBigMatrices, BuildFluxMatrices
    export Initialise_Initial_Condition, Location_Species_To_StateVector, Initial_Constant!, Initial_MaxwellJuttner!, Initial_PowerLaw!, Initial_UnBoostedPowerLaw!, Initial_BlackBody!
    export Initialise_Injection_Condition, Injection_Constant!, Injection_MaxwellJuttner!, Injection_PowerLaw!, Injection_UnBoostedPowerLaw!, Injection_BlackBody!
    export Solve, EulerStruct
    export SolutionFileLoad
    export MaxwellJuttner_Distribution
    export CUDABackend, CPUBackend

    using JLD2
    using DiplodocusCollisions: bounds, location, deltaVector, meanVector, deltaEVector, EmissionFileName, BinaryFileName, EmissionFileLoad_Matrix, BinaryFileLoad_Matrix, DoesConserve
    using LinearAlgebra
    using ProgressMeter
    using Bessels
    using Statistics
    using SparseArrays
    using CUDA
    using CUDA.CUSPARSE

    include("Constants.jl")
    include("Types.jl")
    include("Backends.jl")
    include("StructsAndDictionaries.jl")
    include("PhaseSpaceFactors.jl")
    include("MatrixResizing.jl")
    include("DistributionFunctions.jl")
    include("DistributionMoments.jl")
    include("InitialConditions.jl")
    include("InjectionConditions.jl")
    #include("ValuesOnTheGrid.jl")
    include("DataReading.jl")
    include("GlobalToStateIndices.jl")

    # Collisions
    include("Collisions/BuildBigMatrices.jl")
    include("Collisions/LoadMatrices.jl")
    include("Collisions/EmissionCorrection.jl")

    # Fluxes
    include("Fluxes/BuildFluxMatrices.jl")
    include("Fluxes/FluxFunctionsCoordinate.jl")
    include("Fluxes/FluxFunctionsForces.jl")

    # Stepping 
    include("SteppingStructs.jl")
    include("SteppingMethods.jl")

    # Solver
    include("Solver.jl")

end
