module DiplodocusTransport

    export LoadMatrices, BigMatrices, FluxMatrices
    export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct, OutputStruct, CharacteristicStruct, GridsStruct, ElectroMagneticFieldStruct
    export BinaryStruct, EmiStruct, ForceType
    export BackendType
    export CoordinateType, Cylindrical, Spherical, Cartesian
    export ModeType, Ani, Axi, Iso
    export BoundaryType, Periodic, Open, Closed, Reflective, Escape
    export CoordinateForce, SyncRadReact, GradBInvZDecay
    export BuildBinaryMatrices, BuildEmissionMatrices, BuildFluxMatrices
    export Initialise_Initial_Condition, Location_Species_To_StateVector, Initial_Constant!, Initial_MaxwellJuttner!, Initial_PowerLaw!, Initial_BoostedPowerLaw!, Initial_BlackBody!, Initial_PowerLawExpDecay!
    export Initialise_Injection_Condition, Injection_Constant!, Injection_MaxwellJuttner!, Injection_PowerLaw!, Injection_BoostedPowerLaw!, Injection_BlackBody!, Injection_PowerLawExpDecay!
    export ElectroMagneticField_Constant, ElectroMagneticField_InvZDecay
    export CollisionDomain
    export Solve, ForwardEulerStruct, ForwardSymplecticEulerStruct
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
    using MKL
    using MKLSparse
    using HCubature

    include("Constants.jl")
    include("Types.jl")
    include("Backends.jl")
    include("StructsAndDictionaries.jl")
    include("PhaseSpaceFactors.jl")
    include("MatrixResizing.jl")
    include("NumericalIntegration.jl")
    include("DistributionFunctions.jl")
    include("DistributionMoments.jl")
    include("InitialConditions.jl")
    include("InjectionConditions.jl")
    #include("ValuesOnTheGrid.jl")
    include("DataReading.jl")
    include("GlobalToStateIndices.jl")

    # Fields
    include("ElectroMagneticFieldFunctions.jl")

    # Collisions
    include("Collisions/CollisionDomains.jl")
    include("Collisions/BuildBinaryMatrices.jl")
    include("Collisions/BuildEmissionMatrices.jl")
    include("Collisions/LoadMatrices.jl")
    include("Collisions/EmissionCorrection.jl")

    # Fluxes
    include("Fluxes/BuildFluxMatrices.jl")
    include("Fluxes/FluxFunctionsCoordinate.jl")
    include("Fluxes/FluxFunctionsForces.jl")
    include("Fluxes/FluxFunctionsEmissionForces.jl")

    

    # Stepping 
    include("SteppingStructs.jl")
    include("SteppingMethods.jl")

    # Solver
    include("Solver.jl")

end
