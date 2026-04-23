module DiplodocusTransport

    # Metrics
    export AbstractMetric, MetricComponents!, ChristoffelComponents!, Minkowski, Schwarzschild, Kerr
    # Coordinates
    export AbstractCoordinates, Cartesian, Cylindrical, Spherical, ModifiedSpherical, Paraboloidal
    # Boundary Conditions
    export AbstractBoundaryCondition, Periodic, Open, Closed, Reflective, Escape
    # Modes
    export AbstractMode, Iso, Axi, Ani
    # Spacetime Grids 
    export AbstractSpacetimeGrid, SpacetimeGrid, UniformGrid, Log10Grid, StretchGrid
    # Tetrads
    export AbstractTetrad, TetradComponents!, InverseTetradComponents!, RicciRotationComponents!, StationaryObserverTetrad, UniformElectromagneticFieldTetrad, ParabolicForceFreeFieldTetrad
    # Backends
    export AbstractBackend, CPUBackend, CUDABackend
    # Phasespace
    export PhaseSpaceStruct, MomentumStruct, SpacetimeStruct, TimeStruct, CharacteristicStruct, GridsStruct
    # Forces 
    export AbstractForce, CoordinateForce, AnalyticForce, SpaceVectorForce, SpaceScalarForce
    export FirstOrderGuidingCentre, SyncRadReact, GradBInvZDecay
    # Fluxes
    export BuildFluxMatrices, VolFunction
    # Collisions
    export AbstractInteraction, BinaryInteraction, EmissiveInteraction, InteractionDomain
    export LoadMatrices, BuildBinaryMatrices, BuildEmissionMatrices
    # Initial and Injection Conditions
    export InitialiseInitialCondition, InitialConstant!, InitialMaxwellJuttner!, InitialPowerLaw!, InitialBoostedPowerLaw!, InitialBlackBody!, InitialPowerLawExpDecay!
    export InitialiseInjectionCondition, InjectionConstant!, InjectionMaxwellJuttner!, InjectionPowerLaw!, InjectionBoostedPowerLaw!, InjectionBlackBody!, InjectionPowerLawExpDecay!
    # Solvers 
    export AbstractSteppingMethod, Solve, ForwardEulerStruct, ForwardSymplecticEulerStruct
    export OutputStruct, SolutionFileLoad
    # Utilities
    export GlobalIndicesToStateIndex,LocationSpeciesToStateVector
    # Distribution Moments 
    export FourFlow, HydroFourVelocity, HydroThreeVelocity, ProjectionTensor, StressEnergyTensor, ScalarNumberDensity, ScalarMassDensity, ScalarEnergyDensity, ScalarPressure, ScalarTemperature

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
    using StaticArrays
    using ForwardDiff

    include("Constants.jl")
    include("Backends.jl")
    include("NumericalIntegration.jl")

    # Phasespace
    include("Spacetime/CoordinateStructs.jl")
    include("Spacetime/Grids.jl")
    include("Spacetime/Metrics.jl")
    include("Spacetime/Christoffel.jl")
    include("Spacetime/Tetrads.jl")
    include("Spacetime/RicciRotation.jl")

    
    #include("Types.jl")
    include("PhaseSpaceStructs.jl")
    include("ElectromagneticFields.jl")
    #include("MatrixResizing.jl")
    include("DistributionFunctions.jl")
    include("DistributionMoments.jl")
    include("InitialConditions.jl")
    include("InjectionConditions.jl")
    include("GlobalToStateIndices.jl")

    # Fluxes
    include("Fluxes/FluxStructs.jl")
    include("Fluxes/BuildFluxMatrices.jl")
    include("Fluxes/FluxFunctionsCoordinate.jl")
    include("Fluxes/FluxFunctionsForces.jl")
    include("Fluxes/FluxFunctionsEmissionForces.jl") 

    # Collisions
    include("Collisions/InteractionStructs.jl")
    include("Collisions/InteractionDomains.jl")
    include("Collisions/InteractionPhaseSpaceFactors.jl")
    include("Collisions/BuildBinaryMatrices.jl")
    include("Collisions/BuildEmissionMatrices.jl")
    include("Collisions/LoadMatrices.jl")
    include("Collisions/EmissionCorrection.jl")

    # Solvers
    include("Solvers/SteppingStructs.jl ")
    include("Solvers/SteppingMethods.jl")
    include("Solvers/Solver.jl")
    include("DataReading.jl")

end
