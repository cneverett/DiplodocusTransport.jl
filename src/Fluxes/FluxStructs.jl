"""
    FluxMatricesStruct()

A struct for storing the flux matrices associated with the simulation.
"""
struct FluxMatricesStruct{T<:Union{Float32,Float64}}

    Ap_Flux::Vector{T}                  # Forward Boundary time flux 
    Am_Flux::Vector{T}                  # Backward Boundary time flux
    
    X_Flux::SparseMatrixCSC{T,Int64}    # sum of space fluxes 
    P_Flux::SparseMatrixCSC{T,Int64}    # sum of momentum fluxes 

    Vol::Vector{T}                      # SpaceTime volume element

    B_Flux::SparseMatrixCSC{T,Int64}    # B Flux through x boundaries
    C_Flux::SparseMatrixCSC{T,Int64}    # C Flux through y boundaries
    D_Flux::SparseMatrixCSC{T,Int64}    # D Flux through z boundaries

    I_Flux::SparseMatrixCSC{T,Int64}    # I Flux through px boundaries
    J_Flux::SparseMatrixCSC{T,Int64}    # J Flux through py boundaries
    K_Flux::SparseMatrixCSC{T,Int64}    # K Flux through pz boundaries

end

#========= Force =========#
#=========================#

abstract type AbstractForce end
abstract type SpaceVectorForce <: AbstractForce end
abstract type SpaceScalarForce <: AbstractForce end
abstract type AnalyticForce <: AbstractForce end

struct CoordinateForce <: AbstractForce end

mutable struct SyncRadReact <: AnalyticForce
    mode::AbstractMode
    B::Float64
end

struct FirstOrderGuidingCentre <: SpaceVectorForce end

struct ExB <: AbstractForce
    E0::Float64
    B0::Float64
end

struct GradBInvZDecay <: AbstractForce
end