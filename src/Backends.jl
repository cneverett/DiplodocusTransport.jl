abstract type AbstractBackend end

struct CPUBackend <: AbstractBackend end
struct CUDABackend <: AbstractBackend end

# zeros functions
function Base.zeros(::CPUBackend, T::Type{<:AbstractFloat}, dims::Int...)
    return zeros(T, dims...)
end
function Base.zeros(::CUDABackend, T::Type{<:AbstractFloat}, dims::Int...)
    #@assert T == Float32 "CUDA backend only supports Float32 type."
    return CUDA.zeros(T, dims...)
end

# ones functions
function Base.ones(::CPUBackend, T::Type{<:AbstractFloat}, dims::Int...)
    return ones(T, dims...)
end
function Base.ones(::CUDABackend, T::Type{<:AbstractFloat}, dims::Int...)
    #@assert T == Float32 "CUDA backend only supports Float32 type."
    return CUDA.ones(T, dims...)
end