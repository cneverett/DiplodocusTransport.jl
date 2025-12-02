abstract type BackendType end

struct CPUBackend <: BackendType
end
struct CUDABackend <: BackendType 
end


# zeros functions
function Base.zeros(::CPUBackend, T::Type{<:AbstractFloat}, dims::Int...)
    return zeros(T, dims...)
end
function Base.zeros(::CUDABackend, T::Type{<:AbstractFloat}, dims::Int...)
    @assert T == Float32 "CUDA backend only supports Float32 type."
    return CUDA.zeros(T, dims...)
end