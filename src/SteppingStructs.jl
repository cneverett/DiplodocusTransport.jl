# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct
    BigM::BigMatricesStruct{Matrix{T}}
    FluxM::FluxMatricesStruct{Matrix{T},Vector{T}}

    Implicit::Bool

    M_Bin_Mul_Step::Matrix{T}  # temporary array for matrix multiplication of binary terms
    M_Emi_Step::Matrix{T}      # temporary array for spatial evaluated emission terms
    Jac::Matrix{T}           # jacobian for if Implicit==True
    df::Vector{T}                     # change in distribution function
    df_temp::Vector{T}                 # change in distribution function
    temp::Matrix{T}
    LU::LinearAlgebra.LU{T, Matrix{T}, Vector{Int64}}

    function EulerStruct{T}(f0::fType,PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatricesStruct{AbstractMatrix{T}},Flux_Matrices::FluxMatricesStruct{AbstractMatrix{T},AbstractVector{T}},Implicit::Bool) where T <: AbstractFloat

        self = new()

        self.PhaseSpace = PhaseSpace
        self.BigM = Big_Matrices
        self.FluxM = Flux_Matrices

        self.Implicit = Implicit  

        self.M_Bin_Mul_Step = zeros(T,length(f0),length(f0))
        self.M_Emi_Step = zeros(T,length(f0),length(f0))
        if Implicit
            self.Jac = zeros(T,length(f0),length(f0))
        end
        self.df = zeros(T,length(f0))
        self.df_temp = zeros(T,length(f0))
        self.temp = zeros(T,length(f0),length(f0))
        self.LU = lu(zeros(T,length(f0),length(f0))+I)

        return self
    end

end