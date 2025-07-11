# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct
    BigM::BigMatricesStruct{Matrix{Float32}}
    FluxM::FluxMatricesStruct{Matrix{Float32},Vector{Float32}}

    Implicit::Bool

    M_Bin_Mul_Step::Matrix{Float32}  # temporary array for matrix multiplication of binary terms
    M_Emi_Step::Matrix{Float32}      # temporary array for spatial evaluated emission terms
    Jac::Matrix{Float32}           # jacobian for if Implicit==True
    df::Vector{Float32}                     # change in distribution function
    df_temp::Vector{Float32}                 # change in distribution function
    temp::Matrix{Float32}
    LU64::LinearAlgebra.LU{Float64, Matrix{Float64}, Vector{Int64}}
    LU32::LinearAlgebra.LU{Float32, Matrix{Float32}, Vector{Int64}}

    function EulerStruct(f0::fType,PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatricesStruct{Matrix{Float32}},Flux_Matrices::FluxMatricesStruct{Matrix{Float32},Vector{Float32}},Implicit::Bool)

        self = new()

        self.PhaseSpace = PhaseSpace
        self.BigM = Big_Matrices
        self.FluxM = Flux_Matrices

        self.Implicit = Implicit  

        self.M_Bin_Mul_Step = zeros(Float32,length(f0),length(f0))
        self.M_Emi_Step = zeros(Float32,length(f0),length(f0))
        if Implicit
            self.Jac = zeros(Float32,length(f0),length(f0))
        end
        self.df = zeros(Float32,length(f0))
        self.df_temp = zeros(Float32,length(f0))
        self.temp = zeros(Float32,length(f0),length(f0))
        self.LU32 = lu(zeros(Float32,length(f0),length(f0))+I)
        self.LU64 = lu(zeros(Float64,length(f0),length(f0))+I)

        return self
    end

end