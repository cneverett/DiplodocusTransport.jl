# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct <: SteppingMethod

    PhaseSpace::PhaseSpaceStruct
    BigM::BigMatricesStruct
    FluxM::FluxMatricesStruct

    Implicit::Bool

    M_Bin_Mul_Step::Array{Float32,2}  # temporary array for matrix multiplication of binary terms
    M_Emi_Step::Array{Float32,2}      # temporary array for spatial evaluated emission terms
    Jac::Array{Float32,2}           # jacobian for if Implicit==True
    df::fType                       # change in distribution function
    df_temp::fType                  # change in distribution function
    temp::AbstractArray{Float32,2}
    L::AbstractArray{Float32,2}
    U::AbstractArray{Float32,2}

    function EulerStruct(f0,PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatricesStruct,Flux_Matrices::FluxMatricesStruct,Implicit::Bool)

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
        self.df = fill!(similar(f0),Float32(0))
        self.df_temp = fill!(similar(f0),Float32(0))
        self.temp = zeros(Float32,length(f0),length(f0))
        self.L = similar(self.temp)
        self.U = similar(self.temp)

        return self
    end

end