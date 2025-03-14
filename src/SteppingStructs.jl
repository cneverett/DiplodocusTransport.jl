# Struct for storing the Boltzmann equation and its solution
mutable struct Euler <: SteppingMethod

    t::Float32                  # the last time step time to calculate Δt
    dt::Float32                 # time step

    PhaseSpaceStruct::PhaseSpaceStruct
    BigM::BigMatrices
    FluxM::FluxMatrices

    Implicit::Bool

    A_Binary_Mul::Array{Float32,2}  # temporary array for matrix multiplication of binary terms
    Jac::Array{Float32,2}           # jacobian for if Implicit==True
    df::fType                       # change in distribution function
    df_temp::fType                  # change in distribution function
    temp::AbstractArray{Float32,2}
    L::AbstractArray{Float32,2}
    U::AbstractArray{Float32,2}

    function Euler(f0,PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatrices,Flux_Matrices::FluxMatrices,Implicit::Bool)

        self = new()

        self.PhaseSpace = PhaseSpace
        self.BigM = Big_Matrices
        self.FluxM = Flux_Matrices

        self.Implicit = Implicit  

        self.A_Binary_Mul = zeros(Float32,length(f0),length(f0))
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