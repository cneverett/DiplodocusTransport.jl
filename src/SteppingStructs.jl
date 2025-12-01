# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct
    BigM::BigMatricesStruct{Matrix{T}}
    FluxM::FluxMatricesStruct{Matrix{T},Vector{T}}

    Implicit::Bool

    M_Bin_Mul_Step::Matrix{T}               # temporary array for matrix multiplication of binary terms
    M_Bin_Mul_Step_reshape::Vector{T}       # temporary array for reshaped matrix multiplication of binary terms
    M_Emi_Step::Matrix{T}                   # temporary array for spatial evaluated emission terms
    Jac::Matrix{T}                          # jacobian for if Implicit==True
    df::Vector{T}                           # change in distribution function
    df_Bin::Vector{T}                       # change in distribution function due to binary interactions
    df_Emi::Vector{T}                       # change in distribution function due to emission interactions
    df_Flux::Vector{T}                      # change in distribution function due to fluxes
    df_Inj::Vector{T}                       # change in distribution function due to injection of particles
    temp::Matrix{T}
    LU::LinearAlgebra.LU{T, Matrix{T}, Vector{Int64}}

    function EulerStruct(f0::Vector{T},PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatricesStruct{Matrix{T}},Flux_Matrices::FluxMatricesStruct{Matrix{T},Vector{T}},Implicit::Bool) where T<:AbstractFloat

        self = new{T}()

        self.PhaseSpace = PhaseSpace
        self.BigM = Big_Matrices
        self.FluxM = Flux_Matrices

        self.Implicit = Implicit  

        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list

        n_space = x_num+y_num+z_num
        n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

        if isempty(PhaseSpace.Binary_list) == false
            self.M_Bin_Mul_Step = zeros(T,n_momentum,n_momentum)
            self.M_Bin_Mul_Step_reshape = zeros(T,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        end
        if isempty(PhaseSpace.Emi_list) == false
            self.M_Emi_Step = zeros(T,length(f0),length(f0))
        end
        self.df = zeros(T,length(f0))
        self.df_Bin = zeros(T,length(f0))
        self.df_Emi = zeros(T,length(f0))
        self.df_Flux = zeros(T,length(f0))
        self.temp = zeros(T,length(f0),length(f0))
        if Implicit
            self.Jac = zeros(T,length(f0),length(f0))
            self.LU = lu(zeros(T,length(f0),length(f0))+I)
        end


        return self
    end

end

