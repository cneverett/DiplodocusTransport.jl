# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct

    M_Bin::AbstractMatrix{T}
    M_Emi::AbstractMatrix{T}

    F_Flux::AbstractSparseArray{T,Int64,2}
    Ap_Flux::AbstractVector{T}
    Vol::AbstractVector{T}

    #BigM::BigMatricesStruct{T}
    #FluxM::FluxMatricesStruct{T}

    Implicit::Bool

    M_Bin_Mul_Step::AbstractMatrix{T}               # temporary array for matrix multiplication of binary terms
    M_Bin_Mul_Step_reshape::AbstractVector{T}       # temporary array for reshaped matrix multiplication of binary terms
    f_init::AbstractVector{T}                       # initial distribution function (used by solver to define output struct)
    f::AbstractVector{T}                            # current distribution function
    df::AbstractVector{T}                           # change in distribution function
    df_Bin::AbstractVector{T}                       # change in distribution function due to binary interactions
    df_Emi::AbstractVector{T}                       # change in distribution function due to emission interactions
    df_Flux::AbstractVector{T}                      # change in distribution function due to fluxes
    df_Inj::AbstractVector{T}                       # change in distribution function due to injection of particles
    df_tmp::AbstractVector{T}                       # temporary array the size of f for CFL calculations
    temp::AbstractMatrix{T}
    # if Implicit is true
    Jac::AbstractMatrix{T}                          # Jacobian matrix
    LU::LinearAlgebra.LU{T, AbstractMatrix{T}, AbstractVector{Int64}} # LU factorization of the matrix for implicit solvingue      

    function EulerStruct(f0::Vector{T},PhaseSpace::PhaseSpaceStruct,Big_Matrices::BigMatricesStruct,Flux_Matrices::FluxMatricesStruct;Implicit::Bool=false,Backend::BackendType=CPUBackend(),Verbose::Bool=false) where T<:Union{Float32,Float64}

        self = new{T}()

        self.PhaseSpace = PhaseSpace
        #self.BigM = Big_Matrices
        #self.FluxM = Flux_Matrices

        self.f_init = copy(f0)

        if Backend isa CPUBackend
            self.M_Bin = Big_Matrices.M_Bin
            self.M_Emi = Big_Matrices.M_Emi
            self.F_Flux = Flux_Matrices.F_Flux
            self.Ap_Flux = Flux_Matrices.Ap_Flux
            self.Vol = Flux_Matrices.Vol
            self.f = copy(f0)
        elseif Backend isa CUDABackend
            self.M_Bin = CuArray(Big_Matrices.M_Bin)
            self.M_Emi = CuArray(Big_Matrices.M_Emi)
            self.F_Flux = CuSparseMatrixCSC(Flux_Matrices.F_Flux)
            self.Ap_Flux = CuArray(Flux_Matrices.Ap_Flux)
            self.Vol = CuArray(Flux_Matrices.Vol)
            self.f = CuArray(f0)
        else
            error("Backend type not recognized.")
        end
        self.Implicit = Implicit  

        Momentum = PhaseSpace.Momentum
        Space = PhaseSpace.Space
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list

        n_space = x_num+y_num+z_num
        n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

        if isempty(PhaseSpace.Binary_list) == false
            self.M_Bin_Mul_Step = zeros(Backend,T,n_momentum,n_momentum)
            self.M_Bin_Mul_Step_reshape = reshape(self.M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        end
        self.df = zeros(Backend,T,length(f0))
        self.df_Bin = zeros(Backend,T,length(f0))
        self.df_Emi = zeros(Backend,T,length(f0))
        self.df_Flux = zeros(Backend,T,length(f0))
        self.df_Inj = zeros(Backend,T,length(f0))
        self.df_tmp = zeros(Backend,T,length(f0))
        self.temp = zeros(Backend,T,length(f0),length(f0))
        if Implicit
            self.Jac = zeros(Backend,T,length(f0),length(f0))
            self.LU = lu(zeros(Backend,T,length(f0),length(f0))+I)
        end

        return self
    end

end
