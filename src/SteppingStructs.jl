# Struct for storing the Boltzmann equation and its solution
mutable struct EulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct

    M_Bin::AbstractMatrix{T}
    M_Emi::AbstractMatrix{T}

    F_Flux::AbstractSparseArray{T,<:Integer,2}
    Ap_Flux::AbstractVector{T}
    Vol::AbstractVector{T}

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
    # if Implicit is true
    Jac::AbstractMatrix{T}                          # Jacobian matrix
    LU::LinearAlgebra.LU{T, AbstractMatrix{T}, AbstractVector{<:Integer}} # LU factorization of the matrix for implicit solving     

    function EulerStruct(Initial::Vector{Float64},Injection::Vector{Float64},PhaseSpace::PhaseSpaceStruct,DataDirectory::String;Implicit::Bool=false,Backend::BackendType=CPUBackend(),Precision::T=Float32) where T<:Union{Float32,Float64}

        self = new{Precision}()

        self.PhaseSpace = PhaseSpace

        # Build Big and Flux Matrices
        BigM::BigMatricesStruct = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=false,Precision=Precision);
        FluxM::FluxMatricesStruct = BuildFluxMatrices(PhaseSpace,Precision=Precision);
        
        Initial = convert(Precision,Initial)
        Injection = convert(Precision,Injection)

        self.f_init = copy(Initial)

        if Backend isa CPUBackend
            self.M_Bin = BigM.M_Bin
            self.M_Emi = BigM.M_Emi
            self.F_Flux = FluxM.F_Flux
            self.Ap_Flux = FluxM.Ap_Flux
            self.Vol = FluxM.Vol
            self.f = convert(Precision,Initial)
            self.df_Inj = convert(Precision,Injection)
        elseif Backend isa CUDABackend
            self.M_Bin = CuArray(BigM.M_Bin)
            self.M_Emi = CuArray(BigM.M_Emi)
            self.F_Flux = CuSparseMatrixCSC(FluxM.F_Flux)
            self.Ap_Flux = CuArray(FluxM.Ap_Flux)
            self.Vol = CuArray(FluxM.Vol)
            self.f = CuArray(Initial)
            self.df_Inj = CuArray(Injection)
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
            self.M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            self.M_Bin_Mul_Step_reshape = reshape(self.M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        end
        self.df = zeros(Backend,Precision,length(Initial))
        self.df_Bin = zeros(Backend,Precision,length(Initial))
        self.df_Emi = zeros(Backend,Precision,length(Initial))
        self.df_Flux = zeros(Backend,Precision,length(Initial))
        self.df_tmp = zeros(Backend,Precision,length(Initial))
        if Implicit
            self.Jac = zeros(Backend,Precision,length(Initial),length(Initial))
            self.LU = lu(zeros(Backend,Precision,length(Initial),length(Initial))+I)
        end

        return self
    end

end
