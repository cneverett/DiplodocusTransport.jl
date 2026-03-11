# Struct for storing the Boltzmann equation and its solution
mutable struct ForwardEulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    M_Bin::AbstractMatrix{T}
    Bin_Domain::Union{Vector{Int64},Nothing}

    M_Emi::AbstractMatrix{T}

    F_Flux::AbstractSparseArray{T,<:Integer,2}
    invAp_Flux::AbstractVector{T}                   # inv Ap flux for time stepping
    Vol::AbstractVector{T}

    Adaptive::Bool
    Implicit::Bool
    dt_initial::T

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

    function ForwardEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false)

        Backend = getfield(Main,Symbol("Backend"))
        Precision = getfield(Main,Symbol("Precision"))

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        self = new{Precision}()

        self.Implicit = false
        self.Adaptive = Adaptive
        self.dt_initial = convert(Precision,PhaseSpace.Time.dt_initial)

        self.Binary_Interactions = !isempty(BinM.Binary_list)
        self.Emission_Interactions = !isempty(EmiM.Emission_list)

        self.PhaseSpace = PhaseSpace

        self.Bin_Domain = BinM.Domain

        self.f_init = copy(Initial)

        invAp_Flux = 1 ./ FluxM.Ap_Flux # invert Ap Flux

        if Backend isa CPUBackend
            self.M_Bin = BinM.M_Bin
            self.M_Emi = EmiM.M_Emi
            self.F_Flux = FluxM.X_Flux + FluxM.P_Flux # sum of space and momentum fluxes
            self.invAp_Flux = invAp_Flux
            self.Vol = FluxM.Vol
            self.f = convert(Vector{Precision},Initial)
            self.df_Inj = convert(Vector{Precision},Injection)
        elseif Backend isa CUDABackend
            self.M_Bin = CuArray(BinM.M_Bin)
            self.M_Emi = CuArray(EmiM.M_Emi)
            self.F_Flux = CuSparseMatrixCSC(FluxM.X_Flux + FluxM.P_Flux) # sum of space and momentum fluxes
            self.invAp_Flux = CuArray(invAp_Flux)
            self.Vol = CuArray(FluxM.Vol)
            self.f = CuArray(Initial)
            self.df_Inj = CuArray(Injection)
        else
            error("Backend type not recognized.")
        end  

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

        if !isempty(BinM.Binary_list)
            self.M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            self.M_Bin_Mul_Step_reshape = reshape(self.M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        end
        self.df = zeros(Backend,Precision,length(Initial))
        self.df_Bin = zeros(Backend,Precision,length(Initial))
        self.df_Emi = zeros(Backend,Precision,length(Initial))
        self.df_Flux = zeros(Backend,Precision,length(Initial))
        self.df_tmp = zeros(Backend,Precision,length(Initial))

        return self
    end

end

mutable struct ForwardSymplecticEulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    M_Bin::AbstractMatrix{T}
    Bin_Domain::Union{Vector{Int64},Nothing}

    M_Emi::AbstractMatrix{T}

    X_Flux::AbstractSparseArray{T,<:Integer,2}
    P_Flux::AbstractSparseArray{T,<:Integer,2}
    invAp_Flux::AbstractVector{T}                   # inv Ap flux for time stepping
    Vol::AbstractVector{T}

    Adaptive::Bool
    Implicit::Bool
    dt_initial::T

    M_Bin_Mul_Step::AbstractMatrix{T}               # temporary array for matrix multiplication of binary terms
    M_Bin_Mul_Step_reshape::AbstractVector{T}       # temporary array for reshaped matrix multiplication of binary terms
    f_init::AbstractVector{T}                       # initial distribution function (used by solver to define output struct)
    f::AbstractVector{T}                            # current distribution function
    df::AbstractVector{T}                           # change in distribution function
    df_Momentum::AbstractVector{T}                  # change in distribution function due to momentum fluxes 
    df_Space::AbstractVector{T}                     # change in distribution function due to spatial fluxes
    df_Bin::AbstractVector{T}                       # change in distribution function due to binary interactions
    df_Emi::AbstractVector{T}                       # change in distribution function due to emission interactions
    df_XFlux::AbstractVector{T}                     # change in distribution function due to spatial fluxes
    f_Space::AbstractVector{T}                      # temporary array for f after spatial fluxes
    df_PFlux::AbstractVector{T}                     # change in distribution function due to momentum fluxes
    f_Momentum::AbstractVector{T}                      # temporary array for f after momentum fluxes
    df_Inj::AbstractVector{T}                       # change in distribution function due to injection of particles
    df_tmp::AbstractVector{T}                       # temporary array the size of f for CFL calculations 
    f_tmp::AbstractVector{T}                        # temporary array the size of f for symplectic Euler calculations  

    function ForwardSymplecticEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false)

        Backend = getfield(Main,Symbol("Backend"))
        Precision = getfield(Main,Symbol("Precision"))

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        self = new{Precision}()

        self.Adaptive = Adaptive
        self.Implicit = false
        self.dt_initial = convert(Precision,PhaseSpace.Time.dt_initial)

        self.Binary_Interactions = !isempty(BinM.Binary_list)
        self.Emission_Interactions = !isempty(EmiM.Emission_list)

        self.PhaseSpace = PhaseSpace

        self.Bin_Domain = BinM.Domain

        self.f_init = copy(Initial)

        invAp_Flux = 1 ./ FluxM.Ap_Flux # invert Ap Flux

        if Backend isa CPUBackend
            self.f_init = convert(Vector{Precision},Initial)
            self.M_Bin = BinM.M_Bin
            self.M_Emi = EmiM.M_Emi
            self.X_Flux = FluxM.X_Flux
            self.P_Flux = FluxM.P_Flux
            self.invAp_Flux = invAp_Flux
            self.Vol = FluxM.Vol
            self.f = convert(Vector{Precision},Initial)
            self.df_Inj = convert(Vector{Precision},Injection)
        elseif Backend isa CUDABackend
            self.f_init = CuArray(Initial)
            if BinM.M_Bin isa AbstractSparseMatrix
                self.M_Bin = CuSparseMatrixCSC(BinM.M_Bin)
            else
                self.M_Bin = CuArray(BinM.M_Bin)
            end
            if EmiM.M_Emi isa AbstractSparseMatrix
                self.M_Emi = CuSparseMatrixCSC(EmiM.M_Emi)
            else
                self.M_Emi = CuArray(EmiM.M_Emi)
            end
            self.X_Flux = CuSparseMatrixCSC(FluxM.X_Flux)
            self.P_Flux = CuSparseMatrixCSC(FluxM.P_Flux)
            self.invAp_Flux = CuArray(invAp_Flux)
            self.Vol = CuArray(FluxM.Vol)
            self.f = CuArray(Initial)
            self.df_Inj = CuArray(Injection)
        else
            error("Backend type not recognized.")
        end  

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

        if !isempty(BinM.Binary_list)
            self.M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            self.M_Bin_Mul_Step_reshape = reshape(self.M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        end
        self.df = zeros(Backend,Precision,length(Initial))
        self.df_Momentum = zeros(Backend,Precision,length(Initial))
        self.df_Space = zeros(Backend,Precision,length(Initial))
        self.df_Bin = zeros(Backend,Precision,length(Initial))
        self.df_Emi = zeros(Backend,Precision,length(Initial))
        self.df_XFlux = zeros(Backend,Precision,length(Initial))
        self.df_PFlux = zeros(Backend,Precision,length(Initial))
        self.f_Space = zeros(Backend,Precision,length(Initial))
        self.f_Momentum = zeros(Backend,Precision,length(Initial))
        self.df_tmp = zeros(Backend,Precision,length(Initial))
        self.f_tmp = zeros(Backend,Precision,length(Initial))

        return self
    end

end



# Struct for storing the Boltzmann equation and its solution
mutable struct BackwardEulerStruct{T<:AbstractFloat} <: SteppingMethodType

    PhaseSpace::PhaseSpaceStruct

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    M_Bin::AbstractMatrix{T}
    Bin_Domain::Vector{Int64}
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

    Jac::AbstractMatrix{T}                          # Jacobian matrix
    LU::LinearAlgebra.LU{T, AbstractMatrix{T}, AbstractVector{<:Integer}} # LU factorization of the matrix for implicit solving     

    function BackwardEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Implicit::Bool=false)

        Backend = getfield(Main,Symbol("Backend"))
        Precision = getfield(Main,Symbol("Precision"))

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        self = new{Precision}()

        self.Implicit = true

        self.Binary_Interactions = !isempty(BinM.Binary_list)
        self.Emission_Interactions = !isempty(EmiM.Emission_list)

        self.PhaseSpace = PhaseSpace
        
        Initial = convert(Precision,Initial)
        Injection = convert(Precision,Injection)

        self.Bin_Domain = BinM.Domain

        self.f_init = copy(Initial)

        if Backend isa CPUBackend
            self.M_Bin = BinM.M_Bin
            self.M_Emi = EmiM.M_Emi
            self.F_Flux = FluxM.F_Flux
            self.Ap_Flux = FluxM.Ap_Flux
            self.Vol = FluxM.Vol
            self.f = convert(Precision,Initial)
            self.df_Inj = convert(Precision,Injection)
        elseif Backend isa CUDABackend
            self.M_Bin = CuArray(BinM.M_Bin)
            self.M_Emi = CuArray(EmiM.M_Emi)
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
