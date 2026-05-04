abstract type AbstractSteppingMethod <: Function end

# Struct for storing the Boltzmann equation and its solution
mutable struct ForwardEulerStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},MET<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: AbstractSteppingMethod

    PhaseSpace::PhaseSpaceStruct

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    M_Bin::MBT
    Bin_Domain::BD

    M_Emi::MET

    F_Flux::SMT
    invAp_Flux::VT                  # inv Ap flux for time stepping
    Vol::Vector{T}

    Adaptive::Bool
    Implicit::Bool
    dt0::T
    n_cut::T

    step::Int64

    M_Bin_Mul_Step::MT             # temporary array for matrix multiplication of binary terms
    M_Bin_Mul_Step_reshape::VT     # temporary array for reshaped matrix multiplication of binary terms
    f_init::VT                     # initial distribution function (used by solver to define output struct)
    f::VT                          # current distribution function
    df::VT                         # change in distribution function
    df_Bin::VT                     # change in distribution function due to binary interactions
    df_Emi::VT                     # change in distribution function due to emission interactions
    df_Flux::VT                    # change in distribution function due to fluxes
    df_Inj::VT                     # change in distribution function due to injection of particles
    df_tmp::VT                     # temporary array the size of f for CFL calculations
    f_mask::FD      # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
    df_mask::DFD     # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

    function ForwardEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

        Backend = getfield(Main,Symbol("Backend"))
        Precision = getfield(Main,Symbol("Precision"))

        Momentum = PhaseSpace.Momentum
        Spacetime = PhaseSpace.Spacetime
        x_num = Spacetime.x_num
        y_num = Spacetime.y_num
        z_num = Spacetime.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list

        n_space = x_num+y_num+z_num
        n_momentum = 0
        for i in eachindex(px_num_list)
            n_momentum += px_num_list[i]*py_num_list[i]*pz_num_list[i]
        end

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        Binary_Interactions = !isempty(BinM.Binary_list) && !isnothing(BinM.Domain)
        Emission_Interactions = !isempty(EmiM.Emission_list)

        Bin_Domain = BinM.Domain

        if Binary_Interactions
            M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
        else
            M_Bin_Mul_Step = zeros(Backend,Precision,0,0)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,0)
        end
        df = zeros(Backend,Precision,length(Initial))
        df_Bin = zeros(Backend,Precision,length(Initial))
        df_Emi = zeros(Backend,Precision,length(Initial))
        df_Flux = zeros(Backend,Precision,length(Initial))
        df_tmp = zeros(Backend,Precision,length(Initial))

        Vol = FluxM.Vol

        if Backend isa CPUBackend
            f_init = convert(Vector{Precision},Initial)
            M_Bin = BinM.M_Bin
            M_Emi = EmiM.M_Emi
            F_Flux = FluxM.X_Flux + FluxM.P_Flux # sum of space and momentum fluxes
            invAp_Flux = 1 ./ FluxM.Ap_Flux # invert Ap Flux
            f = convert(Vector{Precision},copy(Initial))
            df_Inj = convert(Vector{Precision},copy(Injection))
        elseif Backend isa CUDABackend
            f_init = CuArray(convert(Vector{Precision},copy(Initial)))
            if BinM.M_Bin isa AbstractSparseArray
                M_Bin = CuSparseMatrixCSC(BinM.M_Bin)
            else
                M_Bin = CuArray(BinM.M_Bin)
            end
            if EmiM.M_Emi isa AbstractSparseArray
                M_Emi = CuSparseMatrixCSC(EmiM.M_Emi)
            else
                M_Emi = CuArray(EmiM.M_Emi)
            end
            F_Flux = CuSparseMatrixCSC(FluxM.X_Flux + FluxM.P_Flux) # sum of space and momentum fluxes
            invAp_Flux = CuArray(1 ./ FluxM.Ap_Flux)
            f = CuArray(convert(Vector{Precision},copy(Initial)))
            df_Inj = CuArray(convert(Vector{Precision},copy(Injection)))
        else
            error("Backend type not recognized.")
        end

        if !isnothing(DistributionDomainMask)
            f_mask = ones(Precision,length(Initial))
            for off_space_idx in DistributionDomainMask
                for species_idx in eachindex(PhaseSpace.name_list)
                LocationSpeciesToStateVector(f_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= Precision(0.0)
                end
            end
            if Backend isa CUDABackend
                f_mask = CuArray(f_mask)
            end
        else
            f_mask = nothing
        end

        if !isnothing(DeltaDistributionDomainMask)
            df_mask = ones(Precision,length(Initial))
            for off_space_idx in DeltaDistributionDomainMask
                for species_idx in eachindex(PhaseSpace.name_list)
                LocationSpeciesToStateVector(df_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= Precision(0.0)
                end
            end
            if Backend isa CUDABackend
                df_mask = CuArray(df_mask)
            end
        else
            df_mask = nothing
        end

        ###### Actually Build the Struct with Concrete Types ######

        self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(M_Emi),typeof(F_Flux),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

        self.Implicit = false
        self.Adaptive = Adaptive
        self.dt0 = PhaseSpace.Spacetime.dt0
        self.n_cut = n_cut
        self.step = 0
        self.Binary_Interactions = Binary_Interactions
        self.Emission_Interactions = Emission_Interactions
        self.PhaseSpace = PhaseSpace
        self.Bin_Domain = BinM.Domain
        self.f_init = f_init
        self.M_Bin = M_Bin
        self.M_Emi = M_Emi
        self.F_Flux = F_Flux
        self.invAp_Flux = invAp_Flux
        self.Vol = Vol
        self.f = f
        self.df_Inj = df_Inj
        self.M_Bin_Mul_Step = M_Bin_Mul_Step
        self.M_Bin_Mul_Step_reshape = M_Bin_Mul_Step_reshape # Thanks to Emma Godden for fixing a bug here
        self.df = df
        self.df_Bin = df_Bin
        self.df_Emi = df_Emi
        self.df_Flux = df_Flux
        self.df_tmp = df_tmp
        self.f_mask = f_mask
        self.df_mask = df_mask

        return self
    end

end

mutable struct ForwardSymplecticEulerStruct{T<:AbstractFloat} <: AbstractSteppingMethod

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
    dt0::T
    p_cut::T

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

    function ForwardSymplecticEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,p_cut::Float64=1e-45)

        Backend = getfield(Main,Symbol("Backend"))
        Precision = getfield(Main,Symbol("Precision"))

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        self = new{Precision}()

        self.Adaptive = Adaptive
        self.Implicit = false
        self.dt0 = convert(Precision,PhaseSpace.Spacetime.dt0)
        self.p_cut = convert(Precision,p_cut)

        self.Binary_Interactions = !isempty(BinM.Binary_list)
        self.Emission_Interactions = !isempty(EmiM.Emission_list)

        self.PhaseSpace = PhaseSpace

        self.Bin_Domain = BinM.Domain

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
        Spacetime = PhaseSpace.Spacetime
        x_num = Spacetime.x_num
        y_num = Spacetime.y_num
        z_num = Spacetime.z_num
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
mutable struct BackwardEulerStruct{T<:AbstractFloat} <: AbstractSteppingMethod

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
