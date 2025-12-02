"""
    Euler(dg,g,t,dt)

Explicit and Implicit Euler time-stepping method for the Boltzmann equation. Explicit/Implicit is defined by the boolean `Implicit` flag when defining the EulerStruct.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````


# Implicit method:


"""
function (Euler::EulerStruct)(dt0,dt,t)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)

    # reset arrays
    fill!(Euler.df,zero(eltype(Euler.df)))
    fill!(Euler.temp,zero(eltype(Euler.temp)))
    if Euler.Implicit
        fill!(Euler.Jac,zero(eltype(Euler.Jac)))
    end
        
    # create df_Bin due to binary interactions (jacobians are added in update_Big_Bin! if implicit)
    if isempty(Euler.PhaseSpace.Binary_list) == false
        update_Big_Bin!(Euler)
        @. Euler.df += Euler.df_Bin
    end
    # create df_Emi due to emission terms  (jacobians are added in update_Big_Emi! if implicit)
    if isempty(Euler.PhaseSpace.Emi_list) == false
        update_Big_Emi!(Euler)
        @. Euler.df += Euler.df_Emi
    end

    # create df_Flux due to space and momentum flux terms
    mul!(Euler.df_Flux,Euler.F_Flux,Euler.f)
    @. Euler.df -= Euler.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation
    if Euler.Implicit
        @. Euler.Jac -= Euler.F_Flux
    end

    # phase space correction for non-uniform time stepping only applied to spatial coordinate fluxes and interactions 
    if Euler.PhaseSpace.Time.t_grid != "u" 
        Euler.df .*= dt / dt0
        if Euler.Implicit
            Euler.Jac .*= dt / dt0
        end
    end

    # df_Flux due to time fluxes TODO: can remove this step if system is stationary therefore Ap=-Am. This will also allow more timestep control
    #mul!(Euler.df_Flux,Euler.FluxM.Am_Flux+Euler.FluxM.Ap_Flux,f)
    #@. Euler.df -= Euler.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation

    if isinf(sum(Euler.df))
        println("overflow in df calculation")
        #@. g.temp = g.temp*(g.temp!=Inf)
    end

    if Euler.Implicit

        # TODO: Update Implicit
        #mul!(Euler.df_temp,Euler.temp,f)
        #println("t = $t")
        #println("cond = $(cond(Euler.FluxM.Ap_Flux .+ Euler.Jac))")
        #lu!(Euler.LU,(Euler.FluxM.Ap_Flux .+ Euler.Jac))
        #@. Euler.LU.ipiv = Euler.LU.ipiv
        #ldiv!(Euler.df,Euler.LU,Euler.df_temp)


    else

        @. Euler.df /= Euler.Ap_Flux # Assumes Ap_flux is diagonal and stored as a vector

        # Cr (CFL) condition check
        @. Euler.df_tmp = Euler.df / Euler.f 
        replace!(Euler.df_tmp,Inf32=>0f0,NaN32=>0f0,-Inf32=>0f0,-NaN32=>0f0)
        Cr = maximum(abs.(Euler.df_tmp))

        if Euler.Verbose
            println("Cr = $Cr")
        elseif Cr > 1.0
            println("Cr = $Cr, system may be unstable")
        end

    end
    
    # update state vector f
    @. Euler.f += Euler.df

    # removing negative values (values less than 1f-28 for better stability)
    @. Euler.f = Euler.f*(Euler.f>=1f-28)
    # hacky fix for inf values
    @. Euler.f = Euler.f*(Euler.f!=Inf)

    #println("$df")
    #error("")

end

function update_Big_Bin!(method::SteppingMethodType)

    f = method.f
    @assert size(method.M_Bin) == (length(f)^2,length(f)) "M_Bin is not the correct size"

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = method.Vol

    # Thanks to Emma Godden for fixing a bug here
    #temp = reshape(method.M_Bin_Mul_Step,length(f)*length(f))

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1 # starts at 0

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view f[start_idx:end_idx]
        df_BinView = @view method.df_Bin[start_idx:end_idx]
        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape
        mul!(df_BinView,method.M_Bin_Mul_Step,fView)

        # multiply by volume element
        df_BinView .*= Vol[off_space+1]

    end

    if method.Implicit
        # assign jacobian elements
        @. method.Jac += 2*method.M_Bin_Mul_Step*Vol[off_space+1]
    end

    return nothing

end

function update_Big_Emi!(method::SteppingMethodType)

    f = method.f

    @assert size(method.M_Emi) == (length(f),length(f)) "M_Emi is not the correct size"

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = method.Vol

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view f[start_idx:end_idx]
        df_EmiView = @view method.df_Emi[start_idx:end_idx]

        M_EmiView = @view method.M_Emi[start_idx:end_idx,start_idx:end_idx] # TODO: make M_Emi block diagonal to save memory 

        mul!(df_EmiView,M_EmiView,fView)

        # multiply by volume element
        df_EmiView .*= Vol[off_space+1]

    end

    # assign jacobian elements
    if method.Implicit
        @. method.Jac += method.M_Emi
    end

    return nothing

end

# to be removed in julia 1.12#

function generic_lufact!(A::AbstractMatrix{T}, pivot::Union{RowMaximum,NoPivot,RowNonZero} = lupivottype(T), ipiv::AbstractVector{LinearAlgebra.BlasInt} = Vector{LinearAlgebra.BlasInt}(undef,min(size(A)...));
                         check::Bool = true, allowsingular::Bool = false) where {T}
    check && LAPACK.chkfinite(A)
    # Extract values
    m, n = size(A)
    minmn = min(m,n)

    # Initialize variables
    info = 0
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if pivot === RowMaximum() && k < m
                amax = abs(A[k, k])
                for i = k+1:m
                    absi = abs(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            elseif pivot === RowNonZero()
                for i = k:m
                    if !iszero(A[i,k])
                        kp = i
                        break
                    end
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    if pivot === NoPivot()
        # Use a negative value to distinguish a failed factorization (zero in pivot
        # position during unpivoted LU) from a valid but rank-deficient factorization
        info = -info
    end
    check && LinearAlgebra._check_lu_success(info, allowsingular)
    return LU{T,typeof(A),typeof(ipiv)}(A, ipiv, convert(LinearAlgebra.BlasInt, info))
end

function lu!(F::LU{<:Any,<:AbstractMatrix}, A; check::Bool = true, allowsingular::Bool = false)
    copyto!(F.factors, Float64.(A))
    return generic_lufact!(F.factors, LinearAlgebra.lupivottype(eltype(A)), F.ipiv; check, allowsingular)
end