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
function (Euler::EulerStruct)(df::Vector{F},f::Vector{F},dt0,dt,t) where F<:AbstractFloat

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)

    # reset arrays
    fill!(df,zero(eltype(df)))
    fill!(Euler.df,zero(eltype(Euler.df)))
    fill!(Euler.temp,zero(eltype(Euler.temp)))
    if Euler.Implicit
        fill!(Euler.Jac,zero(eltype(Euler.Jac)))
    end
        
    # add binary terms to temp (jacobians are added in update_Big_Bin! if implicit)
    if isempty(Euler.PhaseSpace.Binary_list) == false
        update_Big_Bin!(Euler,f)
        @. Euler.temp += Euler.M_Bin_Mul_Step 
    end
    # add emission terms to temp (jacobians are added in update_Big_Emi! if implicit)
    if isempty(Euler.PhaseSpace.Emi_list) == false
        update_Big_Emi!(Euler,f)
        @. Euler.temp += Euler.M_Emi_Step
    end
    # add flux terms to temp
    #  @. Euler.temp -= Euler.FluxM.B_Flux
    #  @. Euler.temp -= Euler.FluxM.C_Flux
    #  @. Euler.temp -= Euler.FluxM.D_Flux
    #@. Euler.temp -= Euler.FluxM.K_Flux 
    #@. Euler.temp -= Euler.FluxM.J_Flux 
    #@. Euler.temp -= Euler.FluxM.I_Flux
    #@. Euler.temp -= Euler.FluxM.I_Flux + Euler.FluxM.J_Flux + Euler.FluxM.K_Flux
    Euler.temp -= Euler.FluxM.F_Flux
    if Euler.Implicit
        #@. Euler.Jac -= Euler.FluxM.K_Flux
        #@. Euler.Jac -= Euler.FluxM.J_Flux
        #@. Euler.Jac -= Euler.FluxM.I_Flux
        @. Euler.Jac -= Euler.FluxM.F_Flux
    end
    # phase space correction for non-uniform time stepping only applied to spatial coordinate fluxes and interactions 
    if Euler.PhaseSpace.Time.t_grid != "u" 
        Euler.temp .*= dt / dt0
        if Euler.Implicit
            Euler.Jac .*= dt / dt0
        end
    end
    # add time fluxes to temp
    #@. Euler.temp -= Euler.FluxM.Ap_Flux
    #@. Euler.temp -= Euler.FluxM.Am_Flux
    @. Euler.temp -= Euler.FluxM.Ap_Flux + Euler.FluxM.Am_Flux

    if isinf(sum(Euler.temp))
        error("overflow in arrays")
        #@. g.temp = g.temp*(g.temp!=Inf)
    end

    if Euler.Implicit
        # TOFIX add fluxes inside mul!
        mul!(Euler.df_temp,Euler.temp,f)
        println("t = $t")
        println("cond = $(cond(Euler.FluxM.Ap_Flux .+ Euler.Jac))")
        lu!(Euler.LU,(Euler.FluxM.Ap_Flux .+ Euler.Jac))
        @. Euler.LU.ipiv = Euler.LU.ipiv
        ldiv!(Euler.df,Euler.LU,Euler.df_temp)
        #@. g.Jac = 2*g.M_Bin_Mul_Step
        #g.df .= (I-dt*g.Jac)\g.df
        #@. df = dt*g.df

        #if t > 1e-16
            #println("$(Euler.LU)")
        #    println("$(Euler.df)")
        #    error("here")
        #end

    else
        #= mul!(Euler.df_temp,Euler.temp,f)
        if isdiag(Euler.FluxM.Ap_Flux)
            ldiv!(Euler.df,factorize(Euler.FluxM.Ap_Flux),Euler.df_temp)
        else
            lu!(Euler.LU,Euler.FluxM.Ap_Flux)
            ldiv!(Euler.df,Euler.LU,Euler.df_temp)
        end =#

        
        #if isdiag(Euler.FluxM.Ap_Flux)
            Euler.df_temp .= diag(Euler.FluxM.Ap_Flux) # ASSUMING Ap_Flux is DIAGONAL, TO BE UPDATED LATER !!!!!!!!!!!
            Euler.temp ./= Euler.df_temp
            #ldiv!(Euler.temp,factorize(Euler.FluxM.Ap_Flux),Euler.temp)
            #println("$(sum(Euler.temp))")
            # CFL ish check
            #println("t = $t")
            #println("max λ = $(maximum(abs.(eigvals(Euler.temp))))")
            #println("min λ = $(minimum(abs.(eigvals(Euler.temp))))")
        #else
        #    lu!(Euler.LU,Euler.FluxM.Ap_Flux)
        #    ldiv!(Euler.temp,Euler.LU,Euler.temp)
        #end
        mul!(Euler.df,Euler.temp,f)

    end
    
    @. df = Euler.df

    #println("$df")
    #error("")

end

function update_Big_Bin!(method::SteppingMethodType,f)

    if size(method.BigM.M_Bin) != (length(f)^2,length(f))
        error("M_Bin is not the correct size")
    end

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    FluxM = method.FluxM 

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = FluxM.Vol

    # Thanks to Emma Godden for fixing a bug here
    temp = reshape(method.M_Bin_Mul_Step,length(f)*length(f))

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

        tempView = @view temp[n_momentum^2*off_space+1:n_momentum^2*(off_space+1)]
        fView = @view f[n_momentum*off_space+1:n_momentum*(off_space+1)]

        mul!(tempView,method.BigM.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        # multiply by volume element
        tempView .*= Vol[off_space+1]

        # assign jacobian elements
        #=if method.Implicit
            JacView = @view method.Jac[n_momentum*off_space+1:n_momentum*(off_space+1),n_momentum*off_space+1:n_momentum*(off_space+1)]
            @. JacView += 2*method.M_Bin_Mul_Step 
        end 
        =#

    end

    if method.Implicit
        # assign jacobian elements
        @. method.Jac += 2*method.M_Bin_Mul_Step 
    end

    return nothing

end

function update_Big_Emi!(method::SteppingMethodType,f)

    if size(method.BigM.M_Emi) != (length(f),length(f))
        error("M_Bin is not the correct size")
    end

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    FluxM = method.FluxM 

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = FluxM.Vol

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

        tempView = @view method.M_Emi_Step[n_momentum*off_space+1:n_momentum*(off_space+1),n_momentum*off_space+1:n_momentum*(off_space+1)]

        @. tempView = method.BigM.M_Emi * Vol[off_space+1]

        # multiply by volume element
        #tempView .*= Vol[off_space+1]

    end

    # assign jacobian elements
    if method.Implicit
        @. method.Jac += method.M_Emi_Step 
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