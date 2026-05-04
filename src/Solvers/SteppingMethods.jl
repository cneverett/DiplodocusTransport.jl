"""
    ForwardEuler(t_start,t_stop,dt,Verbose)

Forward (Explicit) Euler time-stepping method for the Boltzmann equation.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

    # scaling of time stepping

    dt_scale = dt / dt0

    # update momentum and space space using f at time t

    mul!(method.df_Flux,method.F_Flux,method.f)
    @. method.df = -method.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation, also resets df_Flux

    # create df_Emi due to emission terms
    if method.Emission_Interactions
        mul!(method.df_Emi,method.M_Emi,method.f)
        @. method.df += method.df_Emi
    end
        
    # create df_Bin due to binary interactions
    if method.Binary_Interactions
        update_Big_Bin!(method)
        @. method.df += method.df_Bin
    end

    #println("df: $(method.df)")

    @. method.df *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
    # add injection term
    @. method.df += method.df_Inj * dt_scale 

    if !isnothing(method.df_mask)
        @. method.df *= method.df_mask
    end 

    if !isfinite(sum(method.df))
        println("non-finite value in df calculation")
    end

    # Cr (CFL) condition check
    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 3

        Cr_Bin = 0.0
        Cr_Emi = 0.0
        Cr_Flux = 0.0

        sum_f = sum(method.f)
        if sum_f != 0.0

            # Binary CFL
            if method.Binary_Interactions
                @. method.df_tmp = method.df_Bin / method.f 
                @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
                Cr_Bin = -minimum(method.df_tmp) # non-allocating and GPU compatible without CPU fallback
            end

            # Emission CFL
            if method.Emission_Interactions
                @. method.df_tmp = method.df_Emi / method.f 
                @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
                Cr_Emi = -minimum(method.df_tmp) # non-allocating and GPU compatible without CPU fallback
            end

            # Flux CFL
            @. method.df_tmp = -method.df_Flux / method.f 
            @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
            Cr_Flux = -minimum(method.df_tmp) # non-allocating and GPU compatible without CPU fallback

        end   

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt,sigdigits=3)), dt_adapted = $(round(dt * 0.8 / Cr,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt,sigdigits=3)), dt_adapted = $(round(dt * 0.8 / Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt=$dt")
    end
    if Verbose > 0
        flush(stdout)
    end

    # Basic timestep control based on CFL condition
    #=if Cr > 1.0 
        dt = dt * 0.5 
        dt_scale = dt / dt0
        adaptive_factor = 0.5
    elseif Cr < 1/2^9 && Cr > 0.0
        dt = dt * 2.0
        dt_scale = dt / dt0
        adaptive_factor = 2.0        
    else
        adaptive_factor = 1.0
    end=#
    # Slightly better timestep control based on CFL condition
    if Cr != 0.8 
        dt = dt * 0.8 / Cr
        dt_scale = dt / dt0
        adaptive_factor = 0.8 / Cr
    end

    # will we reached the next t_save?

    t_next = t_start + dt
    if abs(t_next - t_stop) <= eps(max(abs(t_next), abs(t_stop)))
        save = true
    elseif t_next >= t_stop
        adaptive_factor *= (t_stop - t_start) / dt # adjust adaptive factor for final time step to ensure we end exactly at t_stop 
        dt = t_stop - t_start
        save = true
    else
        save = false
    end

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # update state vector f with momentum, space and injection updates
    @. method.f += method.df * adaptive_factor 
    #@. method.f += (method.df + method.df_Inj * dt_scale) * adaptive_factor # add injection term and apply adaptive factor to total update 

    # remove masked off domain regions
    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    # removing negative values (values less than 1f-28 for better stability)
    @. method.f = method.f * (method.f>=method.n_cut) * sign(method.f)
    # hacky fix for inf values
    @. method.f = method.f * (method.f!=Inf)

    return dt,save

end

"""
    ForwardSymplecticEuler(t_start,t_stop,dt,Verbose)

Forward Symplectic (Semi-Implicit) Euler time-stepping method for the Boltzmann equation. Symplectic integrator updates momentum space first then physical space.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+\\left(M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right)\\mathcal{V}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardSymplecticEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    dt0 = method.dt0

    # will we reached the next t_save? TODO: adjust this for adaptive time stepping, currently assumes constant time stepping

        if isapprox(t_start + dt, t_stop)
            save = true
        elseif t_start + dt >= t_stop
            dt = t_stop - t_start
            save = true
        else
            save = false
        end

    # scaling of time stepping

        dt_scale = dt / dt0

    # update momentum space using f at time t

        # create df_PFlux due to momentum flux terms
        mul!(method.df_PFlux,method.P_Flux,method.f)
        @. method.df_Momentum = -method.df_PFlux # minus sign as flux terms are on RHS of Boltzmann equation, also resets df_Momentum

        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f)
            @. method.df_Momentum += method.df_Emi
        end
            
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method)
            @. method.df_Momentum += method.df_Bin
        end

        @. method.df_Momentum *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector

        # update f_momentum (f after momentum update)
        @. method.f_Momentum = method.f + method.df_Momentum
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Momentum = method.f_Momentum*(method.f_Momentum>=1f-10)

    # update physical space using f_Momentum (f after momentum update)

        # create df_XFlux due to space flux terms
        mul!(method.df_XFlux,method.X_Flux,method.f_Momentum)
        @. method.df_Space = -method.df_XFlux * method.invAp_Flux * dt_scale # minus sign as flux terms are on RHS of Boltzmann equation, also resets df_Space

        # update f_Space (f after momentum and space updates)
        @. method.f_Space = method.f_Momentum + method.df_Space
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Space = method.f_Space*(method.f_Space>=1f-20)

    # CFL condition check TODO: add adaptive time stepping based on CFL condition 

        if Verbose == 1 || Verbose == 2 || Verbose == 3

            Cr = 0.0
            Cr_Bin = 0.0
            Cr_Emi = 0.0
            Cr_PFlux = 0.0
            Cr_Momentum = 0.0
            Cr_Space = 0.0

            # Cr (CFL) condition check
            if sum(method.f) != 0.0

                if Verbose == 3
            
                    # Binary CFL
                    if method.Binary_Interactions
                        @. method.df_tmp = method.df_Bin * method.invAp_Flux * dt_scale / method.f
                        Cr_Bin = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    end

                    # Emission CFL
                    if method.Emission_Interactions
                        @. method.df_tmp = method.df_Emi * method.invAp_Flux * dt_scale / method.f
                        Cr_Emi = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    end

                    # P Flux CFL
                    @. method.df_tmp = -method.df_PFlux * method.invAp_Flux * dt_scale / method.f
                    Cr_PFlux = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback

                    # Momentum CFL
                    @. method.df_tmp = method.df_Momentum /  method.f
                    Cr_Momentum = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    # Space CFL
                    @. method.df_tmp = method.df_Space / method.f_Momentum
                    Cr_Space = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                end

                # Cr is calculated for the entire time step (momentum and space updates)
                # f_tmp - f = df of the momentum step
                @. method.df_tmp = (method.df_Space + method.df_Momentum) / method.f 
                Cr = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) 

            end   

            if Verbose == 1 && Cr > 1.0
                println("Cr = $(round(Cr, sigdigits=3)), t=$t_start, dt=$dt, system may be unstable")
            elseif Verbose == 2
                println("Cr = $(round(Cr, sigdigits=3)), t=$t_start, dt=$dt")
            elseif Verbose == 3
                println("Cr = $(round(Cr, sigdigits=3)), Cr_Bin = $(round(Cr_Bin, sigdigits=3)), Cr_Emi = $(round(Cr_Emi, sigdigits=3)), Cr_PFlux = $(round(Cr_PFlux, sigdigits=3)), Cr_Momentum = $(round(Cr_Momentum, sigdigits=3)), Cr_Space = $(round(Cr_Space, sigdigits=3)),  t=$t_start, t_save=$t_stop , dt=$dt")
            end
        end

    # update state vector f with momentum, space and injection updates
    
        @. method.f = method.f_Space + method.df_Inj * dt_scale
        # removing negative values (values less than p_cut for better stability) and ensure positivity for CFL calculations
        @. method.f = method.f*(method.f>=method.p_cut)*sign(method.f)
        # hacky fix for inf values
        @. method.f = method.f*(method.f!=Inf)

    # return dt and save

    return dt,save

end

"""
    BackwardsEuler(dg,g,t,dt)

Backwards (Implicit) Euler time-stepping method for the Boltzmann equation. 

"""
function (BackwardEuler::BackwardEulerStruct)(dt0,dt,t;Verbose::Bool=false)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)

    # reset arrays
    fill!(BackwardEuler.df,zero(eltype(BackwardEuler.df)))
    #fill!(Euler.temp,zero(eltype(Euler.temp)))
    if BackwardEuler.Implicit
        fill!(BackwardEuler.Jac,zero(eltype(BackwardEuler.Jac)))
    end
        
    # create df_Bin due to binary interactions (jacobians are added in update_Big_Bin! if implicit)
    if BackwardEuler.Binary_Interactions
        update_Big_Bin!(BackwardEuler)
        @. BackwardEuler.df += BackwardEuler.df_Bin
        if !isfinite(sum(BackwardEuler.df_Bin))
            println("overflow in df_Bin calculation, $(sum(BackwardEuler.df_Bin))")
        end
    end
    # create df_Emi due to emission terms  (jacobians are added in update_Big_Emi! if implicit)
    if BackwardsEuler.Emission_Interactions
        update_Big_Emi!(BackwardEuler)
        @. BackwardEuler.df += BackwardEuler.df_Emi
    end

    # create df_Flux due to space and momentum flux terms
    mul!(BackwardEuler.df_Flux,BackwardEuler.F_Flux,BackwardEuler.f)
    @. BackwardEuler.df -= BackwardEuler.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation
    if !isfinite(sum(BackwardEuler.df_Flux))
        println("overflow in df_Flux calculation, $(sum(BackwardEuler.df_Flux))")
    end
    if BackwardEuler.Implicit
        @. BackwardEuler.Jac -= BackwardEuler.F_Flux
    end

    # Add injection term 
    @. BackwardEuler.df += BackwardEuler.df_Inj
    # phase space correction for non-uniform time stepping only applied to spatial coordinate fluxes and interactions 
    if BackwardEuler.PhaseSpace.Time.t_grid != "u" 
        BackwardEuler.df .*= dt / dt0
        if BackwardEuler.Implicit
            BackwardEuler.Jac .*= dt / dt0
        end
    end

    # df_Flux due to time fluxes TODO: can remove this step if system is stationary therefore Ap=-Am. This will also allow more timestep control
    #mul!(Euler.df_Flux,Euler.FluxM.Am_Flux+Euler.FluxM.Ap_Flux,f)
    #@. Euler.df -= Euler.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation

    if !isfinite(sum(BackwardEuler.df))
        println("overflow in df calculation")
    end


        # TODO: Update Implicit
        #mul!(Euler.df_temp,Euler.temp,f)
        #println("t = $t")
        #println("cond = $(cond(Euler.FluxM.Ap_Flux .+ Euler.Jac))")
        #lu!(Euler.LU,(Euler.FluxM.Ap_Flux .+ Euler.Jac))
        #@. Euler.LU.ipiv = Euler.LU.ipiv
        #ldiv!(Euler.df,Euler.LU,Euler.df_temp)


    # Add injection term 
    if BackwardEuler.PhaseSpace.Time.t_grid != "u" 
        @. BackwardEuler.df += BackwardEuler.df_Inj * dt / dt0
    else
        @. BackwardEuler.df += BackwardEuler.df_Inj
    end
    
    # update state vector f
    @. BackwardEuler.f += BackwardEuler.df
    # removing negative values (values less than 1f-28 for better stability)
    @. BackwardEuler.f = BackwardEuler.f*(BackwardEuler.f>=1f-28)
    # hacky fix for inf values
    @. BackwardEuler.f = BackwardEuler.f*(BackwardEuler.f!=Inf)


end

function update_Big_Bin!(method::AbstractSteppingMethod)
    
    n_momentum=size(method.M_Bin,2)

    for off_space in method.Bin_Domain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view(method.f[start_idx:end_idx])

        df_BinView = @view method.df_Bin[start_idx:end_idx]
        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        @inbounds vol = method.Vol[off_space+1]
        mul!(df_BinView,method.M_Bin_Mul_Step,fView,vol,zero(eltype(method.df_Bin)))

    end

    #=
    # TODO: Update Implicit
    if method.Implicit
        # assign jacobian elements
        @. method.Jac += 2*method.M_Bin_Mul_Step*view(Vol,off_space+1)
    end=#

    return nothing

end

function update_Big_Emi!(method::AbstractSteppingMethod)

    @assert size(method.M_Emi) == (length(f),length(f)) "M_Emi is not the correct size"

    PhaseSpace = method.PhaseSpace
    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
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

        fView = @view(method.f[start_idx:end_idx])
        df_EmiView = @view method.df_Emi[start_idx:end_idx]

        M_EmiView = @view method.M_Emi[start_idx:end_idx,start_idx:end_idx] # TODO: make M_Emi block diagonal to save memory 

        mul!(df_EmiView,M_EmiView,fView)

        # multiply by volume element
        df_EmiView .*= Vol[off_space+1]

    end

    #=
    TODO: Update Implicit
    # assign jacobian elements
    if method.Implicit
        @. method.Jac += method.M_Emi
    end=#

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