"""
    ForwardEuler(dt0,dt,t,Verbose)

Forward (Explicit) Euler time-stepping method for the Boltzmann equation.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````
"""
function (ForwardEuler::ForwardEulerStruct)(dt0,dt,t,Verbose::Int64)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)

    #println(ForwardEuler.f)
    #println("t = $t")

    # reset arrays
    fill!(ForwardEuler.df,zero(eltype(ForwardEuler.df)))
    #fill!(Euler.temp,zero(eltype(Euler.temp)))
        
    # create df_Bin due to binary interactions
    if ForwardEuler.Binary_Interactions
        update_Big_Bin!(ForwardEuler)
        @. ForwardEuler.df += ForwardEuler.df_Bin
        if !isfinite(sum(ForwardEuler.df_Bin))
            println("non-finite value in df_Bin calculation, $(sum(ForwardEuler.df_Bin))")
        end
    end

    # create df_Emi due to emission terms
    if ForwardEuler.Emission_Interactions
        mul!(ForwardEuler.df_Emi,ForwardEuler.M_Emi,ForwardEuler.f)
        @. ForwardEuler.df += ForwardEuler.df_Emi
        if !isfinite(sum(ForwardEuler.df_Emi))
            println("non-finite value in df_Emi calculation, $(sum(ForwardEuler.df_Emi))")
        end
    end

    # create df_Flux due to space and momentum flux terms
    mul!(ForwardEuler.df_Flux,ForwardEuler.F_Flux,ForwardEuler.f)
    @. ForwardEuler.df -= ForwardEuler.df_Flux # minus sign as flux terms are on RHS of Boltzmann equation
    if !isfinite(sum(ForwardEuler.df_Flux))
        println("non-finite value in df_Flux calculation, $(sum(ForwardEuler.df_Flux))")
    end

    @. ForwardEuler.df *= ForwardEuler.invAp_Flux # Assumes Ap_flux is diagonal and stored as a vector

    # Add injection term 
    @. ForwardEuler.df += ForwardEuler.df_Inj

    # phase space correction for non-uniform time stepping only applied to spatial coordinate fluxes and interactions 
    if ForwardEuler.PhaseSpace.Time.t_grid != "u" 
        ForwardEuler.df .*= dt / dt0
    end

    if !isfinite(sum(ForwardEuler.df))
        println("non-finite value in df calculation")
    end

    if Verbose == 1 || Verbose == 2 || Verbose == 3

        Cr = 0.0
        Cr_Bin = 0.0
        Cr_Emi = 0.0
        Cr_Flux = 0.0

        # Cr (CFL) condition check
        if sum(ForwardEuler.f) != 0.0

            if Verbose == 3
        
                # Binary CFL
                if ForwardEuler.Binary_Interactions
                    @. ForwardEuler.df_tmp = ForwardEuler.df_Bin / ForwardEuler.f 
                    Cr_Bin = -minimum(filter(isfinite,ForwardEuler.df_tmp))
                end

                # Emission CFL
                if ForwardEuler.Emission_Interactions
                    @. ForwardEuler.df_tmp = ForwardEuler.df_Emi / ForwardEuler.f 
                    Cr_Emi = -minimum(filter(isfinite,ForwardEuler.df_tmp))
                end

                # Flux CFL
                @. ForwardEuler.df_tmp = -ForwardEuler.df_Flux / ForwardEuler.f 
                Cr_Flux = -minimum(filter(isfinite,ForwardEuler.df_tmp))

            end

            @. ForwardEuler.df_tmp = ForwardEuler.df / ForwardEuler.f
            Cr = -minimum(filter(isfinite,ForwardEuler.df_tmp)) 

        end   

        if Verbose == 1 && Cr > 1.0
            println("Cr = $Cr, t=$t, dt=$dt, system may be unstable")
        elseif Verbose == 2
            println("\rCr = $Cr, t=$t, dt=$dt")
        elseif Verbose == 3
            println("Cr = $Cr,Cr_Bin = $Cr_Bin, Cr_Emi = $Cr_Emi, Cr_Flux = $Cr_Flux, t=$t, dt=$dt")
        end
    end
    
    # update state vector f
    @. ForwardEuler.f += ForwardEuler.df
    # removing negative values (values less than 1f-28 for better stability)
    @. ForwardEuler.f = ForwardEuler.f*(ForwardEuler.f>=1f-30)
    # hacky fix for inf values
    @. ForwardEuler.f = ForwardEuler.f*(ForwardEuler.f!=Inf)


end

"""
    ForwardSymplecticEuler(dt0,dt,t,Verbose)

Forward Symplectic (Semi-Implicit) Euler time-stepping method for the Boltzmann equation. Symplectic integrator updates momentum space first then physical space.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardSymplecticEulerStruct)(dt0,dt,t,Verbose::Int64)

    # scaling of timestepping

        dt_scale = dt / dt0

    # Add injection term

        @. method.f_tmp = method.f + method.df_Inj * dt_scale

    # update momentum space using f at time t

        # create df_PFlux due to momentum flux terms
        mul!(method.df_PFlux,method.P_Flux,method.f_tmp)
        @. method.df_Momentum = -method.df_PFlux # minus sign as flux terms are on RHS of Boltzmann equation, also resets df_Momentum
        if !isfinite(sum(method.df_PFlux)) #TODO: speed up using find first non-finite value instead of summing entire array
            println("non-finite value in df_PFlux calculation, $(sum(method.df_PFlux))")
        end

        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f_tmp)
            @. method.df_Momentum += method.df_Emi
            if !isfinite(sum(method.df_Emi))
                println("non-finite value in df_Emi calculation, $(sum(method.df_Emi))")
            end
        end
            
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method)
            @. method.df_Momentum += method.df_Bin
            if !isfinite(sum(method.df_Bin))
                println("non-finite value in df_Bin calculation, $(sum(method.df_Bin))")
            end
        end

        @. method.df_Momentum *= method.invAp_Flux # Assumes Ap_flux is diagonal and stored as a vector

        if !isfinite(sum(method.df_Momentum))
            println("non-finite value in momentum df calculation")
        end

        # update f_momentum (f after momentum update)
        @. method.f_Momentum = method.f_tmp + method.df_Momentum * dt_scale
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Momentum = method.f_Momentum*(method.f_Momentum>=0f0)

    # update physical space using f_Momentum (f after momentum update)

        # create df_XFlux due to space flux terms
        mul!(method.df_XFlux,method.X_Flux,method.f_Momentum)
        @. method.df_Space = -method.df_XFlux # minus sign as flux terms are on RHS of Boltzmann equation, also resets df_Space
        if !isfinite(sum(method.df_XFlux))
            println("non-finite value in df_XFlux calculation, $(sum(method.df_XFlux))")
        end

        @. method.df_Space *= method.invAp_Flux # Assumes Ap_flux is diagonal and stored as a vector

        if !isfinite(sum(method.df_Space))
            println("non-finite value in space df calculation")
        end

        # update f_Space (f after momentum and space updates)
        @. method.f_Space = method.f_Momentum + method.df_Space * dt_scale
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Space = method.f_Space*(method.f_Space>=1f-28)


    # CFL condition check TODO: add adaptive time stepping based on CFL condition 

        if Verbose == 1 || Verbose == 2 || Verbose == 3

            Cr = 0.0
            Cr_Bin = 0.0
            Cr_Emi = 0.0
            Cr_PFlux = 0.0
            Cr_Momentum = 0.0
            Cr_XFlux = 0.0


            # Cr (CFL) condition check
            if sum(method.f) != 0.0

                if Verbose == 3
            
                    # Binary CFL
                    if method.Binary_Interactions
                        @. method.df_tmp = method.df_Bin * method.invAp_Flux / method.f_tmp * dt_scale
                        Cr_Bin = -minimum(filter(isfinite,method.df_tmp))
                    end

                    # Emission CFL
                    if method.Emission_Interactions
                        @. method.df_tmp = method.df_Emi * method.invAp_Flux / method.f_tmp * dt_scale
                        Cr_Emi = -minimum(filter(isfinite,method.df_tmp))
                    end

                    # P Flux CFL
                    @. method.df_tmp = -method.df_PFlux * method.invAp_Flux / method.f_tmp * dt_scale
                    Cr_PFlux = -minimum(filter(isfinite,method.df_tmp))

                    # Momentum CFL
                    @. method.df_tmp = (method.df_Bin + method.df_Emi - method.df_PFlux) * method.invAp_Flux / method.f_tmp * dt_scale
                    Cr_Momentum = -minimum(filter(isfinite,method.df_tmp))

                    # X Flux CFL
                    @. method.df_tmp = -method.df_XFlux * method.invAp_Flux / method.f_Momentum * dt_scale
                    Cr_XFlux = -minimum(filter(isfinite,method.df_tmp))

                end

                # Cr is calculated for the entire time step (momentum and space updates)
                # f_tmp - f = df of the momentum step
                @. method.df_tmp = (method.df_Space + method.df_Momentum) / method.f_tmp * dt / dt0 
                Cr = -minimum(filter(isfinite,method.df_tmp)) 

            end   

            if Verbose == 1 && Cr > 1.0
                println("Cr = $(round(Cr, sigdigits=3)), t=$t, dt=$dt, system may be unstable")
            elseif Verbose == 2
                println("\rCr = $(round(Cr, sigdigits=3)), t=$t, dt=$dt")
            elseif Verbose == 3
                println("Cr = $(round(Cr, sigdigits=3)), Cr_Bin = $(round(Cr_Bin, sigdigits=3)), Cr_Emi = $(round(Cr_Emi, sigdigits=3)), Cr_PFlux = $(round(Cr_PFlux, sigdigits=3)), Cr_Momentum = $(round(Cr_Momentum, sigdigits=3)), Cr_X_Flux = $(round(Cr_XFlux, sigdigits=3)),  t=$t, dt=$dt")
            end
        end

    # update state vector f with momentum, space and injection updates
    
        @. method.f = method.f_Space
        # removing negative values (values less than 1f-28 for better stability)
        @. method.f = method.f*(method.f>=1f-10)
        # hacky fix for inf values
        @. method.f = method.f*(method.f!=Inf)

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

function update_Big_Bin!(method::SteppingMethodType)
    
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

    @assert size(method.M_Bin) == (n_momentum^2,n_momentum) "M_Bin is not the correct size"

    f = method.f

    Vol = method.Vol

    Domain = method.Bin_Domain

    # Thanks to Emma Godden for fixing a bug here
    #temp = reshape(method.M_Bin_Mul_Step,length(f)*length(f))

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1 # starts at 0

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view f[start_idx:end_idx]

        if isnothing(Domain) || in(off_space,Domain) || isnothing(findfirst(!iszero,fView))

            df_BinView = @view method.df_Bin[start_idx:end_idx]
            mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape
            mul!(df_BinView,method.M_Bin_Mul_Step,fView)

            #println("M_Bin = $(sum(method.M_Bin))")
            #println("M_Bin_Mul_Step_reshape = $(sum(method.M_Bin_Mul_Step_reshape))")
            #println("M_Bin_Mul_Step = $(sum(method.M_Bin_Mul_Step))")

            # multiply by volume element
            df_BinView .*= view(Vol,off_space+1)

        else
            continue
        end

    end

    #=
    # TODO: Update Implicit
    if method.Implicit
        # assign jacobian elements
        @. method.Jac += 2*method.M_Bin_Mul_Step*view(Vol,off_space+1)
    end=#

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