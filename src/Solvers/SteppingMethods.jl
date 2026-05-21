"""
    ForwardEuler(t_start,t_stop,dt,Verbose)

Forward (Explicit) Euler time-stepping method for the transport equation.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # update momentum and space space using f at time t

    mul!(method.df_Flux,method.F_Flux,method.f)
    @. method.df = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # create df_Emi due to emission terms
    if method.Emission_Interactions
        mul!(method.df_Emi,method.M_Emi,method.f)
        @. method.df += method.df_Emi
    end
        
    # create df_Bin due to binary interactions
    if method.Binary_Interactions
        update_Big_Bin!(method,method.f)
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

    # remove df that are small compared to f
    #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)

    # Cr (CFL) condition check
    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
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
                Cr_Bin = -minimum(method.df_tmp)
            end

            # Emission CFL
            if method.Emission_Interactions
                @. method.df_tmp = method.df_Emi / method.f 
                @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
                Cr_Emi = -minimum(method.df_tmp)
            end

            # Flux CFL
            @. method.df_tmp = -method.df_Flux / method.f 
            @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
            Cr_Flux = -minimum(method.df_tmp)

        end   

    end

    # Adaptive timestep control based on CFL (Courant) number `Cr`.
    # Strategy: be conservative when increasing dt (small, gradual growth)
    # but allow immediate reduction when Cr rises. This avoids oscillation
    # caused by aggressive grow/shrink symmetry.
    dt_old = dt
    target_Cr = 1f0        # desired Courant number
    safety = 0.8f0         # safety factor
    max_growth = 1.1f0     # hard cap on growth per step
    max_shrink = 0.67f0    # hard cap on shrink per step
    gentle_growth = 1.1f0  # soft cap for routine increases
    grow_trigger_ratio = 1.1e0  # require Cr_prev/Cr >= this to allow growth

    if !(isfinite(Cr) && Cr > 0.0)
        adaptive_factor = 1.0f0
    else
        raw_factor = safety * (target_Cr / max(Cr, eps(eltype(Cr))))
        clamped = clamp(raw_factor, max_shrink, max_growth)

        if raw_factor > 1.0f0 
            # growth: only if Cr has dropped sufficiently compared to previous Cr
            if method.step == 1
                # no history: do not grow on first step
                adaptive_factor = 1.0f0
            else
                # reduce growth factor smoothly as Cr approaches target: when Cr ≈ safety*target_Cr, factor → 1.0
                target_threshold = safety * target_Cr
                reduction_factor = min(1.0f0, Cr / max(target_threshold, eps(eltype(Cr))))
                clamped = 1.0f0 + (clamped - 1.0f0) * (1.0f0 - reduction_factor)
                adaptive_factor = clamped
                println("1: ", raw_factor, " ", adaptive_factor)
            end
        else
            # shrink: apply immediately (within hard cap)
            adaptive_factor = (Cr_prev+safety*target_Cr)/2Cr#raw_factor
            println("2: ", raw_factor, " ", adaptive_factor)
        end
    end

    if method.step == 1
        adaptive_factor = 1.0
    end

    #=if method.step > 2
    if (Cr < 1e-3)
        adaptive_factor = min(1.1,(Cr_prev/Cr))
        println("5: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    elseif Cr > Cr_prev
        adaptive_factor = max(0.67,0.5Cr_prev/Cr)
        println("7: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    else
        adaptive_factor = 1.0
        println("7: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    end
    end=#

    adaptive_factor = 1.0

    # apply adaptive factor and update scaling
    dt = dt * adaptive_factor
    dt_scale = dt / dt0

    if method.step == 1
        method.Cr = Cr
    else
        method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
    end

    # will we reached the next t_save?

    t_next = t_start + dt
    if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
        save = true
    elseif t_next >= t_stop
        adaptive_factor *= (t_stop - t_start) / dt # adjust adaptive factor for final time step to ensure we end exactly at t_stop 
        dt = t_stop - t_start
        save = true
    else
        save = false
    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)), Cr_adapted = $(round(Cr*adaptive_factor,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # remove df that are small compared to f again after scaling time
    @. method.df *= adaptive_factor
    #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
    #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

    #@. method.df_tmp = method.df / method.f
    #@. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
    #println(-minimum(method.df_tmp) )

    # update state vector f 
    @. method.f += method.df  

    # remove masked off domain regions
    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    # removing negative values (values less than 1f-28 for better stability)
    @. method.f = method.f * (method.f>method.n_cut) * sign(method.f)
    # hacky fix for inf values
    @. method.f = method.f * (method.f!=Inf)

    return dt,save

end

"""
    ForwardSymplecticEuler(t_start,t_stop,dt,Verbose)

Forward Symplectic (Semi-Implicit) Euler time-stepping method for the transport equation. Symplectic integrator updates momentum space first then physical space.

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
        @. method.df_Momentum = -method.df_PFlux # minus sign as flux terms are on RHS of transport equation, also resets df_Momentum

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
        @. method.df_Space = -method.df_XFlux * method.invAp_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Space

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
    Heun(t_start,t_stop,dt,Verbose)

Heun's time-stepping method for the transport equation.
"""
function (method::HeunStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # first step

        mul!(method.df_Flux,method.F_Flux,method.f)
        @. method.df_step = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux
        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f)
            @. method.df_step += method.df_Emi
        end
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f)
            @. method.df_step += method.df_Bin
        end
        @. method.df_step *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
        # add injection term
        @. method.df_step += method.df_Inj * dt_scale 

        if !isfinite(sum(method.df_step))
            println("non-finite value in df_step calculation")
        end
        # update state vector f_step 
        @. method.f_step += method.df_step  

        # removing negative values

        @. method.f_step = method.f_step * (method.f_step>method.n_cut) * sign(method.f_step)

    # second step

        mul!(method.df_Flux,method.F_Flux,method.f_step)
        @. method.df = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux
        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f_step)
            @. method.df += method.df_Emi
        end
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f_step)
            @. method.df += method.df_Bin
        end
        @. method.df *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
        # add injection term
        @. method.df += method.df_Inj * dt_scale 

        if !isfinite(sum(method.df))
            println("non-finite value in df calculation")
        end

    # update state vector df for both steps
    
        @. method.df += method.df_step

        # mask of df regions
        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        if !isfinite(sum(method.df))
            println("non-finite value in df calculation")
        end

    # Cr (CFL) condition check

        Cr = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=method.n_cut, Inf,method.df / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
            Cr = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        if method.step == 1
            method.Cr = Cr
        else
            method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        end

        # remove df that are small compared to f again after scaling time
        @. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

    # will we reached the next t_save?

        t_next = t_start + dt
        if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
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

    # output verbose information 

        if Verbose == 1 && Cr > 1.0
            println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
        elseif Verbose == 2
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)), Cr_adapted = $(round(Cr*adaptive_factor,sigdigits=3))")
        elseif Verbose == 3
            println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        end
        if Verbose > 0
            flush(stdout)
        end

    # update state vector f 

        two = eltype(method.f)(2)
        @. method.f += method.df / two

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing negative values (values less than 1f-28 for better stability)

        @. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

"""
    MPE(t_start,t_stop,dt,Verbose)

Modified Patankar Euler time-stepping method for the transport equation.
"""
function (method::MPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # Create Q matrix for Patankar update

        println("1")
        Q = copy(method.Q) # method.Q is (F+MEmi) with structural zeros for binary interactions 
        println("max diag: ",maximum(diag(Q)))
        println("min diag: ",minimum(diag(Q)))
        println(length(Q.nzval))
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))
        println("2")
        # add M_Bin_Mul_Step to Q due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f,Q)
        end # now Q = (F+MEmi+MBin)
        println(length(Q.nzval))
        println(maximum(abs.(Q)))
        println("max diag: ",maximum(diag(Q)))
        println("min diag: ",minimum(diag(Q)))
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))

    # set up Ay=b system (where A=Q, y=f(t+1), b=A^-f(t)+df_Inj)

        println("4")
        Q = method.invA_Flux * Q # Assumes Ap_flux is diagonal and stored as a vector, also scales Q by dt_scale
        println(length(Q.nzval)) 
        
        println("5")
        @. Q.nzval *= dt_scale

        println(maximum(diag(Q)))

        println("6")
        Q = I - Q 
        println(length(Q.nzval))
        
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))

        println("7")
        method.f_tmp .= method.f + method.df_Inj .* dt_scale #+ method.invA_Flux * (method.MEmi * method.f) .* dt_scale

    # do inversion, gives f(t+1) in f_tmp using GMRES with ILU preconditioner

        println(maximum(diag(Q)))
        println("8")
        Qpre = ilu(Q)
        (x,stats) = gmres(Q,method.f_tmp, M=Qpre,ldiv=true,rtol=eps(eltype(Q)),atol=eps(eltype(Q)),restart=true)
        display(stats)
        method.f_tmp .= x

        #Qfact = lu(Q)
        #ldiv!(Qfact,method.f_tmp)

    # mask of df regions

        println("9")
        @. method.df = method.f_tmp - method.f 

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check

        Cr = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=method.n_cut, Inf,method.df / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
            Cr = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        if method.step == 1
            method.Cr = Cr
        else
            method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        end

        # remove df that are small compared to f again after scaling time
        @. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

    # will we reached the next t_save?

        t_next = t_start + dt
        if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
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

    # output verbose information 

        if Verbose == 1 && Cr > 1.0
            println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
        elseif Verbose == 2
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)), Cr_adapted = $(round(Cr*adaptive_factor,sigdigits=3))")
        elseif Verbose == 3
            println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        end
        if Verbose > 0
            flush(stdout)
        end

    # update state vector f 

        if !isnothing(method.df_mask)
            @. method.f = (method.f_tmp * (method.df_mask)) + method.f * (1-method.df_mask)
        else
            @. method.f = method.f_tmp
        end
        #@. method.f += method.df

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        #@. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

"""
    SymplecticMPE(t_start,t_stop,dt,Verbose)

Symplectic Modified Patankar Euler time-stepping method for the transport equation. Momentum transport is updated first using the implicit MPE scheme then space transport is performed using simple Explicit Forward Euler. 
"""
function (method::SymplecticMPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    # Create Q matrix for Patankar update and solve for momentum transport (includes injection)

        updateMomentumPatankarEuler!(method,method.f,dt_scale)
        if !isnothing(method.df_mask)
            @. method.f_step = (method.f_step * (method.df_mask)) + method.f * (1-method.df_mask)
        end 

    # space transport

        #@. method.f_step = method.f + method.df_Inj * dt_scale # add injection term to f_step for space transport
        mul!(method.df,method.X_Flux,method.f_step)
        mul!(method.df_Flux,method.invA_Flux,method.df,-dt_scale,zero(eltype(method.df))) # minus as X_Flux on RHS of transport equation

        @. method.df = method.df_Flux

    # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check

        Cr = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=0.0, Inf,(method.df+method.f_step- method.f) / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        if method.step == 1
            method.Cr = Cr
        else
            method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        end

        # remove df that are small compared to f again after scaling time
        @. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

    # will we reached the next t_save?

        t_next = t_start + dt
        if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
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

    # output verbose information 

        if Verbose == 1 && Cr > 1.0
            println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
        elseif Verbose == 2
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)), Cr_adapted = $(round(Cr*adaptive_factor,sigdigits=3))")
        elseif Verbose == 3
            println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        end
        if Verbose > 0
            flush(stdout)
        end

    # update state vector f 

        #println(sum(method.f[1:288]), " ", sum(method.df[1:288]), " ", sum(method.f_step[1:288]))
        #println(sum(method.f[289:888]), " ", sum(method.df[289:888]), " ", sum(method.f_step[289:888]))
        println(sum(method.f), " ", sum(method.df), " ", sum(method.f_step))

        @. method.f = method.df + method.f_step

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        #@. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

"""
    BackwardsEuler(dg,g,t,dt)

Backwards (Implicit) Euler time-stepping method for the transport equation. 

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
    @. BackwardEuler.df -= BackwardEuler.df_Flux # minus sign as flux terms are on RHS of transport equation
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
    #@. Euler.df -= Euler.df_Flux # minus sign as flux terms are on RHS of transport equation

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

function update_Big_Bin!(method::ExplicitSteppingMethod,f::AbstractVector)
    
    n_momentum=size(method.M_Bin,2)

    for off_space in method.Bin_Domain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view(f[start_idx:end_idx])

        df_BinView = @view method.df_Bin[start_idx:end_idx]
        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        @inbounds vol = method.Vol[off_space+1]
        mul!(df_BinView,method.M_Bin_Mul_Step,fView,vol,zero(eltype(method.df_Bin)))

    end

    return nothing

end

function update_Big_Bin!(method::ImplicitSteppingMethod,f::AbstractVector,Q::AbstractSparseMatrix)
    
    n_momentum=size(method.M_Bin,2)

    for off_space in method.Bin_Domain

        #println("Updating Bin for space index: ", off_space)
        
        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view(f[start_idx:end_idx])
        #println("Bin 1")

        @inbounds vol = method.Vol[off_space+1]
        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(eltype(fView))) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        #println("Bin 2")
        @. @view(Q[start_idx:end_idx,start_idx:end_idx]) += method.M_Bin_Mul_Step

    end

    return nothing

end

function updateMomentumPatankarEuler!(method::SymplecticMPEStruct,f::AbstractVector{T},dt_scale::T) where T

    Spacetime = method.PhaseSpace.Spacetime
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    n_space = x_num*y_num*z_num
    
    n_momentum=size(method.Q,2)

    I_matrix = Diagonal(ones(method.Precision, n_momentum))
    Dl = Diagonal(ones(Float64, n_momentum))
    Dr = Diagonal(ones(Float64, n_momentum))
    Qs = zeros(Float64, size(method.Q))
    Qtmp = zeros(Float64, size(method.Q))
    Qflux = zeros(Float64, size(method.Q))
    Qbin = zeros(Float64, size(method.Q))
    fstep = zeros(eltype(method.Precision), n_momentum)
    df_Inj = similar(fstep)
    ftmp = similar(fstep)
    fs = zeros(Float64,size(fstep))
    sigma = Diagonal(zeros(Float64, n_momentum))
    tmpdiag = Diagonal(zeros(Float64, n_momentum))

    for off_space in 0:n_space-1 

        fill!(method.Q,zero(method.Precision)) # reset Q

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        P_Flux_View = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx])
        df_Inj .= @view(method.df_Inj[start_idx:end_idx])
        @. sigma.diag = ifelse(fstep > 0.0, 1.0,0.0)

        if method.Emission_Interactions
            M_Emi_View = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
            @. Qflux = P_Flux_View #- M_Emi_View # plus as P_Flux on RHS of transport equation and so is Q, therefore minus on M_Emi which is on LHS
        else
            @. Qflux = P_Flux_View # plus as P_Flux on RHS of transport equation and so is Q
        end
        
        if method.Binary_Interactions && off_space in method.Bin_Domain
        
            @inbounds vol = method.Vol[off_space+1]

            #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            # loss
            fill!(Qbin,zero(eltype(Qbin))) # reset Qbin
            mul!(method.M_Bin_Mul_Step_reshape,method.Liij,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            mul!(tmpdiag.diag,method.M_Bin_Mul_Step,fstep)

            @. tmpdiag.diag *= ifelse(fstep > 0.0, 1/fstep,0.0) 
            Qbin .-= tmpdiag # minus as M_Bin is on LHS of transport equation and Q is on RHS

            # gain 
            mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            Qbin .-= method.M_Bin_Mul_Step * sigma# minus as M_Bin is on LHS of transport equation and Q is on RHS
            
        end

        invA_Flux_View = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx])

        Qtmp .= invA_Flux_View * (Qbin + Qflux * sigma) * dt_scale

        println("dt scale: $dt_scale")
        colsum = sum(Qtmp,dims=1)
        rowsum = sum(Qtmp,dims=2)
        println("max col sum: ", maximum(colsum), " max row sum: ", maximum(rowsum), " min col sum: ", minimum(colsum), " min row sum: ", minimum(rowsum))

        cutoff_offdiagonal_to_diag_relative!(Qtmp, 1e-6)

        println("sum diag: ", sum(diag(Qtmp)), " sum off diag: ", sum(Qtmp) - sum(diag(Qtmp)))
        display(sparse(Qtmp))

        #fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        Qtmp .+= I_matrix # add identity for Patankar update

        kwo = cond(Qtmp)
        println("off space: $off_space, condition without equilibration: $kwo")

        norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=3)
        Qs .= Dl * Qtmp * Dr

        kwo = cond(Qs)
        println("off space: $off_space, condition with equilibration: $kwo")

        if kwo > sqrt(1/eps(eltype(Qtmp)))
            k = kwo
            step = zero(eltype(method.Precision))
            fstep .= fstep
            dstep = one(eltype(method.Precision))
            a = one(eltype(method.Precision))
            while step < one(eltype(method.Precision))
                k = method.Precision(Inf)
                while k > sqrt(1/eps(eltype(Qtmp)))

                    s = svdvals(Qtmp)
                    println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

                    if method.Binary_Interactions && off_space in method.Bin_Domain   
            
                        @inbounds vol = method.Vol[off_space+1]

                        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                        @. Qbin = -method.M_Bin_Mul_Step # minus as M_Bin is on LHS of transport equation and Q is on RHS

                        Qtmp .= invA_Flux_View * (Qbin + Qflux) * dt_scale * a + I_matrix # add identity for Patankar update
                
                    else

                        Qtmp .= invA_Flux_View * Qflux * dt_scale * a + I_matrix # add identity for Patankar update

                    end

                    norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=3)
                    Qs .= Dl * Qtmp * Dr

                    k = cond(Qs)

                    println("off space: $off_space, cond without: $kwo, condition with right and left: $k, a = $a, step = $step, dstep = $dstep")

                    if k > sqrt(1/eps(eltype(Qtmp)))
                        a *= method.Precision(1/2)
                        dstep *= method.Precision(1/2)
                    end
                    
                end

                if method.Emission_Interactions
                    ftmp .= fstep + df_Inj .* dt_scale * a + invA_Flux_View * M_Emi_View * fstep .* dt_scale * a # add M_Emi contribution to ftmp for Patankar update
                else
                    ftmp .= fstep + df_Inj .* dt_scale * a
                end

                s = svdvals(Qs)
                println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

                fs .= Dl * ftmp
                    println("debug ftmp: sum=", sum(ftmp), " max=", maximum(ftmp), " sum(Dl*ftmp)=", sum(Dl * ftmp))

                QLU = lu(Qs)
                ldiv!(QLU,fs)
                fstep .= Dr * fs
                # diagnostics: check mass conservation and residual for this block
                mass_before = sum(@view(f[start_idx:end_idx]))
                mass_after = sum(fstep)
                res = norm(Qs * fs - Dl * ftmp) / max(one(eltype(fs)), norm(Dl * ftmp))
                println("off_space=$off_space mass_before=", mass_before, " mass_after=", mass_after, " Δ=", mass_after-mass_before, " res=", res)
                @. sigma.diag = ifelse(fstep != 0.0, 1.0,0.0)

                #gmres!(method.GMRESWorkspace,Qs,fs;rtol = 1f-8,atol = 0f0)
                #fstep .= Dr * method.GMRESWorkspace.x

                step += dstep

            end

        else

            if method.Emission_Interactions
                ftmp .= fstep + df_Inj .* dt_scale + invA_Flux_View * M_Emi_View * fstep * dt_scale # add M_Emi contribution to ftmp for Patankar update
            else
                ftmp .= fstep + df_Inj .* dt_scale
            end

            #norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
            #Qs .= Dl * Qtmp * Dr

            #k = cond(Qs)

            #println("off space: $off_space, cond without: $kwo, condition with right and left: $k")

            s = svdvals(Qs)
            println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")
             
            fs .= Dl * ftmp
                println("debug ftmp: sum=", sum(ftmp), " max=", maximum(ftmp), " sum(Dl*ftmp)=", sum(Dl * ftmp))

            QLU = lu(Qs)
            ldiv!(QLU,fs)
            fstep .= Dr * fs

            # diagnostics: check mass conservation and residual for this block
            mass_before = sum(@view(f[start_idx:end_idx]))
            mass_after = sum(fstep)
            res = norm(Qs * fs - Dl * ftmp) / max(one(eltype(fs)), norm(Dl * ftmp))
            println("off_space=$off_space mass_before=", mass_before, " mass_after=", mass_after, " Δ=", mass_after-mass_before, " res=", res)

            #gmres!(method.GMRESWorkspace,Qtmp,ftmp;rtol = 1f-10,atol = 0f0)
            #fstep .= method.GMRESWorkspace.x

        end

        @view(method.f_step[start_idx:end_idx]) .= fstep

        # GMRES solve for momentum update
        #lu!(method.QLU,method.Q)
        #Qs = Dl * method.Q * Dr
        #fs = Dl * f_tmpView
        #gmres!(method.GMRESWorkspace,Qs,fs,fs;rtol = 1f-6,atol = 0f0,itmax = n_momentum,history = true)

        #relres = norm(method.Q * Dr * method.GMRESWorkspace.x - f_tmpView) / norm(f_tmpView)
        #println("off_space $off_space residual = ", relres)

        #@view(method.f_step[start_idx:end_idx]) .= Dr * method.GMRESWorkspace.x

        #gmres!(method.GMRESWorkspace,method.Q,f_tmpView;M=Dl,N=Dr,ldiv=false)
        #pre = diag(method.Q) # diagonal preconditioner
        #pre = Diagonal(pre)
        #gmres!(method.GMRESWorkspace,method.Q,f_tmpView)

        #relres = norm(method.Q*method.GMRESWorkspace.x - f_tmpView) / norm(f_tmpView)
        #println("off_space $off_space residual = ", relres)

        #display(method.GMRESWorkspace.stats)
        #@view(method.f_step[start_idx:end_idx]) .= method.GMRESWorkspace.x

        # direct LU solve for momentum update
        #lu!(method.QLU,method.Q)
        #ldiv!(QLU,f_tmpView)
        #@view(method.f_step[start_idx:end_idx]) .= f_tmpView

        # direct QR solve for momentum update
        #QQR = qr(method.Q)
        #ldiv!(QQR,f_tmpView)
        #@view(method.f_step[start_idx:end_idx]) .= f_tmpView

        # direct SVD solve for momentum update
        #QSVD = svd(method.Q)
        #ldiv!(QSVD,f_tmpView)

    end

    return nothing

end

function pow2_equilibration_matrices(A::AbstractMatrix{T}; iters=3) where {T<:AbstractFloat}
    n, m = size(A)

    dl = ones(T, n)
    dr = ones(T, m)

    As = copy(A)

    for k in 1:iters
        # Left scaling: rows
        for i in 1:n
            s = maximum(abs, @view As[i, :])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[i, :] .*= α
                dl[i] *= α
            end
        end

        # Right scaling: columns
        for j in 1:m
            s = maximum(abs, @view As[:, j])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[:, j] .*= α
                dr[j] *= α
            end
        end
    end

    return Diagonal(dl), Diagonal(dr)
end

function pow2_equilibration_matrices!(Dl::Diagonal{T,Vector{T}}, Dr::Diagonal{T,Vector{T}}, A::AbstractMatrix{T}; iters=3) where {T<:AbstractFloat}

    n, m = size(A)

    Dl.diag .= ones(T)
    Dr.diag .= ones(T)

    As = copy(A)

    for k in 1:iters
        # Left scaling: rows
        for i in 1:n
            s = maximum(abs, @view As[i, :])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[i, :] .*= α
                Dl.diag[i] *= α
            end
        end

        # Right scaling: columns
        for j in 1:m
            s = maximum(abs, @view As[:, j])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[:, j] .*= α
                Dr.diag[j] *= α
            end
        end
    end

    return nothing
end

function norm_equilibration_matrices(A::AbstractMatrix{T}; iters=5) where {T<:AbstractFloat}
    n, m = size(A)

    dl = ones(T, n)
    dr = ones(T, m)

    As = copy(A)

    for k in 1:iters
        # Row norm scaling
        for i in 1:n
            s = norm(@view As[i, :])
            if s > zero(T)
                α = one(T) / T(s)
                As[i, :] .*= α
                dl[i] *= α
            end
        end

        # Column norm scaling
        for j in 1:m
            s = norm(@view As[:, j])
            if s > zero(T)
                α = one(T) / T(s)
                As[:, j] .*= α
                dr[j] *= α
            end
        end
    end

    return Diagonal(dl), Diagonal(dr)
end

function norm_equilibration_matrices!(Dl::Diagonal{T,Vector{T}}, Dr::Diagonal{T,Vector{T}}, A::AbstractMatrix{T}; iters=5) where {T<:AbstractFloat}
    n, m = size(A)

    Dl.diag .= ones(T)
    Dr.diag .= ones(T)
    As = copy(A)

    for k in 1:iters
        # Row norm scaling
        for i in 1:n
            s = norm(@view As[i, :])
            if s > zero(T)
                α = one(T) / T(s)
                As[i, :] .*= α
                Dl.diag[i] *= α
            end
        end

        # Column norm scaling
        for j in 1:m
            s = norm(@view As[:, j])
            if s > zero(T)
                α = one(T) / T(s)
                As[:, j] .*= α
                Dr.diag[j] *= α
            end
        end
    end

    return nothing

end

function cutoff_offdiagonal_to_diag_relative!(A, frac)
    n, m = size(A)
    @assert n == m "A must be square"

    for j in 1:n
        colmax = maximum(abs(A[i, j]) for i in 1:n if i != j)
        cutoff = eps(colmax)

        for i in 1:n
            if i != j && abs(A[i, j]) < cutoff
                removed = A[i, j]
                A[i, j] = zero(eltype(A))
                A[j, j] -= removed
            end
        end
    end

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

# to be removed in julia 1.13#

function generic_lufact!(A::AbstractMatrix{T}, pivot::Union{RowMaximum,NoPivot,RowNonZero} = lupivottype(T), ipiv::AbstractVector{LinearAlgebra.BlasInt} = Vector{LinearAlgebra.BlasInt}(undef,min(size(A)...));check::Bool = true, allowsingular::Bool = false) where {T}
    
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
