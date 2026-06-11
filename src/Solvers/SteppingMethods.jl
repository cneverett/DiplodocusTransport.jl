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

    @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    # Create Q matrix for Patankar update and solve for momentum transport (includes injection)

        updateMomentumPatankarEuler!(method,method.f,dt_scale)
        if !isnothing(method.df_mask)
            @. method.f_step = (method.f_step * method.df_mask) + method.f * (1-method.df_mask)
        end 
        if !isnothing(method.f_mask)
            @. method.f_step *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f_step = ifelse(method.f_step<=method.n_cut,zero(eltype(method.f_step)),method.f_step)

    # Cr (CFL) condition check for momentum update

        Cr_momentum = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=0.0, Inf,(method.f_step - method.f) / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_momentum = -minimum(method.df_tmp)
        end

    # space transport

        #@. method.f_step = method.f + method.df_Inj * dt_scale # add injection term to f_step for space transport
        mul!(method.df_Flux,method.X_Flux,method.f_step)
        method.df .= -method.invA_Flux * method.df_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check space update

        Cr_space = 0.0
        dt_old = dt
        sum_f = sum(method.f_step)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f_step<=0.0, Inf,(method.df / method.f_step))
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_space = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        #if method.step == 1
        #    method.Cr = Cr
        #else
        #    method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        #end

        # remove df that are small compared to f again after scaling time
        #@. method.df *= adaptive_factor
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
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr_momentum = $(round(Cr_momentum,sigdigits=3)), Cr_space = $(round(Cr_space,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
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

function (method::SymplecticSymmetricMPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    # Create Q matrix for Patankar update and solve for momentum transport (includes injection)

        updateMomentumPatankarEulerSymmetric!(method,method.f,dt_scale)
        if !isnothing(method.df_mask)
            @. method.f_step = (method.f_step * method.df_mask) + method.f * (1-method.df_mask)
        end 
        if !isnothing(method.f_mask)
            @. method.f_step *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f_step = ifelse(method.f_step<=method.n_cut,zero(eltype(method.f_step)),method.f_step)

    # Cr (CFL) condition check for momentum update

        Cr_momentum = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=0.0, Inf,(method.f_step - method.f) / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_momentum = -minimum(method.df_tmp)
        end

    # space transport

        #@. method.f_step = method.f + method.df_Inj * dt_scale # add injection term to f_step for space transport
        mul!(method.df_Flux,method.X_Flux,method.f_step)
        method.df .= -method.invA_Flux * method.df_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check space update

        Cr_space = 0.0
        dt_old = dt
        sum_f = sum(method.f_step)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f_step<=0.0, Inf,(method.df / method.f_step))
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_space = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        #if method.step == 1
        #    method.Cr = Cr
        #else
        #    method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        #end

        # remove df that are small compared to f again after scaling time
        #@. method.df *= adaptive_factor
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
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr_momentum = $(round(Cr_momentum,sigdigits=3)), Cr_space = $(round(Cr_space,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
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
function (method::BackwardEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)
        
    # create df_Bin due to binary interactions (jacobians are added in update_Big_Bin! if implicit)
    if method.Binary_Interactions
        update_Big_Bin!(method,dt_scale)
    end


    
    # removing negative values
    @. method.fstep = method.fstep * (method.fstep>=method.n_cut)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    return dt,save

end

"""
    ES2(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. 

"""
function (method::ES2Struct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = zeros(eltype(method.f),length(method.f))

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # half momentum update

        #=mul!(ftmp,method.invImMP,method.fstep+method.df_Inj * dt_scale/2) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end =#
        
    # binary update

        if method.Binary_Interactions
            update_momentum!(method,dt_scale)
        end

    # half momentum update

        #=mul!(ftmp,method.invImMP,method.fstep+method.df_Inj * dt_scale/2) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end=#

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end

function update_momentum!(method::ES2Struct,dt_scale::T) where T

    Precision = method.Precision
    
    n_momentum=method.PhaseSpace.Grids.n_momentum
    n_space=method.PhaseSpace.Grids.n_space
    
    fold = zeros(Precision,n_momentum)
    fxdiv2= zeros(Precision,n_momentum)
    fx = zeros(Precision,n_momentum)
    fz = zeros(Precision,n_momentum)
    fnonzero = Diagonal(zeros(Precision,n_momentum))
    J = method.J

    n_species = length(method.PhaseSpace.name_list)

    M = zeros(Precision,(n_species+1),n_momentum)

    #LinearFlux = zeros(Precision,n_momentum,n_momentum)

    @view(M[1:n_species,:]) .= method.N 
    @view(M[n_species+1,:]) .= method.E

    #A = [ zeros(Precision,n_momentum,n_momentum) zeros(Precision,n_momentum) ; zeros(Precision,n_momentum)' zero(Precision) ]

    consv = zeros(n_species+1)
    consv[end] = one(Precision)

    for off_space in 0:n_space-1

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        fin = zeros(Precision,n_momentum)
        fout = zeros(Precision,n_momentum)
        df_Inj = @view(method.df_Inj[start_idx:end_idx])
        @inbounds vol = method.Vol[off_space+1]
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        A_Flux = @view(method.A_Flux[start_idx:end_idx,start_idx:end_idx])

        @. fnonzero.diag = ifelse(fold > zero(Precision), one(Precision),zero(Precision))

        #LinearFlux .= @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        #if method.Emission_Interactions
        #    LinearFlux .-= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        #end

        if method.Binary_Interactions && off_space in method.Bin_Domain

            fold .= @view(method.f[start_idx:end_idx]) 

            #if sum(fold) == Precision(0)
            #    continue
            #end

            α = Precision(1.0)
            accurate = false
            ΔE = zero(Precision)
            step = zero(Float64)
            dstep = one(Float64)
            ϵ = eps(Precision)

            while !accurate || step < 1.0

                #fold .+= df_Inj * α * dt_scale

            # linear flux update (should conserve energy directly) so do single exponential step
            #LinearFlux .= invA_Flux * LinearFlux * fnonzero

            #display(reshape(real.(eigvals(-LinearFlux .* dt_scale)) .* fnonzero,37,24))
            #Println("max 1/τ: ", maximum(-real.(eigvals(-LinearFlux .* dt_scale) .* fnonzero) ))
            #display(fnonzero)

            #fold .= exp(-LinearFlux .* dt_scale) * fold

            #@. fnonzero.diag = ifelse(fold > 0.0, 1.0,0.0)

            # Second order in time ES2 method
                # x/2 step

                    @time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fold

                    #s = zeros(Precision,n_momentum)
                    #D = Diagonal(s)
                    #s .= max.(abs.(fin), 0.5*dt_scale*α * abs.(method.M_Bin_Mul_Step * fin), ϵ)
                    #D = Diagonal(100 .*s)

                    #println("max D: ", maximum(D), " min D : ", minimum(D))
                    #println("max D: ", maximum(s), " min s : ", minimum(s))

                    #display(D)

                    #println("max J: ", maximum(method.M_Bin_Mul_Step), " min J : ", minimum(method.M_Bin_Mul_Step))

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])

                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    A = [J (D \ df_Inj) ; zeros(Precision,n_momentum+1)']

                    #println("max J: ", maximum(J), " min J : ", minimum(J))

                    #display(sparse(A))

                    fin .= D \ fin 

                    f = [fin ; one(Precision)]

                    #fxdiv2 .= exp(method.M_Bin_Mul_Step .* dt_scale) * fin
                    #fxdiv2 .= expv(0.5*dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)
                    f .= expv(dt_scale*α,A,f;tol=Precision(1e-12),m=500)

                    fxdiv2 .= D * @view(f[1:n_momentum])

                    #display(reshape(fxdiv2,37,24))
                    #display(reshape(df_Inj*α*dt_scale,37,24))

                    #F = schur(method.M_Bin_Mul_Step)
                    #w = F.Z' * fin
                    #fxdiv2 .= F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    #println(maximum(real.(F.values)))



                    #F = schur(J)
                    #w = F.Z' * fin
                    #fxdiv2 .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    #display(F)
                    #display(sum(F.values))
                    #println(maximum(real.(F.values)))

                    fxdiv2 .= invA_Flux * fxdiv2

                    @. fxdiv2 = ifelse(fxdiv2<=method.n_cut,zero(eltype(fxdiv2)),fxdiv2)
                
                # z step 

                    #=@time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fxdiv2,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fold 

                    s = max.(abs.(fin), dt_scale*α * abs.(method.M_Bin_Mul_Step*fin), ϵ)
                    D = Diagonal(s)

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])

                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    #fz .= expv(dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)

                    fin .= D \ fin

                    fz .= D * expv(dt_scale*α,J,fin;tol=Precision(1e-8),m=100)

                    #F = schur(J)
                    #w = F.Z' * fin
                    #fz .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    @. fz = ifelse(fz<=method.n_cut,zero(eltype(fz)),fz)

                    fz .= invA_Flux * fz=#

                # x step 

                    #=@time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fz,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fxdiv2

                    s = max.(abs.(fin), 0.5*dt_scale*α * abs.(method.M_Bin_Mul_Step*fin), ϵ)
                    D = Diagonal(s)

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])  
                    
                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    #fx .= expv(0.5*dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)

                    fin .= D \ fin

                    fx .= D * expv(0.5*dt_scale*α,J,fin;tol=Precision(1e-8),m=100)

                    #F = schur(J)
                    #w = F.Z' * fin
                    #fx .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    @. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)

                    fx .= invA_Flux * fx=#

                # combine 
               
                #@. fx = (fx + fz) / 2

                #@. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)

                @. fx = ifelse(fxdiv2<=method.n_cut,zero(eltype(fxdiv2)),fxdiv2)

                #println(sum(fx), " ", sum(fold), " ", sum(fx-fold), " ", sum(df_Inj * α * dt_scale) )

                #println(sum(fx[1:288]), " ", sum(fold[1:288]), " ", sum((fx-fold)[1:288]), " ", sum((df_Inj * α * dt_scale)[1:288]))

                # energy correction  
                ΔE = dot(method.E, fx) - dot(method.E, fold .+ df_Inj * α *dt_scale)

                Eold = dot(method.E, fold .+ df_Inj * α *dt_scale)
                Eerr = abs(dot(method.E, fx) / Eold - one(Precision))

                #=if Eerr < 1e-3
                    accurate = true
                    step += dstep * α
                    # error correction step to better improve accuracy``
                    #ϵ = eps(Float64(method.n_cut))
                    #Winv = method.E .* max.(fx, ϵ)
                    #q = Winv .* method.E
                    #q ./= dot(method.E, q)
                    #@. fx -= ΔE * q
  
                    #=ϵ = eps(Float64(method.n_cut))
                    Winv = method.E .* max.(fx, ϵ)
                    B = M * Diagonal(Winv) * M'
                    # precondition B with diagonal of B
                    D = Diagonal(1.0 ./ sqrt.(diag(B)))
                    Bs = D * B * D
                    q = Diagonal(Winv) * (M' * D * (Bs \ (D * consv)))
                    @. fx -= ΔE*q
                    @. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)=#
                    fold .= fx
                    println("step: ", step, " α: ", α, " Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fx))
                    if Eerr < 1e-4 
                        if step + dstep * 4.0 * α < 1.0 
                            # can fit two 2x steps
                            α *= 2.0
                        end
                    end 
                else
                    accurate = false
                    println("Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fx), " reducing α: ", α)
                    α *= 0.5
                end=#

                    ϵ = eps(Precision(method.n_cut))

                    Winv = #=method.E .*=# max.(fx, ϵ)

                    MW = M .* Winv'

                    L = MW * M'

                    display(M)
                    display(M * fx)
                    display([dot(method.N[1,:],fx) ; dot(method.N[2,:],fx) ; Eold])

                    r = M * fx - [dot(method.N[1,:],fx) ; dot(method.N[2,:],fx) ; Eold]
                    r = [0f0 ; 0f0 ; ΔE]

                    println("r: ", r, " ΔE: ", ΔE)

                    λ = L \ r

                    δu = Winv .* (M' * λ)

                    fx .-= δu

                    accurate = true
                    step += dstep * α

            end # while not accurate

            @. fstep = fx

        else #

            # linear flux update (should conserve energy directly) so do single exponential step
            #@view(A[1:n_momentum,1:n_momentum]) .= -invA_Flux * LinearFlux

            #fin = [fold ; one(Precision)]
            #fout = similar(fin)

            #fout .= expv(dt_scale,A,fin;tol=Precision(1e-6),m=50)

            #fstep .= @view(fout[1:n_momentum])

        end

    end

    return nothing

end

"""
    ExponentialRosenbrock(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. 

"""
function (method::ExponentialRosenbrockStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = zeros(eltype(method.f),length(method.f))

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end 

    # binary update

        update_momentum!(method,dt_scale)


    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end

function update_momentum!(method::ExponentialRosenbrockStruct,dt::T) where T

    Precision = method.Precision
    
    n_momentum=method.PhaseSpace.Grids.n_momentum
    n_space=method.PhaseSpace.Grids.n_space
    
    fold::Vector{Precision} = zeros(Precision,n_momentum)
    fout::Vector{Precision} = zeros(Precision,n_momentum)
    J = method.J
    Jtmp = copy(J)
    F::Vector{Precision} = zeros(Precision,n_momentum)
    J64 = zeros(Float64,n_momentum,n_momentum)
    F64 = zeros(Float64,n_momentum)
    Ftmp = copy(F)
    ϕ = zeros(Precision,n_momentum,2) # for storing ϕ functions
    D = Diagonal(1 ./ method.E)
    Dinv = Diagonal(method.E)
    #D.diag .= 1f0
    #Dinv.diag .= 1f0
    δ::Vector{Precision} = zeros(Precision,n_momentum)

    m = 100#200#30 # dimension of Krylov subspace for exponential action approximation, can adjust based on problem size and desired accuracy

    #Ks64s = KrylovSubspace{Float64}(n_momentum, m)
    #Ks32s = KrylovSubspace{Float64}(n_momentum, m)
    Ks64 = KrylovSubspace{Float64}(n_momentum, m)
    Ks32 = KrylovSubspace{Float64}(n_momentum, m)
    #Ks64n = KrylovSubspace{Float64}(n_momentum, m)
    #Ks32n = KrylovSubspace{Float64}(n_momentum, m)

    ϕcache = ExponentialUtilities.PhivCache(ϕ,m,1)

    ηtarget = 1f-3
    ϵ = 1f-20

    n_species = length(method.PhaseSpace.name_list)

    M = zeros(Precision,(n_species+1),n_momentum)

    #LinearFlux = zeros(Precision,n_momentum,n_momentum)

    @view(M[1:n_species,:]) .= method.N 
    @view(M[n_species+1,:]) .= method.E

    #A = [ zeros(Precision,n_momentum,n_momentum) zeros(Precision,n_momentum) ; zeros(Precision,n_momentum)' zero(Precision) ]

    consv = zeros(n_species+1)
    consv[end] = one(Precision)

    for off_space in method.ActiveDomain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        df_Inj = @view(method.df_Inj[start_idx:end_idx])

        if sum(df_Inj) == zero(Precision) && sum(fstep) == zero(Precision) 
            continue
        end

        @inbounds vol = method.Vol[off_space+1]
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        A_Flux = @view(method.A_Flux[start_idx:end_idx,start_idx:end_idx])
        M_Emi = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        P_Flux = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        invA = invA_Flux[1,1] # assumes the diagonal is just constant, should be true for a single spacial cell
        A = A_Flux[1,1]

        #LinearFlux .= @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        #if method.Emission_Interactions
        #    LinearFlux .-= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        #end

        if method.Binary_Interactions && off_space in method.Bin_Domain

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            while t < 1.0

                dtldt0 = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision)) 
                    J .= 2 .* method.M_Bin_Mul_Step 
                    J .+= M_Emi
                    #J .-= P_Flux
                    lmul!(dtldt0,J)
                    # Form F
                    method.M_Bin_Mul_Step .+= M_Emi
                    #method.M_Bin_Mul_Step .-= P_Flux
                    mul!(Ftmp,method.M_Bin_Mul_Step,fold)
                    Ftmp .+= df_Inj #* dtldt0
                    lmul!(A*dtldt0,Ftmp)
                    # remove small changes in J and F 
                    #threshold_Jrows!(J;rtol = eps(Precision))
                    #threshold_Jcols!(J;rtol = eps(Precision))
                    #J[abs.(J) .< 1f-20] .= zero(Precision) # remove small entries in D2 to improve stability of second step

                    #threshold_F!(F,fold,k;rtol=1f-4)

                    # Form Krylov Subspace
                    #Dn , Jn = balance_similarity(J; maxiter=5, tol=1e-3)
                    #println("cond Dn: ", cond(Dn))

                    #D.diag .= 1f0
                    mul!(Jtmp,J,D)
                    mul!(J,Dinv,Jtmp)
                    #@time J .= D \ J * D
                    mul!(F,Dinv,Ftmp)
                    #@time F .= D \ F
                    

                    #threshold_Jrows!(J;rtol = eps(Precision))
                    #threshold_Jcols!(J;rtol = eps(Precision))
                    #threshold_F!(F,D \ fold,k;rtol=1f-2)

                    #vals = abs.(J)

                    #println("minimum abs(nonzeros(J)): ", minimum(vals[vals .> 0]))
                    #println("maximum abs(nonzeros(J)): ", maximum(vals))

                    #Jn = Dn \ J * Dn
                    #Fn = Dn \ F

                    J64 .= Float64.(J)
                    F64 .= Float64.(F)

                    #F64n = Float64.(Fn)
                    #J64n = Float64.(Jn)
                    #println("cond J: ", cond(I -J))
     
                    arnoldi!(Ks64,J64,F64;m=m,reorthogonalize=true)
 
                    #arnoldi!(Ks32n, Jn, Fn;m=m)
                    #arnoldi!(Ks64n, J64n, F64n;m=m)

                    arnoldi!(Ks32, J, F;m=m,reorthogonalize=true)
                    #arnoldi!(Ks64, J64, F64;m=m)

                    #V32 = ExponentialUtilities.getV(Ks32)[:,1:end-1]
                    #H32 = ExponentialUtilities.getH(Ks32)[1:end-1,1:end]
                    V64 = ExponentialUtilities.getV(Ks64)[:,1:end-1]
                    H64 = ExponentialUtilities.getH(Ks64)[1:end-1,1:end]
                    #V32n = ExponentialUtilities.getV(Ks32n)[:,1:end-1]
                    #H32n = ExponentialUtilities.getH(Ks32n)[1:end-1,1:end]
                    #V64n = ExponentialUtilities.getV(Ks64n)[:,1:end-1]
                    #H64n = ExponentialUtilities.getH(Ks64n)[1:end-1,1:end]

                    #norm_scaled = norm(V32s' * V32s - I)
                    norm_unscaled = norm(V64' * V64 - I)
                    #norm_new = norm(V32n' * V32n - I)
                    #println("|V'V - I| Float64 scaled: ", norm(V64s' * V64s - I))  # TODO: scale switch between for when norm changes
                    #println("|V'V - I| Float32 scaled: ", norm_scaled)
                    #println("|V'V - I| Float64: ", norm(V64' * V64 - I))
                    println("|V'V - I| Float32: ", norm_unscaled)
                    #println("|V'V - I| Float64 balanced: ", norm(V64n' * V64n - I))
                    #println("|V'V - I| Float32 balanced: ", norm_new)
                    #println("cond((I - dt_scale*H64)) scaled: ", cond((I - dt*H64s)))
                    #println("cond((I - dt_scale*H32)) scaled: ", cond((I - dt*H32s)))
                    #println("cond((I - dt_scale*H64)): ", cond((I - dt*H64)))
                    println("cond((I - dt_scale*H64)): ", cond((I - dt*H64)))
                    #println("cond((I - dt_scale*H64)) balanced: ", cond((I - dt*H64n)))
                    #println("cond((I - dt_scale*H32)) balanced: ", cond((I - dt*H32n)))

                    if cond((I - dt*H64)) > 1f8 ||  norm(V64' * V64 - I) > 1f-5 
                        dt_local = dt_local * 0.5
                        continue
                    end
                    #cond_scaled = cond((I - dt*H64s))
                    #cond_unscaled = cond((I - dt*H64))

                    #if cond_scaled < cond_unscaled
                    #    println("Using scaled Krylov subspace for error estimation")
                    #    Ks = Ks64s
                    #    V = V64s
                    #    H = H64s
                    #    F = Fs
                    #    J = Js
                    #    D = D
                    #else
                        #println("Using unscaled Krylov subspace for error estimation")
                        #V = V64
                        #H = H64
                        #F = F
                        #J = J
                        #D.diag .= 1.0
                    #end
                    #println("size V : ", size(V))
                    #println("size H : ", size(H))

                    #println(abs(ExponentialUtilities.getH(Ks)[m+1,m]))
                    #β = norm(F)
                    #e1 = zeros(Precision,size(V,2))
                    #e1[1] = one(Precision)

                    # Compute φ functions of H
                    phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=false) # TODO: This allocates
                    mul!(δ,D,@view(ϕ[:,2]),k*invA,zero(Precision))
                    #@time δ .= invA_Flux * δ # TODO: this allocates
                    fout .= fold .+ δ
                    
                    #@. fstar = ifelse(fstar <= method.n_cut, zero(Precision),fstar)

                    #@. fstar = max(fstar, zero(Precision)) # remove negative values before error calculation, can cause issues with error estimation and stability of second step

                    # adaptive time stepping

                        # second step accuracy prediction

                        #=# form D2 with δ only need A(δ,δ)
                        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,δ,vol*dtldt0,zero(Precision))
                        mul!(D2,method.M_Bin_Mul_Step,δ)

                        D2 .= D \ D2

                        D2[abs.(D2) .< eps(Precision)] .= zero(Precision) # remove small entries in D2 to improve stability of second step=#

                        # compute second step ϕ functions
                        #=Ks2 = KrylovSubspace{Float64}(n_momentum, m)
                        arnoldi!(Ks2, J, D2;m=m) 
                        V2 = ExponentialUtilities.getV(Ks2)[:,1:end-1]
                        H2 = ExponentialUtilities.getH(Ks2)[1:end-1,1:end]
                        #println("cond((I - dt_old*H2)): ", cond((I - dt_old*H2)))
                        #println("norm(V2'V2 - I): ", norm(V2' * V2 - I))
                        β = norm(D2)
                        e1 = zeros(Precision,size(V2,2))
                        e1[1] = one(Precision)
                        ϕ3 = 2 * k * D * V2 * phiv(k,H2,β*e1,3)[:,4]=#
                        
                        # combine steps for final solution
                        
                        #@. fout = fold + δ #+ ϕ3

                        #@. fout = max(fout, zero(Precision)) # remove negative values before error calculation, can cause issues with error estimation and stability of second step

                        #println(minimum(fout), " ", maximum(fout))

                        #@. fout = max(fout, zero(Precision)) # remove negative values before error calculation, can cause issues with error estimation and stability of second step
                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        ΔE = dot(method.E, fout) - dot(method.E, fold .+ df_Inj * k*dtldt0)
                        Eold = dot(method.E, fold .+ df_Inj * k * dtldt0)
                        Eerr = abs(dot(method.E, fout) / Eold - one(Precision))
                        println("Energy error before adaptive: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))
                        #println("Energy in δ: ", dot(method.E, δ), " Energy in ϕ3: ", dot(method.E, ϕ3))

                        ηE = Eerr

                        @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)


                        #atol = 1e-5 * maximum(fold)
                        #rtol = 1e-4
                        #err = fout .- fstar
                        #scale = atol .+ rtol .* max.(abs.(fold),abs.(fout))

                        #η = sqrt(mean((err ./ scale).^2))
                        #println("η for second step: ", norm(η))

                    # Form F with fstar
                        #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fstar,vol * dtldt0,zero(Precision))
                        #J2 = 2 .* method.M_Bin_Mul_Step 
                        #J2 .+= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx]) * dtldt0
                        #threshold_Jcols!(J2;rtol = 1f-7)
                    # Form F
                        #method.M_Bin_Mul_Step .+= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx]) * dtldt0
                        #mul!(F2,method.M_Bin_Mul_Step,fstar)
                        #F2 .+= df_Inj * dtldt0
                        #threshold_F!(F2,fstar,dt_local;rtol=1f-2)
                        # error calculation for first step and adaptivity
                        #ηR = norm(fstar - fold -k*F2)/(norm(δ)+eps(Precision)) # want this to be about <10^-3 for accuracy
                        #ηJ = norm(J2-J)/norm(J)
                        #ηδ = norm(δ)/norm(fold)
                        #println("ηR: ", ηR)
                        #println("ηJ: ", ηJ)
                        #println("ηδ: ", ηδ)

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        k = (1e-4/(ηE+eps(Float64)))^(1.0/(order+1)) 
                        #k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        dt_local = min(1.2*dt_old,k*dt_old)
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space, "t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space, "t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old
                        #k = 1.0

                        #println("k: ", k)

                        #if dt_local < 1e-5
                        #    error("small step")
                        #end


                    # second step with adapted dt
                    #k=0.99

                    #β = norm(F)
                    #e1 = zeros(Precision,size(V,2))
                    #e1[1] = one(Precision)

                    if k != 1.0
                        phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=false)
                        mul!(δ,D,@view(ϕ[:,2]),k*invA,zero(Precision))
                        fout .= fold .+ δ
                        @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    end

                    # form D2 with δ only need A(δ,δ)
                    #=mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,δ,vol*dtldt0,zero(Precision))
                    mul!(D2,method.M_Bin_Mul_Step,δ)

                    D2 .= D \ D2

                    D2[abs.(D2) .< eps(Precision)] .= zero(Precision) # remove small entries in D2 to improve stability of second step

                    # compute second step ϕ functions
                    Ks2 = KrylovSubspace{Float64}(n_momentum, m)
                    arnoldi!(Ks2, J, D2;m=m) 
                    V2 = ExponentialUtilities.getV(Ks2)[:,1:end-1]
                    H2 = ExponentialUtilities.getH(Ks2)[1:end-1,1:end]

                    β = norm(D2)
                    e1 = zeros(Precision,m)
                    e1[1] = one(Precision)
                    ϕ3 = 2 * k * D * V2 * phiv(k,H2,β*e1,3)[:,4]=#
                    
                    # combine steps for final solution
                    
                    #@. fout = fold + δ #+ ϕ3

                    #println(minimum(fout), " ", maximum(fout))

                    #@. fout = max(fout, zero(Precision)) # remove negative values before error calculation, can cause issues with error estimation and stability of second step

                    # energy correction  
                    #ΔE = dot(method.E, fout) - dot(method.E, fold .+ df_Inj * k * dtldt0)
                    #Eold = dot(method.E, fold .+ df_Inj * k * dtldt0)
                    #Eerr = abs(dot(method.E, fout) / Eold - one(Precision))
                    #println("Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))
                    #println("Energy in δ: ", dot(method.E, δ), " Energy in ϕ3: ", dot(method.E, ϕ3))

                    #@. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fout = max(fout, zero(Precision)) # remove negative values after error calculation, can cause issues with stability if negative values are used in next step

                    # energy correction  
                    ΔE = dot(method.E, fout) - dot(method.E, fold .+ df_Inj * k * dtldt0)
                    Eold = dot(method.E, fold .+ df_Inj * k * dtldt0)
                    Eerr = abs(dot(method.E, fout) / Eold - one(Precision))
                    #println("Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                    if Eerr > 1e-4
                        @warn "Energy error in EXPRB32 step is large $Eerr, may be unstable, consider reducing time step or adjusting ηtarget"
                    end

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        else #

            fold .= @view(method.fstep[start_idx:end_idx]) 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            while t < 1.0

                dtldt0 = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    J .= M_Emi
                    #J .-= P_Flux
                    lmul!(dtldt0,J)
                    # Form F
                    mul!(Ftmp,J,fold)
                    Ftmp .+= df_Inj
                    lmul!(A * dtldt0,Ftmp)
                    # remove small changes in J and F 
                    #threshold_Jrows!(J;rtol = eps(Precision))
                    #threshold_Jcols!(J;rtol = eps(Precision))
                    #J[abs.(J) .< 1f-20] .= zero(Precision) # remove small entries in D2 to improve stability of second step

                    println("max J: ", maximum(J), " min J: ", minimum(J))
                    println("max F: ", maximum(Ftmp), " min F: ", minimum(Ftmp))
                    println("max fold: ", maximum(fold), " min fold: ", minimum(fold))

                    #threshold_F!(F,fold,k;rtol=1f-4)

                    # Form Krylov Subspace
                    #Dn , Jn = balance_similarity(J; maxiter=5, tol=1e-3)
                    #println("cond Dn: ", cond(Dn))

                    #D.diag .= 1f0
                    mul!(Jtmp,J,D)
                    mul!(J,Dinv,Jtmp)
                    #@time J .= D \ J * D
                    mul!(F,Dinv,Ftmp)
                    #@time F .= D \ F

                    #threshold_Jrows!(J;rtol = eps(Precision))
                    #threshold_Jcols!(J;rtol = eps(Precision))
                    #threshold_F!(F,D \ fold,k;rtol=1f-2)

                    #vals = abs.(J)

                    #println("minimum abs(nonzeros(J)): ", minimum(vals[vals .> 0]))
                    #println("maximum abs(nonzeros(J)): ", maximum(vals))

                    #Jn = Dn \ J * Dn
                    #Fn = Dn \ F

                    #J64 = Float64.(J)
                    #F64 = Float64.(F)

                    #F64n = Float64.(Fn)
                    #J64n = Float64.(Jn)
                    #println("cond J: ", cond(I -J))
     
                    arnoldi!(Ks64,J,F;m=m,reorthogonalize=true)
 
                    #arnoldi!(Ks32n, Jn, Fn;m=m)
                    #arnoldi!(Ks64n, J64n, F64n;m=m)

                    arnoldi!(Ks32, J, F;m=m,reorthogonalize=true)
                    #arnoldi!(Ks64, J64, F64;m=m)

                    #V32 = ExponentialUtilities.getV(Ks32)[:,1:end-1]
                    #H32 = ExponentialUtilities.getH(Ks32)[1:end-1,1:end]
                    V64 = ExponentialUtilities.getV(Ks64)[:,1:end-1]
                    H64 = ExponentialUtilities.getH(Ks64)[1:end-1,1:end]
                    #V32n = ExponentialUtilities.getV(Ks32n)[:,1:end-1]
                    #H32n = ExponentialUtilities.getH(Ks32n)[1:end-1,1:end]
                    #V64n = ExponentialUtilities.getV(Ks64n)[:,1:end-1]
                    #H64n = ExponentialUtilities.getH(Ks64n)[1:end-1,1:end]

                    #norm_scaled = norm(V32s' * V32s - I)
                    norm_unscaled = norm(V64' * V64 - I)
                    if isnan(norm_unscaled) # no J or F means no Krylov subspace generated, so just accept step as is and move on,
                        break
                    end
                    #norm_new = norm(V32n' * V32n - I)
                    #println("|V'V - I| Float64 scaled: ", norm(V64s' * V64s - I))  # TODO: scale switch between for when norm changes
                    #println("|V'V - I| Float32 scaled: ", norm_scaled)
                    #println("|V'V - I| Float64: ", norm(V64' * V64 - I))
                    println("|V'V - I| Float32: ", norm_unscaled)
                    #println("|V'V - I| Float64 balanced: ", norm(V64n' * V64n - I))
                    #println("|V'V - I| Float32 balanced: ", norm_new)
                    #println("cond((I - dt_scale*H64)) scaled: ", cond((I - dt*H64s)))
                    #println("cond((I - dt_scale*H32)) scaled: ", cond((I - dt*H32s)))
                    #println("cond((I - dt_scale*H64)): ", cond((I - dt*H64)))
                    println("cond((I - dt_scale*H64)): ", cond((I - dt*H64)))
                    #println("cond((I - dt_scale*H64)) balanced: ", cond((I - dt*H64n)))
                    #println("cond((I - dt_scale*H32)) balanced: ", cond((I - dt*H32n)))

                    if cond((I - dt*H64)) > 1f8 ||  norm(V64' * V64 - I) > 1f-5 
                        dt_local = dt_local * 0.5
                        continue
                    end
                    #cond_scaled = cond((I - dt*H64s))
                    #cond_unscaled = cond((I - dt*H64))

                    #if cond_scaled < cond_unscaled
                    #    println("Using scaled Krylov subspace for error estimation")
                    #    Ks = Ks64s
                    #    V = V64s
                    #    H = H64s
                    #    F = Fs
                    #    J = Js
                    #    D = D
                    #else
                        #println("Using unscaled Krylov subspace for error estimation")
                        #V = V64
                        #H = H64
                        #F = F
                        #J = J
                        #D.diag .= 1.0
                    #end
                    #println("size V : ", size(V))
                    #println("size H : ", size(H))

                    #println(abs(ExponentialUtilities.getH(Ks)[m+1,m]))
                    #β = norm(F)
                    #e1 = zeros(Precision,size(V,2))
                    #e1[1] = one(Precision)

                    # Compute φ functions of H
                    phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true) # TODO: This allocates
                    mul!(δ,D,@view(ϕ[:,2]),k*invA,zero(Precision))
                    #@time δ .= invA_Flux * δ # TODO: this allocates
                    fout .= fold .+ δ

                    # adaptive time stepping

                        @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                        # energy error
                        ΔE = dot(method.E, fout) - dot(method.E, fold .+ df_Inj * k*dtldt0)
                        Eold = dot(method.E, fold .+ df_Inj * k * dtldt0)
                        Eerr = abs(dot(method.E, fout) / Eold - one(Precision))
                        println("Energy error before adaptive: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))
                        #println("Energy in δ: ", dot(method.E, δ), " Energy in ϕ3: ", dot(method.E, ϕ3))

                        ηE = Eerr

                        #atol = 1e-5 * maximum(fold)
                        #rtol = 1e-4
                        #err = fout .- fstar
                        #scale = atol .+ rtol .* max.(abs.(fold),abs.(fout))

                        #η = sqrt(mean((err ./ scale).^2))
                        #println("η for second step: ", norm(η))

                    # Form F with fstar
                        #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fstar,vol * dtldt0,zero(Precision))
                        #J2 = 2 .* method.M_Bin_Mul_Step 
                        #J2 .+= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx]) * dtldt0
                        #threshold_Jcols!(J2;rtol = 1f-7)
                    # Form F
                        #method.M_Bin_Mul_Step .+= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx]) * dtldt0
                        #mul!(F2,method.M_Bin_Mul_Step,fstar)
                        #F2 .+= df_Inj * dtldt0
                        #threshold_F!(F2,fstar,dt_local;rtol=1f-2)
                        # error calculation for first step and adaptivity
                        #ηR = norm(fstar - fold -k*F2)/(norm(δ)+eps(Precision)) # want this to be about <10^-3 for accuracy
                        #ηJ = norm(J2-J)/norm(J)
                        #ηδ = norm(δ)/norm(fold)
                        #println("ηR: ", ηR)
                        #println("ηJ: ", ηJ)
                        #println("ηδ: ", ηδ)

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        k = (1e-5/(ηE+eps(Float64)))^(1.0/(order+1)) 
                        k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        #dt_local = min(1.1*dt_old,k*dt_old)
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space, "t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space, "t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                        #println("k: ", k)

                        #if dt_local < 1e-5
                        #    error("small step")
                        #end


                    # second step with adapted dt
                    #k=0.99

                    #β = norm(F)
                    #e1 = zeros(Precision,size(V,2))
                    #e1[1] = one(Precision)

                    if k != 1.0
                        phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true)
                        mul!(δ,D,@view(ϕ[:,2]),k*invA,zero(Precision))
                        fout .= fold .+ δ
                        @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    end

                    # energy correction  
                    ΔE = dot(method.E, fout) - dot(method.E, fold .+ df_Inj * k * dtldt0)
                    Eold = dot(method.E, fold .+ df_Inj * k * dtldt0)
                    Eerr = abs(dot(method.E, fout) / Eold - one(Precision))
                    #println("Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                    if Eerr > 1e-4
                        @warn "Energy error in EXPRB32 step is large $Eerr, may be unstable, consider reducing time step or adjusting ηtarget"
                    end

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        end

    end

    return nothing

end

@inline function arnoldi_orthogonalize!(
    y,
    V,
    H,
    j::Integer,
    lo::Integer,
    hi::Integer,
    U,
)
    @inbounds for i in lo:hi
        α = ExponentialUtilities.coeff(U, dot(@view(V[:, i]), y))
        H[i, j] += α
        axpy!(-α, @view(V[:, i]), y)
    end
    return nothing
end

function arnoldi_step!(
    j::Integer, iop::Integer, A::AT,
    V::AbstractMatrix{T}, H::AbstractMatrix{U},
    n::Int = -1, p::Int = -1;
    reorthogonalize::Bool = false,
) where {AT, T, U}

    x, y = @view(V[:, j]), @view(V[:, j + 1])
    ExponentialUtilities.applyA!(y, A, x, V, j, n, p)

    lo = max(1, j - iop + 1)

    # First orthogonalization pass
    arnoldi_orthogonalize!(y, V, H, j, lo, j, U)

    # Second pass (full reorthogonalization over the same window)
    if reorthogonalize
        arnoldi_orthogonalize!(y, V, H, j, lo, j, U)
    end

    β = H[j + 1, j] = norm(y)
    if !iszero(β)
        @. y /= β
    end
    return β
end

function arnoldi!(
    Ks::KrylovSubspace{T1, U}, A::AT, b;
    tol::Real = 1.0e-7,
    m::Int = min(Ks.maxiter, size(A, 1)),
    ishermitian::Bool = LinearAlgebra.ishermitian(A isa Tuple ? first(A) : A),
    opnorm = nothing,
    iop::Int = 0,
    reorthogonalize::Bool = false,
    init::Int = 0,
    t::Number = NaN,
    mu::Number = NaN,
    l::Int = -1,
) where {T1 <: Number, U <: Number, AT}

    Ks.wasbreakdown = false

    ishermitian && return lanczos!(Ks, A, b; tol = tol, m = m, init = init, t = t, mu = mu, l = l)

    m > Ks.maxiter ? resize!(Ks, m) : Ks.m = m

    @inbounds V, H = ExponentialUtilities.getV(Ks), ExponentialUtilities.getH(Ks)

    b′, b_aug, n, p = ExponentialUtilities.checkdims(A, b, V)
    if iszero(init)
        isaugmented = AT <: Tuple
        isaugmented ? ExponentialUtilities.firststep!(Ks::KrylovSubspace, V, H, b′, b_aug, t, mu, l) :
                      ExponentialUtilities.firststep!(Ks::KrylovSubspace, V, H, b)
        init = 1
    end

    iszero(Ks.beta) && return Ks
    iszero(iop) && (iop = m)

    for j in init:m
        beta = arnoldi_step!(j, iop, A, V, H, n, p; reorthogonalize = reorthogonalize)
        if beta < tol
            Ks.m = j
            Ks.wasbreakdown = true
            break
        end
    end

    return Ks
end

"""
    threshold_Jrows!(J; rtol=1e-14)

Removes weak coupling in the rows of the jacobian matrix J, i.e. elements of the equation for each state variable that are small compared to the largest element of that equation.
"""
function threshold_Jrows!(J; rtol=1e-14)
    m, n = size(J)

    for i in 1:m
        rowmax = maximum(abs, @view J[i,:])

        if rowmax == 0
            continue
        end

        thresh = rtol * rowmax

        for j in 1:n
            if i!=j && sign(J[i,j]) == 1 && abs(J[i,j]) < thresh # doesn't remove diagonal entries 
                J[i,j] = zero(eltype(J))
            end
        end
    end

    return nothing
end


"""
    threshold_Jcols!(J; rtol=1e-14)

Removes weak coupling in the columns of the jacobian matrix J, i.e. elements of the equation for each state variable that are small compared to the largest element of that equation.
"""
function threshold_Jcols!(J; rtol=1e-14)

    m, n = size(J)

    for j in 1:n

        colmax = maximum(abs, @view J[:,j])

        if colmax == 0
            continue
        end

        thresh = rtol * colmax

        for i in 1:m
            if i != j && sign(J[i,j]) == 1 && abs(J[i,j]) < thresh # postive and small
                J[i,j] = zero(eltype(J))
            end
        end
    end

    return J
end

"""
    threshold_F!(F; rtol=1e-14)

Removes weak transition rates in F, i.e. terms df_i/dt = F_i that are small and by small I mean f^{t+1}=f^{t}+ΔtF so if ΔtF < eps(f^{t})*rtol.
"""
function threshold_F!(F,f,dt; rtol=1e-2)
    
    @. F = ifelse(dt*abs(F) < f*rtol, zero(eltype(F)), F)
    #@. F = ifelse(abs(F) < rtol, zero(eltype(F)), F)

    return nothing
end

"""
    balance_similarity(J; maxiter=50, tol=1e-3)

Compute a diagonal scaling D such that

    Js = D \\ J * D

has approximately balanced row and column 2-norms.

Returns (D, Js).
"""
function balance_similarity(J; maxiter=50, tol=1e-3)

    n = size(J,1)
    @assert size(J,2) == n

    d = ones(eltype(J), n)

    Js = copy(J)

    for k = 1:maxiter

        maxerr = 0.0

        for i = 1:n

            r = norm(view(Js, i, :))
            c = norm(view(Js, :, i))

            if r == 0 || c == 0
                continue
            end

            ratio = c / r

            ratio = clamp(ratio, 1e-16, 1e16) # prevent extreme scaling

            # desired scaling
            s = sqrt(ratio)

            d[i] *= s
            d[i] = clamp(d[i], 1e-8, 1e8) # prevent extreme scaling

            # similarity update
            Js[:,i] .*= s
            Js[i,:] ./= s

            maxerr = max(maxerr, abs(log(ratio)))
        end

        maxerr < tol && break
    end

    return Diagonal(d), Js
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

function update_Big_Bin!(method::BackwardEulerStruct,dt_scale::T) where T
    
    n_momentum=size(method.M_Bin,2)
    
    f_iter = zeros(T,n_momentum)
    f_guess = zeros(T,n_momentum)
    f_tmp = zeros(T,n_momentum)
    df_tmp = zeros(T,n_momentum)
    J = method.J
    F = method.F

    tol = T(1e-6)
    maxiter = 10
    damping = false

    for off_space in method.Bin_Domain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        f = @view(method.f[start_idx:end_idx])
        fstep = @view(method.fstep[start_idx:end_idx])
        @inbounds vol = method.Vol[off_space+1] * dt_scale
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        
        f_guess .= f

        for λ in [T(1e-12), T(1e-10), T(1e-8), T(1e-6), T(1e-4), T(1e-2), T(1.0)]

            f_iter .= f_guess

            # Newton iteration
            for it in 1:maxiter

                println("Newton iteration: $it, λ: $λ")

                build_residual_jacobian!(method, f_iter, f, invA_Flux, λ * vol)

                if norm(F, Inf) < tol
                    f_guess .= f_iter
                    break
                end

                println("Residual norm: $(norm(F, Inf))")

                println("cond J: $(cond(J))")

                δ = J \ F

                δ .= δ.* (abs.(δ) .> 0.0) # remove small values
                
                display(reshape(δ, 37,24))
                display(reshape(δ ./ f_iter, 37,24))

                f_iter .-= δ

                f_iter .= f_iter .* (f_iter .> 0.0) # enforce positivity

                # Simple backtracking line search on residual norm
                #=α = one(T)
                Fn = norm(F, 2)

                while α > T(1e-12)
                    trial = f_guess .+ α .* δ
                    build_residual_jacobian!(method, trial, f, invA_Flux, vol)

                    if norm(F, 2) <= (one(T) - T(1e-4) * α) * Fn
                        f_iter .= trial
                        break
                    end

                    α *= T(0.5)
                end

                if α <= T(1e-12)
                    error("Newton line search failed to reduce the residual.")
                end=#

            end

            #error("Newton did not converge in $maxiter iterations.")

            f_guess .= f_iter

        end

        fstep .= f_iter

    end

    return nothing

end

function build_residual_jacobian!(method::BackwardEulerStruct, f_iter, f, invA_Flux, vol)

    J = method.J
    F = method.F

    I_matrix = Diagonal(ones(method.Precision, length(F)))

    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,f_iter,vol,zero(eltype(f_iter))) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

    mul!(J,invA_Flux,method.M_Bin_Mul_Step,-2.0,0.0)

    @. J += I_matrix

    mul!(F,method.M_Bin_Mul_Step,f_iter)

    tmp = f_iter .- f .- (invA_Flux * F) 

    F .= tmp

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
    fstep = zeros(Float64, n_momentum)
    df_Inj = zeros(Float64, n_momentum)
    ftmp = zeros(Float64, n_momentum)
    fs = zeros(Float64,size(fstep))
    sigma = Diagonal(zeros(Float64, n_momentum))
    fdivsigma = Diagonal(zeros(Float64, n_momentum))
    ftildedivsigma = Diagonal(zeros(Float64, n_momentum))
    tmpdiag = Diagonal(zeros(Float64, n_momentum))
    onesvec = ones(Float64, n_momentum)
    source = Diagonal(zeros(Float64, n_momentum))

    for off_space in 0:n_space-1 

        fill!(method.Q,zero(method.Precision)) # reset Q
        fill!(tmpdiag.diag,zero(eltype(tmpdiag.diag))) # reset tmpdiag
        fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        fill!(Qflux,zero(eltype(Qflux))) # reset Qflux

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)
        println("off space: $off_space, start idx: $start_idx, end idx: $end_idx")

        P_Flux_View = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        invA_Flux_View = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        if method.Emission_Interactions
            M_Emi_View = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        end
        df_Inj .= @view(method.df_Inj[start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx]) 
        @. fdivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        #fstep .+= (df_Inj + invA_Flux_View * M_Emi_View * fstep) * dt_scale
        #fstep[fstep .< 1e16eps(0.0)] .= 0.0 # remove negative values for stability, also sets small values to zero which allows for more sparsity and faster solves in the Patankar update
        @. sigma.diag = max(fstep,eps(Float64))
        @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0) #fstep ./ sigma.diag # this is just fstep but with a cutoff to avoid division by zero, also allows for more sparsity and faster solves in the Patankar update
        fill!(Qbin,zero(eltype(Qbin))) # reset Qbin NEEDED otherwise binary interactions from previous off_space will carry over and contaminate current off_space update
        @. sigma.diag = ifelse(fstep > 0.0, fstep,eps(Float64))
        Qflux .= P_Flux_View

        # first do advection update as binary interaction can be be large in comparision.

        #=    mul!(tmpdiag.diag,transpose(Qflux),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
            Qtmp .= invA_Flux_View * (Qflux * ftildedivsigma - tmpdiag) * dt_scale

            #Qtmp .-= invA_Flux_View * M_Emi_View * ftildedivsigma * dt_scale # minus as M_Emi is on LHS of transport equation and Q is on RHS

            #mul!(source.diag, M_Emi_View, fstep)
            #@. source.diag /= sigma.diag
            #Qbin .-= invA_Flux_View * source * dt_scale # minus as M_Bin is on LHS of transport equation and Q is on RHS

            #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            Qtmp .+= I_matrix # add identity for Patankar update 

            ftmp .= fstep

            QLU = lu(Qtmp)
            ldiv!(QLU,ftmp)

            fstep .= ftmp
            =#
            # TODO: if dt too big this leads to poor energy conservation due to excess gain of photons and fast cooling of electrons 

        # now do binary interaction update
        
        if method.Binary_Interactions && off_space in method.Bin_Domain

            @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        
            @inbounds vol = method.Vol[off_space+1]

            #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            # loss
            #mul!(ftmp,method.Lij,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape
            #tmpdiag.diag .= ftmp .* ftildedivsigma.diag
            #Qbin .-= tmpdiag # minus as M_Bin is on LHS of transport equation and Q is on RHS
            mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstep,vol,zero(method.Precision)) 
            mul!(tmpdiag.diag,transpose(method.M_Bin_Mul_Step),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag

            Qbin .= tmpdiag # minus as M_Bin is on LHS of transport equation and Q is on RHS

            # gain 
            #mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            Qbin .-= method.M_Bin_Mul_Step * ftildedivsigma # minus as M_Bin is on LHS of transport equation and Q is on RHS


            cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            colsum = sum(Qbin,dims=1)
            rowsum = sum(Qbin,dims=2)
            println("Qbin - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))

            Qtmp .= invA_Flux_View * (Qbin) * dt_scale

            Qtmp .+= I_matrix # add identity for Patankar update

            norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
            Qs .= Dl * Qtmp * Dr

            ftmp .= fstep

            fs .= Dl * ftmp

            QLU = lu(Qs)
            ldiv!(QLU,fs)

            s = svdvals(Qs)
            println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

            ftmp .= Dr * fs

            #fstep .= patankar_energy_fix!(method,ftmp,fstep)

            #fstep .= patankar_energy_fix_quadratic!(method,ftmp,fstep)

            #fstep .= patankar_energy_fix_linear!(method,ftmp,fstep)

            fstep .= patankar_energy_fix_linear_filtered!(method,ftmp,fstep)

            #fstep .= ftmp

            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")
            
        end

       #= if method.Emission_Interactions
            mul!(source.diag, M_Emi_View, fstep)
            @. source.diag /= sigma.diag # divide by fstep for source term to be in correct units for Patankar update, also sets source to zero where fstep is zero which allows for more sparsity and faster solves in the Patankar update
            println("max source diag: ", findmax(source.diag), " min source diag: ", findmin(source.diag), " sum source diag: ", sum(source.diag))
        else
            fill!(source.diag, zero(eltype(source.diag)))
        end=#
 



        #=println("dt scale: $dt_scale")
        Qtmp .+= invA_Flux_View * (Qbin) * dt_scale
        colsum = sum(Qtmp,dims=1)
        rowsum = sum(Qtmp,dims=2)
        println("Qtmp - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
        display(Qtmp[115:125,120]')

        colsum = sum(Qtmp,dims=1)
        rowsum = sum(Qtmp,dims=2)
        println("Qtmp - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
        display(Qtmp[115:125,120]')

        cutoff_offdiagonal_to_diag_relative!(Qtmp, 1e-6)
        display(Qtmp[115:125,120]')

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

        #=if kwo > sqrt(1/eps(eltype(Qtmp)))
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

        else=#

            if method.Emission_Interactions
                ftmp .= fstep #+ df_Inj .* dt_scale + invA_Flux_View * (M_Emi_View * fstep) .* dt_scale # add M_Emi contribution to ftmp for Patankar update
                #println("min ftmp: ", findmin(ftmp), " max ftmp: ", findmax(ftmp), " sum ftmp: ", sum(ftmp))
                #ftmp[ftmp .< 1e16eps(0.0)] .= 0.0
            else
                ftmp .= fstep #+ df_Inj .* dt_scale
            end

            #norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
            #Qs .= Dl * Qtmp * Dr

            #k = cond(Qs)

            #println("off space: $off_space, cond without: $kwo, condition with right and left: $k")


             
            fs .= Dl * ftmp

            QLU = lu(Qs)
            ldiv!(QLU,fs)
            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")

            #QLU = lu(Qtmp)
            #ldiv!(QLU,ftmp)

            fstep .= Dr * fs
            #fstep .= ftmp

            #gmres!(method.GMRESWorkspace,Qtmp,ftmp;rtol = 1f-10,atol = 0f0)
            #fstep .= method.GMRESWorkspace.x

        #end
        =#

        @view(method.f_step[start_idx:end_idx]) .= fstep

    end

    return nothing

end

function updateMomentumPatankarEulerSymmetric!(method::SymplecticSymmetricMPEStruct,f::AbstractVector{T},dt_scale::T) where T

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
    fstep = zeros(Float64, n_momentum)
    df_Inj = zeros(Float64, n_momentum)
    ftmp = zeros(Float64, n_momentum)
    fs = zeros(Float64,size(fstep))
    sigma = Diagonal(zeros(Float64, n_momentum))
    fdivsigma = Diagonal(zeros(Float64, n_momentum))
    ftildedivsigma = Diagonal(zeros(Float64, n_momentum))
    tmpdiag = Diagonal(zeros(Float64, n_momentum))
    onesvec = ones(Float64, n_momentum)
    source = Diagonal(zeros(Float64, n_momentum))

    for off_space in 0:n_space-1 

        fill!(method.Q,zero(method.Precision)) # reset Q
        fill!(tmpdiag.diag,zero(eltype(tmpdiag.diag))) # reset tmpdiag
        fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        fill!(Qflux,zero(eltype(Qflux))) # reset Qflux

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)
        println("off space: $off_space, start idx: $start_idx, end idx: $end_idx")

        P_Flux_View = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        invA_Flux_View = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        if method.Emission_Interactions
            M_Emi_View = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        end
        df_Inj .= @view(method.df_Inj[start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx]) 
        @. fdivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        #fstep .+= (df_Inj + invA_Flux_View * M_Emi_View * fstep) * dt_scale
        #fstep[fstep .< 1e16eps(0.0)] .= 0.0 # remove negative values for stability, also sets small values to zero which allows for more sparsity and faster solves in the Patankar update
        @. sigma.diag = max(fstep,eps(Float64))
        @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0) #fstep ./ sigma.diag # this is just fstep but with a cutoff to avoid division by zero, also allows for more sparsity and faster solves in the Patankar update
        fill!(Qbin,zero(eltype(Qbin))) # reset Qbin NEEDED otherwise binary interactions from previous off_space will carry over and contaminate current off_space update
        @. sigma.diag = ifelse(fstep > 0.0, fstep,eps(Float64))
        Qflux .= P_Flux_View

        # first do advection update as binary interaction can be be large in comparision.

            #=mul!(tmpdiag.diag,transpose(Qflux),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
            Qtmp .= invA_Flux_View * (Qflux * ftildedivsigma - tmpdiag) * dt_scale

            #Qtmp .-= invA_Flux_View * M_Emi_View * ftildedivsigma * dt_scale # minus as M_Emi is on LHS of transport equation and Q is on RHS

            #mul!(source.diag, M_Emi_View, fstep)
            #@. source.diag /= sigma.diag
            #Qbin .-= invA_Flux_View * source * dt_scale # minus as M_Bin is on LHS of transport equation and Q is on RHS

            #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            Qtmp .+= I_matrix # add identity for Patankar update 

            ftmp .= fstep

            QLU = lu(Qtmp)
            ldiv!(QLU,ftmp)

            fstep .= ftmp=#
            
            # TODO: if dt too big this leads to poor energy conservation due to excess gain of photons and fast cooling of electrons 

        # now do binary interaction update
        
        if method.Binary_Interactions && off_space in method.Bin_Domain

            fstepE = fstep #.* method.E

            @. ftildedivsigma.diag = ifelse(fstepE > 0.0, 1.0,0.0)
            @. sigma.diag = ifelse(fstepE > 0.0, 1/fstepE,1/eps(Float64))
        
            @inbounds vol = method.Vol[off_space+1]

            #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            # loss ij
                fill!(Qbin,zero(eltype(Qbin))) # reset Qbin for loss calculation

                mul!(method.M_Bin_Mul_Step_reshape,method.Lijk,fstepE,vol,zero(method.Precision))
                mul!(tmpdiag.diag,method.M_Bin_Mul_Step,fstepE) 

                println(isdiag(method.M_Bin_Mul_Step))

                Qbin .-= tmpdiag * sigma #./2 #* sigma #=+ method.M_Bin_Mul_Step * ftildedivsigma)/2=# # minus as M_Bin is on LHS of transport equation and Q is on RHS

                #mul!(method.M_Bin_Mul_Step_reshape,method.Lijk,fstepE,vol,zero(method.Precision))
                #mul!(tmpdiag.diag,method.M_Bin_Mul_Step,fstepE) 

                #Qbin .-= sigma * tmpdiag ./2 #* sigma #=+ method.M_Bin_Mul_Step * ftildedivsigma)/2=# # minus as M_Bin is on LHS of transport equation and Q is on RHS

                # gain 
                mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstepE,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                Qbin .-= (method.M_Bin_Mul_Step * ftildedivsigma) #./2 #+ ftildedivsigma * method.M_Bin_Mul_Step) / 2 # minus as M_Bin is on LHS of transport equation and Q is on RHS

                #mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,ftildedivsigma.diag,vol,zero(method.Precision))
                #Qbin .-= (method.M_Bin_Mul_Step * diagm(fstepE)) ./2
                # together
                #mul!(method.M_Bin_Mul_Step_reshape,method.Aijk,fstepE,vol,zero(method.Precision))
                #Qbin .= method.M_Bin_Mul_Step
                #Qbin .= -Qbin * ftildedivsigma 

                #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e0)

                #mul!(tmpdiag.diag,transpose(Qflux),onesvec)
                #tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
                #Qbin .+=  (Qflux * ftildedivsigma - tmpdiag)

                display(sparse(Qbin))

                println("has_positive_offdiag(Qbin): $(has_positive_offdiag(Qbin))")

                Qtmp .= invA_Flux_View * (Qbin) * dt_scale

                colsum = sum(Qtmp,dims=1)
                rowsum = sum(Qtmp,dims=2)
                println("Qbin - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
                println("Qbin - max: ", findmax(Qtmp[:,findmax(colsum)[2][2]]), " min: ", findmin(Qtmp[:,findmax(colsum)[2][2]]))
                #display(reshape(Qtmp[:,findmax(colsum)[2][2]],37,24))
                display(reshape(Qtmp[:,150],37,24))
                println("sum: $(sum(Qtmp[:,150])), Esum: $(sum(method.E .* Qtmp[:,150]))")

                println("sumQf ", sum(Qtmp * fstep), " sumEQf ", sum(method.E .*(Qtmp * fstep)), " maxQf ", maximum(Qtmp * fstep), " minQf ", minimum(Qtmp * fstep))
                println("sumQf ele ", sum(Qtmp[1:288,:] * fstep), " sumEQf ele ", sum(method.E[1:288] .*(Qtmp[1:288,:] * fstep)), "sumQf pho ", sum(Qtmp[289:end,:] * fstep), "  sumEQf pho ", sum(method.E[289:end] .*(Qtmp[289:end,:] * fstep)))

                println("max timescale: ", maximum(replace!(abs.(max.(fstep,eps(Float64)) ./ (Qtmp * fstep)),NaN=>0.0,0.0=>Inf)), " min timescale: ", minimum(replace(abs.(fstep ./ (Qtmp * fstep)),NaN=>Inf,0.0=>Inf)))

                Qtmp .+= I_matrix # add identity for Patankar update

                norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
                Qs .= Dl * Qtmp * Dr

                println("max 1TQ ", maximum(transpose(onesvec) * Qtmp - transpose(onesvec)), " min 1TQ ", minimum(transpose(onesvec) * Qtmp - transpose(onesvec)))
                println("max ETQ ", maximum(transpose(method.E) * Qtmp - transpose(method.E))," min ETQ ", minimum(transpose(method.E) * Qtmp - transpose(method.E)))

                println("has_positive_offdiag(Qs): $(has_positive_offdiag(Qs))")

                ftmp .= fstepE

                fs .= Dl * ftmp

                QLU = lu(Qs)
                ldiv!(QLU,fs)

                s = svdvals(Qs)
                println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")
                println("cond = $(cond(Qs))")

                ftmp .= (Dr * fs) #./ method.E

            println("num before $(sum(fstep)), num after $(sum(ftmp)), Δ = $(sum(ftmp)-sum(fstep))")
            println("eng before $(sum(method.E .* fstep)), eng after $(sum(method.E .* ftmp)), Δ = $(sum(method.E .* ftmp) - sum(method.E .* fstep))")

            fstep .= ftmp

            println("any negative fstep: ", any(fstep .< 0.0))
            println("sum of negative fstep: ", sum((fstep .< 0.0) .* fstep))
            println("sum of positive fstep: ", sum((fstep .> 0.0) .* fstep))

            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")
            
        end


        @view(method.f_step[start_idx:end_idx]) .= fstep

    end

    return nothing

end

has_positive_offdiag(M) = any((i != j && M[i, j] > 0) for i in axes(M, 1), j in axes(M, 2))

function patankar_energy_fix!(method::SymplecticMPEStruct,fstep::AbstractVector{T},fprev::AbstractVector{T}) where T

    E = method.E
    df = similar(fstep)
    fill!(df,zero(eltype(df)))
    ftmp = copy(fstep)

    w = zeros(eltype(df), length(df))
    df_active = zeros(eltype(df), length(df))
    df_inactive = zeros(eltype(df), length(df))

    active = trues(length(fstep))
    nonpositive = true

    while nonpositive

        fill!(w, zero(eltype(w)))
        fill!(df_active, zero(eltype(df_active)))
        fill!(df_inactive, zero(eltype(df_inactive)))

    # fixed inactive correction
    df_inactive .= -ftmp .* (.!active)

    # remaining defects caused by clipping
    ΔM = sum(ftmp[.!active])

    ΔE = sum(E .* fprev) - sum(E .* (ftmp .+ df_inactive))

    # weights on active set
    eta = E
    w[active] .= ftmp[active] #.* eta

    S0 = sum(w)
    Ebar = sum(E .* w) / S0

    Ẽ = E .- Ebar

    S2 = sum((Ẽ.^2) .* w)

    λ = -2ΔM / S0
    μ = -2*(ΔE - Ebar*ΔM) / S2
    # active correction ONLY
    df_active[active] .= -0.5 .* w[active] .* (λ .+ μ .* Ẽ[active])

    # total correction
    df = df_inactive + df_active

    println(sum(df_active))
    println("mass defect = ", ΔM)
    println("energy defect = ", ΔE)
    println("constraint mass = ", sum(df_active))
    println("constraint energy = ", sum(E .* df_active))
    println("minimum corrected state = ", minimum(ftmp .+ df))

    ucorr = ftmp .+ df

    new_active = ucorr .> 1e-10

    active .&= new_active

    nonpositive = any(ucorr .< 0.0)

    println("sum(df_inactive) = ", sum(df_inactive))
    println("sum(df_active)   = ", sum(df_active))
    println("total correction = ", sum(df))

    end

    println("mass defect = ", sum(ftmp .+ df) - sum(fprev))
    println("energy defect = ", sum(E .* (ftmp .+ df)) - sum(E .* fprev))


    return ftmp .+ df

end

function patankar_energy_fix_linear!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,) where T

    E = method.E

    #
    # Uncorrected Patankar state
    #
    ftmp = copy(fstep)

    N = length(ftmp)

    #
    # Corrections
    #
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    #
    # Active-set mask
    #
    active = ftmp .> tol

    #
    # Weights
    #
    w = zeros(T, N)

    for iter in 1:maxiter

        old_active = copy(active)

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))

        #
        # Inactive states are clipped to zero
        #
        df_inactive .= -ftmp .* (.!active)

        #
        # Remaining conservation defects
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Active weights
        #
        w[active] .= ftmp[active]

        #
        # Weighted energy centroid
        #
        W0 = sum(w)

        if W0 <= tol
            error("No active states remain in energy correction.")
        end

        Ebar = sum(E .* w) / W0

        #
        # Centered energies for conditioning
        #
        Ẽ = E .- Ebar

        #
        # Weighted moments
        #
        S0 = sum(w)
        S2 = sum((Ẽ.^2) .* w)

        #
        # Solve:
        #
        # δu_i = w_i (a + b Ẽ_i)
        #
        # Constraints:
        #
        # Σ δu_i      = ΔM
        # Σ E_i δu_i = ΔE
        #
        # Because Ẽ is centered:
        #
        # Σ w_i Ẽ_i = 0
        #
        # giving a diagonalized system.
        #

        a = ΔM / S0

        #
        # Energy defect relative to centroid
        #
        ΔẼ = ΔE - Ebar * ΔM

        b = ΔẼ / S2

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ)

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        #
        # Corrected state
        #
        ucorr = ftmp .+ df

        #
        # Update active set
        #
        active .&= (ucorr .> 0.0)

        #
        # Convergence test
        #
        if all(active .== old_active)
            break
        end

        if iter == maxiter
            @warn "Linear active-set energy projection reached maxiter."
        end
    end

    fcorrected = ftmp .+ df

    #
    # Diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
end

function patankar_energy_fix_linear_filtered!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,
    energy_window_orders = 4,) where T

    E = method.E

    #
    # Uncorrected Patankar state
    #
    ftmp = copy(fstep)

    N = length(ftmp)

    #
    # Correction vectors
    #
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    #
    # Working arrays
    #
    w      = zeros(T, N)
    filter = zeros(T, N)

    #
    # Permanent inactive mask
    #
    inactive = falses(N)

    #
    # Initial positivity check
    #
    inactive .|= (ftmp .<= tol)

    for iter in 1:maxiter

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))
        fill!(filter,      zero(T))

        #
        # Active set is complement of inactive set
        #
        active = .!inactive

        #
        # Permanently clipped states
        #
        df_inactive .= -ftmp .* inactive

        #
        # Remaining defects AFTER clipping
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Determine occupied energy range
        #
        occupied =
            (ftmp .+ df_inactive) .> tol

        if !any(occupied)
            error("No occupied states remain.")
        end

        Emax_occ = maximum(E[occupied])

        #
        # Smooth logarithmic filter
        #
        for i in eachindex(E)

            if E[i] <= 0 || Emax_occ <= 0
                filter[i] = zero(T)
            else

                x =
                    log10(E[i] / Emax_occ) +
                    energy_window_orders

                filter[i] =
                    clamp(x, zero(T), one(T))
            end
        end

        #
        # Active correction weights
        #
        w[active] .=
            ftmp[active] .* filter[active]

        #
        # Prevent singular systems
        #
        S0 = sum(w)

        if S0 <= tol
            error("Correction weights vanished.")
        end

        #
        # Center energies
        #
        Ebar = sum(E .* w) / S0

        Ẽ = E .- Ebar

        #
        # Weighted variance
        #
        S2 = sum((Ẽ.^2) .* w)

        if S2 <= tol
            error("Energy variance too small.")
        end

        #
        # Linear correction:
        #
        # δu_i = w_i (a + b Ẽ_i)
        #
        a = ΔM / S0

        ΔẼ = ΔE - Ebar * ΔM

        b = ΔẼ / S2

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ)

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        #
        # Corrected state
        #
        ucorr = ftmp .+ df

        #
        # Detect NEW negative states
        #
        newly_inactive =
            (ucorr .<= tol) .& (.!inactive)

        #
        # Diagnostics
        #
        println("iter = ", iter)

        println("ΔM target         = ", ΔM)
        println("ΔE target         = ", ΔE)

        println("constraint mass   = ", sum(df))
        println("constraint energy = ", sum(E .* df))

        println("minimum corrected = ", minimum(ucorr))

        #
        # Converged?
        #
        if !any(newly_inactive)
            break
        end

        #
        # Permanently accumulate clipped states
        #
        inactive .|= newly_inactive

        if iter == maxiter
            @warn "Energy projection reached maxiter."
        end
    end

    #
    # Final corrected state
    #
    fcorrected = ftmp .+ df

    #
    # Final diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
end

function patankar_energy_fix_quadratic!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,) where T

    E = method.E

    # unconstrained Patankar result
    ftmp = copy(fstep)

    N = length(fstep)

    # correction vectors
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    # active set
    active = ftmp .> tol

    # weights
    w = zeros(T, N)

    for iter in 1:maxiter

        old_active = copy(active)

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))

        #
        # Inactive states are clipped to zero
        #
        df_inactive .= -ftmp .* (.!active)

        #
        # Remaining conservation defects
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Active-set weights
        #
        w[active] .= ftmp[active]

        #
        # Center energies for conditioning
        #
        W0 = sum(w)

        if W0 <= tol
            error("No active states remain in correction.")
        end

        Ebar = sum(E .* w) / W0

        Ẽ = E .- Ebar

        #
        # Build weighted moments
        #
        M0 = sum(w)
        M1 = sum(Ẽ .* w)
        M2 = sum((Ẽ.^2) .* w)
        M3 = sum((Ẽ.^3) .* w)
        M4 = sum((Ẽ.^4) .* w)

        #
        # We seek:
        #
        # δu_i = w_i (a + b Ẽ_i + c Ẽ_i^2)
        #
        # subject to:
        #
        # Σ δu_i = ΔM
        # Σ Ẽ_i δu_i = ΔẼ
        #
        # and minimizing weighted norm.
        #
        # This gives the normal equations:
        #
        # [M0 M1 M2] [a]   [ΔM ]
        # [M1 M2 M3] [b] = [ΔẼ]
        # [M2 M3 M4] [c]   [0   ]
        #
        # Last row enforces minimum-norm quadratic correction.
        #

        ΔẼ = ΔE - Ebar * ΔM

        A = [
            M0 M1 M2
            M1 M2 M3
            M2 M3 M4
        ]

        rhs = T[ΔM, ΔẼ, 0]

        coeffs = A \ rhs

        a, b, c = coeffs

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ .+ c .* (Ẽ.^2))

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        ucorr = ftmp .+ df

        #
        # Update active set
        #
        active .&= (ucorr .> 0.0)

        #
        # Converged?
        #
        if all(active .== old_active)
            break
        end

        if iter == maxiter
            @warn "Active-set projection reached maxiter."
        end
    end

    #
    # Final corrected state
    #
    fcorrected = ftmp .+ df

    #
    # Diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
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

function cutoff_offdiagonal_to_diag_relative!(A, scale)
    n, m = size(A)
    @assert n == m "A must be square"

    for j in 1:n
        colmax = maximum(A[i, j] for i in 1:n)
        cutoff = scale *eps(colmax)

        for i in 1:n
            if i != j && abs(A[i, j]) < cutoff
                removed = A[i, j]
                A[i, j] = zero(eltype(A))
                A[j, j] -= removed  # Add the value to diagonal
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
