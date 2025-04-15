function Solve(f1D0::fType,method::SteppingMethodType;save_steps::Int=1,progress=false)

    PhaseSpace = method.PhaseSpace
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids

    f0 = copy(f1D0)
    tr = Grids.tr
    dt0 = tr[2] - tr[1]
    t_num = length(tr)-1

    n_save = ceil(Int64,t_num/save_steps+1)

    tmp = copy(f0)
    dtmp = similar(f0)

    output = OutputStruct(f0,n_save)

    # save initial state (step 1)
    output.f[1] = copy(tmp)
    output.t[1] = tr[1]
    save_count = 1

    # progress bar
    if progress
        p = Progress(t_num)
    end

    for i in 1:t_num

        dt = tr[i+1] - tr[i]

        method(dtmp,tmp,dt0,dt)
        @. tmp += dtmp

        # removing values less than 0f0
        @. tmp = tmp*(tmp>=0f0)
        # hackey fix for inf values
        @. tmp = tmp*(tmp!=Inf)

        # saving state
        if (i)%save_steps == 0
            save_count += 1
            t = tr[i+1]
            output.f[save_count] = copy(tmp)
            output.t[save_count] = t
        end
        
        if progress
            next!(p)
        end
    end
    if progress
        finish!(p)
    end

    return output

end

