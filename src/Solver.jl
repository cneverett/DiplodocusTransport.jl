function Solve(f1D0::fType,method::SteppingMethodType;save_steps::Int=1,progress=false)

    PhaseSpace = method.PhaseSpace
    Time = PhaseSpace.Time

    f0 = copy(f1D0)
    t_low = Time.t_low
    t_up = Time.t_up
    t_num = Time.t_num

    t_steps = range(t_low,t_up,length=t_num+1)
    n_save = ceil(Int64,t_num/save_steps+1)

    tmp = copy(f0)
    dtmp = similar(f0)

    t = t_low

    output = OutputStruct(f0,n_save)

    # save initial state (step 1)
    output.f[1] = copy(tmp)
    output.t[1] = t
    save_count = 1

    # progress bar
    if progress
        p = Progress(t_num)
    end

    for i in 1:t_num

        t = t_steps[i]
        dt = t_steps[i+1] - t_steps[i]

        method(dtmp,tmp,t,dt)
        @. tmp += dtmp

        # removing values less than 0f0
        @. tmp = tmp*(tmp>=0f0)
        # hackey fix for inf values
        @. tmp = tmp*(tmp!=Inf)

        # saving state
        if (i)%save_steps == 0
            save_count += 1
            t = t_steps[i+1]
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

