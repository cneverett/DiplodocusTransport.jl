function Solve(f1D0::Vector{F},method::SteppingMethodType;save_steps::Int=1,progress=false,fileName::String=nothing,fileLocation::String=cwd()) where F<:AbstractFloat

    PhaseSpace = method.PhaseSpace
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids

    f0 = copy(f1D0)
    tr = Grids.tr
    dt0 = tr[2] - tr[1]
    t_num = length(tr)-1

    save_idx = unique([1 ; range(save_steps,t_num,step=save_steps)]) # saves first time step and then on an interval of save_steps (no duplicates)

    n_save = length(save_idx)+1 # +1 for the initial state

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

    #filter_vector = zeros(Bool, length(tmp))

    for i in 1:t_num

        dt = tr[i+1] - tr[i]
        t = tr[i]

        method(dtmp,tmp,dt0,dt,t)
        @. tmp += dtmp

        #Filter_Distribution!(tmp,PhaseSpace,filter_vector)

        # removing negative values
        @. tmp = tmp*(tmp>=0f0)
        # hackey fix for inf values
        @. tmp = tmp*(tmp!=Inf)

        # saving state
        if (i in save_idx)
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

    if !isnothing(fileName)
        # save output to file
        println("Saving Simulation output")

        SolutionFileSave(output,PhaseSpace,fileLocation,fileName)

    end


    return output

end

