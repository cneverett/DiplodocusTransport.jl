function Solve(method::SteppingMethodType;save_steps::Int=1,progress::Bool=false,fileName::String=nothing,fileLocation::String=pwd(),Verbose::Bool=false)

    if isdir(fileLocation) == false
        mkpath(fileLocation)
    end

    PhaseSpace = method.PhaseSpace
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids

    f = copy(method.f_init)
    tr = Grids.tr
    dt0 = tr[2] - tr[1]
    t_num = length(tr)-1

    save_idx = unique([1 ; range(save_steps,t_num,step=save_steps)]) # saves first time step and then on an interval of save_steps (no duplicates)

    n_save = length(save_idx)+1 # +1 for the initial state

    output = OutputStruct(f,n_save)

    # save initial state (step 1)
    output.f[1] = copy(f)
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

        method(dt0,dt,t,Verbose=Verbose)
        #@. tmp += dtmp # now done in stepping method

        #Filter_Distribution!(tmp,PhaseSpace,filter_vector)

        # removing negative values (values less than 1f-28 for better stability) 
        # now done in stepping method
        #@. tmp = tmp*(tmp>=1f-28)
        # hacky fix for inf values
        #@. tmp = tmp*(tmp!=Inf)

        # saving state
        if (i in save_idx)
            save_count += 1
            t = tr[i+1]
            copyto!(f,method.f)
            #println("$(method.f)")
            output.f[save_count] = copy(f)
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

