function Solve(f1D0::fType,method::SteppingMethod;save_steps::Int=1,progress=false)

    f0 = copy(f1D0)
    t_low = method.SpaceTime.t_low
    t_up = method.SpaceTime.t_up
    dt = method.SpaceTime.dt

    tsteps = t_low:dt:t_up
    nsteps = length(tsteps)
    nsave = round(Int64,nsteps/save_steps)+1

    tmp = copy(f0)
    dtmp = similar(f0)

    t = t_low

    output = SolutionOutput(f0,t,nsave)

    # save initial state
    output.f[1] = copy(tmp)
    output.t[1] = t
    save_count = 2

    # progress bar
    if progress
        p = Progress(nsteps-1)
    end

    for i in 2:nsteps

        t = tsteps[i]

        method(dtmp,tmp,t,dt)
        @. tmp += dtmp

        # removing values less than 0f0
        @. tmp = tmp*(tmp>=0f0)

        # saving state
        if (i-1)%save_steps == 0
            output.f[save_count] = copy(tmp)
            output.t[save_count] = t
            save_count += 1
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

function update_ΔS_Sync!(g::BoltzmannEquation,Matrices_Synchrotron,f)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfS_list,Float32(0))

    for i in eachindex(g.interaction_list_Sync)
        interaction = g.interaction_list_Sync[i]
        SMatrix = Matrices_Synchrotron[interaction]
        name1_loc = findfirst(==("Pho"),g.name_list)
        name2 = interaction[1]
        name2_loc = findfirst(==(name2),g.name_list)

        mul!(g.ΔfS_list_temp.x[name1_loc],SMatrix,f.x[name2_loc])
  
        g.ΔfS_list.x[name1_loc] .+= g.ΔfS_list_temp.x[name1_loc]
        
    end

end
