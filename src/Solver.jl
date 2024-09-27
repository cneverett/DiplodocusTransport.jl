function BoltzmannEquationSolver(f1D0,timespan,Lists,dt; method=Euler() #==ImplicitEuler(autodiff=false)=#,save_dt=dt)

    deriv_cpu = BoltzmannEquation(f1D0,Lists);

    prob = ODEProblem(deriv_cpu,f1D0,timespan)

    @time solution = solve(prob,method,tstops=timespan[1]:dt:timespan[2],saveat=save_dt, maxiters = 1e5,progress=true,progress_steps = round(timespan[2]/save_dt)#=, isoutofdomain = (u,p,t)->any(x->x<0,u)=#)

    return solution

end

function (g::BoltzmannEquation)(df1D,f1D,p,t)

    #=
    |---- . ----|---- . ----|
    f1D0    uR1   u1    uA1   u2
    =#

    # limit u to be positive
    @. f1D = f1D*(f1D>0f0)
    fill!(df1D,Float32(0))

    #= # 1/2 T -> 1 S -> 1/2 T
    if t == 0 # set retarded point
        dt = 1f-1/10000
        update_ΔSΔT!(g,CollisionMatriciesBinary,u)
        @. du = g.ΔfS_list - g.ΔfT_list
        @. g.uR = (2*u + du * dt)/2
    else
        dt = (t - g.t)
        update_ΔT!(g,CollisionMatriciesBinary,u)
        # half time step forward with T 
        @. g.uA = u - g.ΔfT_list * dt/2
        @. g.uA = g.uA*(g.uA>0f0)
        # average with previous time step
        @. g.uA = (g.uR + g.uA)/2
        # advance by full time step with S
        update_ΔS!(g,CollisionMatriciesBinary,g.uA) 
        @. g.uA = g.uA + g.ΔfS_list * dt
        # average with previous time step
        @. g.uR = (u + g.uA)/2 # also moves uA
        # advance by half step with T
        update_ΔT!(g,CollisionMatriciesBinary,g.uR)
        utest = g.uR - g.ΔfT_list * dt/2
        @. utest = utest*(utest>0f0) # prevents -ve u in sol
        @. du = (utest-u)/dt
    end
    =#
    
    
    
    #= # 1/2 S -> 1/2 T -> 1/2 S -> 1/2 T
    if t == 0 # set retarded point
        dt = 1f-1/1000
        update_ΔSΔT!(g,CollisionMatriciesBinary,u)
        @. du = g.ΔfS_list - g.ΔfT_list
        @. g.uR = (2*u + du * dt)/2
    else
        dt = (t - g.t)
        update_ΔS!(g,CollisionMatriciesBinary,u)
        # half time step forward with S 
        @. g.uA = u + g.ΔfS_list * dt/2
        # average with previous time step to return to initial time
        @. g.uA = (g.uR + g.uA)/2
        # half time step forward with T
        update_ΔT!(g,CollisionMatriciesBinary,g.uA) 
        @. g.uA = g.uA - g.ΔfT_list * dt/2
        @. g.uR = g.uA*(g.uA>0f0) # sets g.uR for next time step
        # half time ste forward with S
        update_ΔS!(g,CollisionMatriciesBinary,g.uR)
        @. g.uA = g.uR + g.ΔfS_list * dt/2
        # average with previous time step to return to initial half step
        @. g.uA = (u + g.uA)/2
        # half time step forward with T
        update_ΔT!(g,CollisionMatriciesBinary,g.uA)
        utest = g.uR - g.ΔfT_list * dt/2
        @. utest = utest*(utest>0f0) # prevents -ve u in sol
        @. du = (utest-u)/dt
    end
    =#

    #= 1 S -> 1 T 
    if t == 0 # set retarded point
        dt = 1f-1/1000
        update_ΔSΔT!(g,CollisionMatriciesBinary,u)
        @. du = g.ΔfS_list - g.ΔfT_list
        @. g.uR = u
    else
        dt = (t - g.t)
        update_ΔS!(g,CollisionMatriciesBinary,u)
        # time step forward with S 
        @. g.uA = u + g.ΔfS_list * dt
        # average with previous time step to return to initial time
        @. g.uA = (g.uR + g.uA)/2
        # time step forward with T
        update_ΔT!(g,CollisionMatriciesBinary,g.uA) 
        @. g.uA = g.uA - g.ΔfT_list * dt
        @. g.uA = g.uA*(g.uA>0f0)
        @. du = (g.uA-u)/dt

        @. g.uR = u
    end
    =#
    
    #=# 1 S+T -> 1 T for new values
    if t == 0 # set retarded point
        dt = 1f-2/1000
    else
        dt = (t - g.t)
    end
    update_ΔSΔT!(g,CollisionMatriciesBinary,u)
    @. du = g.ΔfS_list - g.ΔfT_list
    @. g.uR = g.ΔfS_list*(u<=0f0)*dt/2 # extra values at initial time by interpolate
    # update T for extra values
    update_ΔT!(g,CollisionMatriciesBinary,g.uR)
    #@. g.uA = (du - g.ΔfT_list)*dt
    #@. g.uA = g.uA*(g.uA>0f0) # prevents -ve u in sol
    @. du = du - g.ΔfT_list
    =#
   
    
    

    # update dt
    g.t = t

    #@. g.uA = u

    #@. u = (g.uA+g.uR)/2
    
    # update changes due to Binary S and T interactions
    if isempty(g.interaction_list_Binary) == false # no sync interactions to load
        update_ΔSΔT_Binary!(g,CollisionMatriciesBinary,f1D)
        @. df1D += g.ΔfS_list - g.ΔfT_list
    end

    if isempty(g.interaction_list_Sync) == false
        update_ΔS_Sync!(g,CollisionMatriciesSync,u)
        @. df1D += g.ΔfS_list
    end

    # update uR as point from this step
    #@. g.uR = u

    #t=t-dt/2
    #update t
    
    #g.state = !(g.state)


end

function update_ΔSΔT_Binary!(g::BoltzmannEquation,CollisionMatriciesBinary,f1D)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfS_list,Float32(0))
    fill!(g.ΔfT_list,Float32(0))

    for i in eachindex(g.interaction_list_Binary)
        interaction = g.interaction_list_Binary[i]
        matricies = CollisionMatriciesBinary[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)
        name3 = interaction[3]
        name3_loc = findfirst(==(name3),g.name_list)
        name4 = interaction[4]
        name4_loc = findfirst(==(name4),g.name_list)

        if (name1 == name2) && (name3 == name4)
            SMatrix = matricies[1]
            TMatrix = matricies[2]
  
            @tullio threads=false g.ΔfS_list.x[name3_loc][i] += SMatrix[i,j,k] * f1D.x[name2_loc][j] * f1D.x[name1_loc][k] 
            g.ΔfT_list.x[name1_loc] .+= f1D.x[name2_loc] .* (TMatrix * f1D.x[name1_loc])
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matricies[1]
            TMatrix1 = matricies[2]
            TMatrix2 = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
    
        if (name1 != name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix1 = matricies[3]
            TMatrix2 = matricies[4]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end=#

    end

end


function update_ΔS_Binary!(g::BoltzmannEquation,CollisionMatriciesBinary,f1D)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfS_list,Float32(0))

    for i in eachindex(g.interaction_list_Binary)
        interaction = g.interaction_list_Binary[i]
        matricies = CollisionMatriciesBinary[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)
        name3 = interaction[3]
        name3_loc = findfirst(==(name3),g.name_list)
        name4 = interaction[4]
        name4_loc = findfirst(==(name4),g.name_list)

        if (name1 == name2) && (name3 == name4)
            SMatrix = matricies[1]
  
            @tullio threads=false g.ΔfS_list.x[name3_loc][i] += SMatrix[i,j,k] * f1D.x[name2_loc][j] * f1D.x[name1_loc][k] 
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matricies[1]
            TMatrix1 = matricies[2]
            TMatrix2 = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
    
        if (name1 != name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix1 = matricies[3]
            TMatrix2 = matricies[4]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end=#

    end

end

function update_ΔT_Binary!(g::BoltzmannEquation,CollisionMatriciesBinary,f1D)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfT_list,Float32(0))

    for i in eachindex(g.interaction_list_Binary)
        interaction = g.interaction_list_Binary[i]
        matricies = CollisionMatriciesBinary[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)


        if (name1 == name2) && (name3 == name4)
            TMatrix = matricies[2]
  
            g.ΔfT_list.x[name1_loc] .+= f1D.x[name2_loc] .* (TMatrix * f1D.x[name1_loc])
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matricies[1]
            TMatrix1 = matricies[2]
            TMatrix2 = matricies[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
    
        if (name1 != name2) && (name3 != name4)
            SMatrix3 = matricies[1]
            SMatrix4 = matricies[2]
            TMatrix1 = matricies[3]
            TMatrix2 = matricies[4]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end=#

    end

end

function update_ΔS_Sync!(g::BoltzmannEquation,CollisionMatriciesSync,f1D)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfS_list,Float32(0))

    for i in eachindex(g.interaction_list_Sync)
        interaction = g.interaction_list_Sync[i]
        matrix = CollisionMatriciesSync[interaction]
        name1_loc = findfirst(==("Pho"),g.name_list)
        name2 = interaction[1]
        name2_loc = findfirst(==(name2),g.name_list)

        SMatrix = matrix
  
        g.ΔfS_list.x[name1_loc] .+= (SMatrix * f1D.x[name2_loc])
        
    end

end
