function Solver(f0,timespan,Lists,BigM::BigMatrices,dt; method=Euler() #= method=ImplicitEuler(autodiff=false) =#,save_dt=dt)

    deriv_cpu = BoltzmannEquation(f0,Lists,BigM,dt);

    prob = ODEProblem(deriv_cpu,f0,timespan)

    @time solution = solve(prob,method,tstops=timespan[1]:dt:timespan[2],saveat=save_dt, maxiters = 1e5,progress=true,progress_steps = round(timespan[2]/save_dt)#=, isoutofdomain = (u,p,t)->any(x->x<0,u)=#)

    return solution

end

function (g::BoltzmannEquation)(df,f,p,t)

    #=
    |---- . ----|---- . ----|
    f0    uR1   u1    uA1   u2
    =#

    # limit u to be positive
    @. f = f*(f>0f0)
    fill!(df,Float32(0))

    # reset arrays
    fill!(g.Δf,Float32(0))
    fill!(g.J,Float32(0))
    
    # update changes due to Binary S and T interactions
    if isempty(g.interaction_list_Binary) == false # no binary interactions to load
    #    update_ΔSΔT_Binary!(g,Matrices_BinaryInteraction,f)
    #    @. df += g.ΔfS_list - g.ΔfT_list
        update_Big_Binary!(g,g.BigM,f)
        @. df = g.Δf
    end

    if isempty(g.interaction_list_Emi) == false
        update_ΔS_Sync!(g,Matrices_Synchrotron,f)
        @. df += g.ΔfS_list
    end

    # explicit Stepping
    #@. df = g.Δf

    # implicit Stepping 
    df .= (I-g.dt.*g.J)\g.Δf

    # leapfrog in time 
        if t==0
            @. g.f1DR = f
        end

        @. df *= 2
        @. df += g.f1DR - f
        # set f i-1 for next time step   
        @. g.f1DR = f


end

function update_Big_Binary!(g::BoltzmannEquation,BigM::BigMatrices,f)
  
    mul!(g.A_Binary_Reshape,BigM.A_Binary,f)
    temp = reshape(g.A_Binary_Reshape,size(BigM.A_Binary,2),size(BigM.A_Binary,2))
    g.J += 2*temp # assign jacobian elements
    mul!(g.Δf_temp,temp,f)
    g.Δf += g.Δf_temp 

end

function update_ΔSΔT_Binary!(g::BoltzmannEquation,Matrices_BinaryInteraction,f)

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
        matrices = Matrices_BinaryInteraction[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)
        name3 = interaction[3]
        name3_loc = findfirst(==(name3),g.name_list)
        name4 = interaction[4]
        name4_loc = findfirst(==(name4),g.name_list)

        if (name1 == name2) && (name3 == name4)
            SMatrix = matrices[1]
            TMatrix = matrices[2]
  
            # 3D Smatrix
                #@tullio threads=false g.ΔfS_list.x[name3_loc][i] += SMatrix[i,j,k] * f.x[name1_loc][j] * f.x[name2_loc][k] 

            # 2D SMatrix
            mul!(g.ΔfS_mul_step.x[i],SMatrix,f.x[name2_loc])
            temp_reshape = reshape(g.ΔfS_mul_step.x[i],g.nump_list[name3_loc]*g.numt_list[name3_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
            mul!(g.ΔfS_list_temp.x[name3_loc],temp_reshape,f.x[name1_loc])
            g.ΔfS_list.x[name3_loc] .+= g.ΔfS_list_temp.x[name3_loc]
            

            # TMatrix
            mul!(g.ΔfT_list_temp.x[name1_loc],TMatrix,f.x[name2_loc])
            g.ΔfT_list.x[name1_loc] .+= f.x[name1_loc] .* g.ΔfT_list_temp.x[name1_loc]
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matrices[1]
            SMatrix4 = matrices[2]
            TMatrix = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matrices[1]
            TMatrix1 = matrices[2]
            TMatrix2 = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
        =#
        if (name1 != name2) && (name3 != name4)

            SMatrix3 = matrices[1]
            SMatrix4 = matrices[2]
            TMatrix1 = matrices[3]
            TMatrix2 = matrices[4]

            # 3D SMatrix
                # @tullio doesn't seem to work with different sized PartitionedArrays
                #@tullio threads=false g.ΔfS_list.x[name3_loc][i] += SMatrix3[i,j,k] * f.x[name1_loc][j] * f.x[name2_loc][k]
                #@tullio threads=false g.ΔfS_list.x[name4_loc][i] += SMatrix4[i,j,k] * f.x[name2_loc][j] * f.x[name1_loc][k] 
                #=@turbo for i in axes(SMatrix3,1), j in axes(SMatrix3,2), k in axes(SMatrix3,3)
                    g.ΔfS_list.x[name3_loc][i] += SMatrix3[i,j,k] * f.x[name1_loc][j] * f.x[name2_loc][k]
                end
                @turbo for i in axes(SMatrix4,1), j in axes(SMatrix4,2), k in axes(SMatrix4,3)
                    g.ΔfS_list.x[name4_loc][i] += SMatrix4[i,j,k] * f.x[name1_loc][j] * f.x[name2_loc][k]
                end
                g.ΔfT_list.x[name1_loc] .+= f.x[name1_loc] .* (TMatrix1 * f.x[name2_loc])
                g.ΔfT_list.x[name2_loc] .+= f.x[name2_loc] .* (TMatrix2 * f.x[name1_loc])
            =#

            # 2D SMatrix
                # SMatrix3
                mul!(g.ΔfS_mul_step.x[i][1],SMatrix3,f.x[name2_loc])
                temp_reshape = reshape(g.ΔfS_mul_step.x[i][1],g.nump_list[name3_loc]*g.numt_list[name3_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
                mul!(g.ΔfS_list_temp.x[name3_loc],temp_reshape,f.x[name1_loc])
                g.ΔfS_list.x[name3_loc] .+= g.ΔfS_list_temp.x[name3_loc]
                # SMatrix4
                mul!(g.ΔfS_mul_step.x[i][2],SMatrix4,f.x[name2_loc])
                temp_reshape = reshape(g.ΔfS_mul_step.x[i][2],g.nump_list[name4_loc]*g.numt_list[name4_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
                mul!(g.ΔfS_list_temp.x[name4_loc],temp_reshape,f.x[name1_loc])
                g.ΔfS_list.x[name4_loc] .+= g.ΔfS_list_temp.x[name4_loc]
                # TMatrix1
                mul!(g.ΔfT_list_temp.x[name1_loc],TMatrix1,f.x[name2_loc])
                g.ΔfT_list.x[name1_loc] .+= f.x[name1_loc] .* g.ΔfT_list_temp.x[name1_loc]
                # TMatrix2
                mul!(g.ΔfT_list_temp.x[name2_loc],TMatrix2,f.x[name1_loc])
                g.ΔfT_list.x[name2_loc] .+= f.x[name2_loc] .* g.ΔfT_list_temp.x[name2_loc]
        end

    end

end


function update_ΔS_Binary!(g::BoltzmannEquation,Matrices_BinaryInteraction,f)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfS_list,Float32(0))

    for i in eachindex(g.interaction_list_Binary)
        interaction = g.interaction_list_Binary[i]
        matrices = Matrices_BinaryInteraction[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)
        name3 = interaction[3]
        name3_loc = findfirst(==(name3),g.name_list)
        name4 = interaction[4]
        name4_loc = findfirst(==(name4),g.name_list)

        if (name1 == name2) && (name3 == name4)
            SMatrix = matrices[1]
  
            mul!(g.ΔfS_mul_step.x[i],SMatrix,f.x[name2_loc])
            temp_reshape = reshape(g.ΔfS_mul_step.x[i],g.nump_list[name3_loc]*g.numt_list[name3_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
            mul!(g.ΔfS_list_temp.x[name3_loc],temp_reshape,f.x[name1_loc])
            g.ΔfS_list.x[name3_loc] .+= g.ΔfS_list_temp.x[name3_loc]
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matrices[1]
            SMatrix4 = matrices[2]
            TMatrix = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matrices[1]
            TMatrix1 = matrices[2]
            TMatrix2 = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
        =#
    
        if (name1 != name2) && (name3 != name4)
            SMatrix3 = matrices[1]
            SMatrix4 = matrices[2]
            # SMatrix3
            mul!(g.ΔfS_mul_step.x[i][1],SMatrix3,f.x[name2_loc])
            temp_reshape = reshape(g.ΔfS_mul_step.x[i][1],g.nump_list[name3_loc]*g.numt_list[name3_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
            mul!(g.ΔfS_list_temp.x[name3_loc],temp_reshape,f.x[name1_loc])
            g.ΔfS_list.x[name3_loc] .+= g.ΔfS_list_temp.x[name3_loc]
            # SMatrix4
            mul!(g.ΔfS_mul_step.x[i][2],SMatrix4,f.x[name2_loc])
            temp_reshape = reshape(g.ΔfS_mul_step.x[i][2],g.nump_list[name4_loc]*g.numt_list[name4_loc],g.nump_list[name1_loc]*g.numt_list[name1_loc])
            mul!(g.ΔfS_list_temp.x[name4_loc],temp_reshape,f.x[name1_loc])
        end

    end

end

function update_ΔT_Binary!(g::BoltzmannEquation,Matrices_BinaryInteraction,f)

    #f_list = g.f_list
    #interaction_list = g.interaction_list
    #name_list = g.name_list
    #ΔfS_list = g.ΔfS_list
    #ΔfT_list = g.ΔfT_list

    # reset arrays
    fill!(g.ΔfT_list,Float32(0))

    for i in eachindex(g.interaction_list_Binary)
        interaction = g.interaction_list_Binary[i]
        matrices = Matrices_BinaryInteraction[interaction]
        name1 = interaction[1]
        name1_loc = findfirst(==(name1),g.name_list)
        name2 = interaction[2]
        name2_loc = findfirst(==(name2),g.name_list)


        if (name1 == name2) && (name3 == name4)
            TMatrix = matrices[2]
  
            mul!(g.ΔfT_list_temp.x[name1_loc],TMatrix,f.x[name2_loc])
            g.ΔfT_list.x[name1_loc] .+= f.x[name1_loc] .* g.ΔfT_list_temp.x[name1_loc]
            
        end
    
        #=if (name1 == name2) && (name3 != name4)
            SMatrix3 = matrices[1]
            SMatrix4 = matrices[2]
            TMatrix = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix3 * f_list[name2_loc])' * f_list[name1_loc]
            ΔS_list[name4_loc] .= (SMatrix4 * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix * f_list[name2_loc]) .* f_list[name1_loc]
        end
    
        if (name1 != name2) && (name3 == name4)
            SMatrix = matrices[1]
            TMatrix1 = matrices[2]
            TMatrix2 = matrices[3]
            ΔS_list[name3_loc] .= (SMatrix * f_list[name2_loc])' * f_list[name1_loc]
            ΔT_list[name1_loc] .= (TMatrix1 * f_list[name2_loc]) .* f_list[name1_loc]
            ΔT_list[name2_loc] .= (TMatrix2 * f_list[name1_loc]) .* f_list[name2_loc]
        end
        =#
    
        if (name1 != name2) && (name3 != name4)
            TMatrix1 = matrices[3]
            TMatrix2 = matrices[4]
            # TMatrix1
            mul!(g.ΔfT_list_temp.x[name1_loc],TMatrix1,f.x[name2_loc])
            g.ΔfT_list.x[name1_loc] .+= f.x[name1_loc] .* g.ΔfT_list_temp.x[name1_loc]
            # TMatrix2
            mul!(g.ΔfT_list_temp.x[name2_loc],TMatrix2,f.x[name1_loc])
            g.ΔfT_list.x[name2_loc] .+= f.x[name2_loc] .* g.ΔfT_list_temp.x[name2_loc]
        end

    end

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
