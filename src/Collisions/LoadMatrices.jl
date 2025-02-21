function LoadMatrices_Binary(Matrices_BinaryInteraction::Dict{Vector{String},Tuple},DataDirectory::String,Lists::ListStruct;mode="AXI")

    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list
    interaction_list_Binary = Lists.interaction_list_Binary

    if isempty(interaction_list_Binary) # no binary interactions to load
        return
    end

    for i in eachindex(interaction_list_Binary)

        interaction = interaction_list_Binary[i]

        name1 = interaction[1]
        name2 = interaction[2]
        name3 = interaction[3]
        name4 = interaction[4]

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        name3_loc = findfirst(==(name3),name_list)
        name4_loc = findfirst(==(name4),name_list)

        p1_grid = p_grid_list[name1_loc]
        p2_grid = p_grid_list[name2_loc]
        p3_grid = p_grid_list[name3_loc]
        p4_grid = p_grid_list[name4_loc]
        u1_grid = u_grid_list[name1_loc]
        u2_grid = u_grid_list[name2_loc]
        u3_grid = u_grid_list[name3_loc]
        u4_grid = u_grid_list[name4_loc]

        p1_num = p_num_list[name1_loc]
        p2_num = p_num_list[name2_loc]
        p3_num = p_num_list[name3_loc]
        p4_num = p_num_list[name4_loc]
        u1_num = u_num_list[name1_loc]
        u2_num = u_num_list[name2_loc]
        u3_num = u_num_list[name3_loc]
        u4_num = u_num_list[name4_loc]

        p1_num_st = string(p_num_list[name1_loc])
        p2_num_st = string(p_num_list[name2_loc])
        p3_num_st = string(p_num_list[name3_loc])
        p4_num_st = string(p_num_list[name4_loc])
        u1_num_st = string(u_num_list[name1_loc])
        u2_num_st = string(u_num_list[name2_loc])
        u3_num_st = string(u_num_list[name3_loc])
        u4_num_st = string(u_num_list[name4_loc])

        p1_low = p_low_list[name1_loc]
        p2_low = p_low_list[name2_loc]
        p3_low = p_low_list[name3_loc]
        p4_low = p_low_list[name4_loc]

        p1_low_st = string(p_low_list[name1_loc])
        p2_low_st = string(p_low_list[name2_loc])
        p3_low_st = string(p_low_list[name3_loc])
        p4_low_st = string(p_low_list[name4_loc])

        p1_up = p_up_list[name1_loc]
        p2_up = p_up_list[name2_loc]
        p3_up = p_up_list[name3_loc]
        p4_up = p_up_list[name4_loc]

        p1_up_st = string(p_up_list[name1_loc])
        p2_up_st = string(p_up_list[name2_loc])
        p3_up_st = string(p_up_list[name3_loc])
        p4_up_st = string(p_up_list[name4_loc])

        p4_r::Vector{Float64} = BCI.bounds(p4_low,p4_up,p4_num,p4_grid)
        p3_r::Vector{Float64} = BCI.bounds(p3_low,p3_up,p3_num,p3_grid)
        p2_r::Vector{Float64} = BCI.bounds(p2_low,p2_up,p2_num,p2_grid)
        p1_r::Vector{Float64} = BCI.bounds(p1_low,p1_up,p1_num,p1_grid)
        u4_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u4_num,u4_grid)
        u3_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u3_num,u3_grid)
        u2_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u2_num,u2_grid)
        u1_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u1_num,u1_grid)

        filename = name1*name2*name3*name4*"#"*p1_low_st*"-"*p1_up_st*p1_grid*p1_num_st*"#"*p2_low_st*"-"*p2_up_st*p2_grid*p2_num_st*"#"*p3_low_st*"-"*p3_up_st*p3_grid*p3_num_st*"#"*p4_low_st*"-"*p4_up_st*p4_grid*p4_num_st*"#"*u1_grid*u1_num_st*"#"*u2_grid*u2_num_st*"#"*u3_grid*u3_num_st*"#"*u4_grid*u4_num_st*".jld2"

        #Parameters = (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num)

        println(filename)

        if mode=="ISO" # ISO
            Parameters = BCI.fload_Matrix_ISO(DataDirectory,filename)[1] # 1 is Parameters
            matrices = BCI.fload_Matrix_ISO(DataDirectory,filename)[2:end] # 1 is Parameters
        elseif mode=="AXI" # AXI
            Parameters = BCI.fload_Matrix(DataDirectory,filename)[1] # 1 is Parameters
            matrices = BCI.fload_Matrix(DataDirectory,filename)[2:end] # 1 is Parameters
        else
            println("Mode not recognized")
        end

        if interaction[1] == interaction[2] && interaction[3] == interaction[4]
            # print conversion statistic
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            println("Scorrected:")
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],TMatrix1=matrices[2])
            
            #SMatrix = SixDtoThreeD(Float32.(matrices[1]))
            SMatrix = SixDtoTwoD(Float32.(matrices[1]))
            TMatrix = FourDtoTwoD(Float32.(matrices[2]))

            if mode=="ISO"
                SMatrix ./= 2*u1_num
                TMatrix ./= 2*u1_num
            end
            Matrices_BinaryInteraction[interaction] = (SMatrix,TMatrix)
        end
    
        if interaction[1] == interaction[2] && interaction[3] != interaction[4]

            # print conversion statistic
            BCI.DoesConserve(matrices[1],matrices[2],matrices[3],zeros(size(matrices[3])),Parameters)
            #SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3])
            
            # 3D matrices
            #SMatrix3 = SixDtoThreeD(Float32.(matrices[1]))
            #SMatrix4 = SixDtoThreeD(Float32.(matrices[2]))

            # 2D matrices
            SMatrix3 = SixDtoTwoD(Float32.(matrices[1]))
            SMatrix4 = SixDtoTwoD(Float32.(matrices[2]))
            TMatrix = FourDtoTwoD(Float32.(matrices[3]))
            Matrices_BinaryInteraction[interaction] = (SMatrix3,SMatrix4,TMatrix)
        end
    
        if interaction[1] != interaction[2] && interaction[3] == interaction[4]

            # print conversion statistic
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],matrices[3],Parameters)
            #SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],TMatrix1=matrices[2],TMatrix2=matrices[3])

            SMatrix = SixDtoTwoD(Float32.(matrices[1]))
            TMatrix1 = FourDtoTwoD(Float32.(matrices[2]))
            TMatrix2 = FourDtoTwoD(Float32.(matrices[3]))
            Matrices_BinaryInteraction[interaction] = (SMatrix,TMatrix1,TMatrix2)
        end
    
        if interaction[1] != interaction[2] && interaction[3] != interaction[4]
            BCI.DoesConserve(matrices[1],matrices[2],matrices[3],matrices[4],Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3],TMatrix2=matrices[4])

            SMatrix3 = SixDtoTwoD(Float32.(matrices[1]))
            SMatrix4 = SixDtoTwoD(Float32.(matrices[2]))
            TMatrix1 = FourDtoTwoD(Float32.(matrices[3]))
            TMatrix2 = FourDtoTwoD(Float32.(matrices[4]))
            Matrices_BinaryInteraction[interaction] = (SMatrix3,SMatrix4,TMatrix1,TMatrix2)
        end

    end # for

end


function SCorrection!(SMatrix3::Array{T,6},SMatrix4::Array{T,6},TMatrix1::Array{T,4},TMatrix2::Array{T,4},Parameters) where T

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num) = Parameters

    p1_r = BCI.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = BCI.deltaVector(p1_r);
    p1_d_full = [p1_d; BCI.deltaVector([p1_r[end]; 2*p1_r[end]])];
    u1_r = BCI.bounds(BCI.u_low,BCI.u_up,u1_num,u1_grid);
    u1_d = BCI.deltaVector(u1_r);

    p2_r = BCI.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = BCI.deltaVector(p2_r);
    p2_d_full = [p2_d; BCI.deltaVector([p2_r[end]; 2*p2_r[end]])];
    u2_r = BCI.bounds(BCI.u_low,BCI.u_up,u2_num,u2_grid);
    u2_d = BCI.deltaVector(u2_r);

    p3_r = BCI.bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = BCI.deltaVector(p3_r);
    p3_d_full = [p3_d; BCI.deltaVector([p3_r[end]; 2*p3_r[end]])];
    u3_r = BCI.bounds(BCI.u_low,BCI.u_up,u3_num,u3_grid);
    u3_d = BCI.deltaVector(u3_r);

    p4_r = BCI.bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = BCI.deltaVector(p4_r);
    p4_d_full = [p4_d; BCI.deltaVector([p4_r[end]; 2*p4_r[end]])];
    u4_r = BCI.bounds(BCI.u_low,BCI.u_up,u4_num,u4_grid);
    u4_d = BCI.deltaVector(u4_r);

    for k in axes(SMatrix3,3), l in axes(SMatrix3, 4), m in axes(SMatrix3,5), n in axes(SMatrix3,6)
        SsumN3 = zero(T)
        SsumN4 = zero(T)
        TsumN1 = zero(T)
        TsumN2 = zero(T)
        for i in axes(SMatrix3,1), j in axes(SMatrix3,2) 
        SsumN3 += SMatrix3[i,j,k,l,m,n]*p3_d_full[i]*u3_d[j]
        end
        for i in axes(SMatrix4,1), j in axes(SMatrix4,2) 
            SsumN4 += SMatrix4[i,j,k,l,m,n]*p4_d_full[i]*u4_d[j]
            end
        TsumN1 += TMatrix1[k,l,m,n]*p1_d_full[k]*u1_d[l]
        TsumN2 += TMatrix2[m,n,k,l]*p2_d_full[m]*u2_d[n]
        SCor = (TsumN1+TsumN2)/(SsumN3+SsumN4)
        @view(SMatrix3[:,:,k,l,m,n]) .*= SCor
        @view(SMatrix4[:,:,k,l,m,n]) .*= SCor
    end


end


function LoadMatrices_Emi(Matrices_Synchrotron::Dict{Vector{String},Array{Float32,2}},DataDirectory::String,Lists::ListStruct;mode="AXI")

    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list
    interaction_list_Emi = Lists.interaction_list_Emi

    if isempty(interaction_list_Emi) # no sync interactions to load
        return
    end

    for i in eachindex(interaction_list_Emi)

        interaction = interaction_list_Emi[i]

        name1 = "Pho"
        name2 = interaction[1]

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)

        p1_grid = p_grid_list[name1_loc]
        p2_grid = p_grid_list[name2_loc]
        u1_grid = u_grid_list[name1_loc]
        u2_grid = u_grid_list[name2_loc]

        p1_num = p_num_list[name1_loc]
        p2_num = p_num_list[name2_loc]
        u1_num = u_num_list[name1_loc]
        u2_num = u_num_list[name2_loc]

        p1_num_st = string(p_num_list[name1_loc])
        p2_num_st = string(p_num_list[name2_loc])
        u1_num_st = string(u_num_list[name1_loc])
        u2_num_st = string(u_num_list[name2_loc])

        p1_low = p_low_list[name1_loc]
        p2_low = p_low_list[name2_loc]

        p1_low_st = string(p_low_list[name1_loc])
        p2_low_st = string(p_low_list[name2_loc])

        p1_up = p_up_list[name1_loc]
        p2_up = p_up_list[name2_loc]

        p1_up_st = string(p_up_list[name1_loc])
        p2_up_st = string(p_up_list[name2_loc])

        p2_r::Vector{Float64} = BCI.bounds(p2_low,p2_up,p2_num,p2_grid)
        p1_r::Vector{Float64} = BCI.bounds(p1_low,p1_up,p1_num,p1_grid)
        u2_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u2_num,u2_grid)
        u1_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u1_num,u1_grid)

        filename = "sync"*name2*"#"*p1_low_st*"-"*p1_up_st*p1_grid*p1_num_st*"#"*p2_low_st*"-"*p2_up_st*p2_grid*p2_num_st*"#"*u1_grid*u1_num_st*"#"*u2_grid*u2_num_st*".jld2";

        println(filename)

        if mode=="ISO" # ISO
            Parameters = BCI.fload_Matrix_SyncISO(DataDirectory,filename)[1] # 1 is Parameters
            matrix = BCI.fload_Matrix_SyncISO(DataDirectory,filename)[2] # 1 is Parameters
        elseif mode=="AXI" # AXI
            Parameters = BCI.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
            matrix = BCI.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
            matrix = BCI.fload_Matrix_Sync(DataDirectory,filename) # remove later
        else
            println("Mode not recognized")
        end
            
        # some SMatrix values are greater than float32 precision!
        PhaseSpaceFactors_Sync_Undo!(matrix,p1_r,u1_r,p2_r,u2_r)
        SMatrix = FourDtoTwoD(Float32.(matrix)) / 3f8 # NOTE THIS FACTOR OF C IS UNCONFIRMED

        if mode=="ISO"
            SMatrix ./= 2*u1_num
        end

        Matrices_Synchrotron[interaction] = SMatrix

    end # for

end


#= ======================================================== =#

function LoadMatrices_Binary_Struct(BigM::BigMatrices,DataDirectory::String,Lists::ListStruct;mode="AXI")

    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list
    interaction_list_Binary = Lists.interaction_list_Binary

    if isempty(interaction_list_Binary) # no binary interactions to load
        return
    end

    fill!(BigM.A_Binary,0f0);

    for i in eachindex(interaction_list_Binary)

        interaction = interaction_list_Binary[i]

        name1_loc = findfirst(==(interaction[1]),name_list)
        name2_loc = findfirst(==(interaction[2]),name_list)
        name3_loc = findfirst(==(interaction[3]),name_list)
        name4_loc = findfirst(==(interaction[4]),name_list)

        name1 = interaction[1]
        name2 = interaction[2]
        name3 = interaction[3]
        name4 = interaction[4]

        # ele pos swap for compton
        if interaction == ["Pos","Pho","Pos","Pho"]
            name1 = "Ele"
            name2 = "Pho"
            name3 = "Ele"
            name4 = "Pho"
        end

        p1_grid = p_grid_list[name1_loc]
        p2_grid = p_grid_list[name2_loc]
        p3_grid = p_grid_list[name3_loc]
        p4_grid = p_grid_list[name4_loc]
        u1_grid = u_grid_list[name1_loc]
        u2_grid = u_grid_list[name2_loc]
        u3_grid = u_grid_list[name3_loc]
        u4_grid = u_grid_list[name4_loc]

        p1_num = p_num_list[name1_loc]
        p2_num = p_num_list[name2_loc]
        p3_num = p_num_list[name3_loc]
        p4_num = p_num_list[name4_loc]
        u1_num = u_num_list[name1_loc]
        u2_num = u_num_list[name2_loc]
        u3_num = u_num_list[name3_loc]
        u4_num = u_num_list[name4_loc]

        p1_num_st = string(p_num_list[name1_loc])
        p2_num_st = string(p_num_list[name2_loc])
        p3_num_st = string(p_num_list[name3_loc])
        p4_num_st = string(p_num_list[name4_loc])
        u1_num_st = string(u_num_list[name1_loc])
        u2_num_st = string(u_num_list[name2_loc])
        u3_num_st = string(u_num_list[name3_loc])
        u4_num_st = string(u_num_list[name4_loc])

        p1_low = p_low_list[name1_loc]
        p2_low = p_low_list[name2_loc]
        p3_low = p_low_list[name3_loc]
        p4_low = p_low_list[name4_loc]

        p1_low_st = string(p_low_list[name1_loc])
        p2_low_st = string(p_low_list[name2_loc])
        p3_low_st = string(p_low_list[name3_loc])
        p4_low_st = string(p_low_list[name4_loc])

        p1_up = p_up_list[name1_loc]
        p2_up = p_up_list[name2_loc]
        p3_up = p_up_list[name3_loc]
        p4_up = p_up_list[name4_loc]

        p1_up_st = string(p_up_list[name1_loc])
        p2_up_st = string(p_up_list[name2_loc])
        p3_up_st = string(p_up_list[name3_loc])
        p4_up_st = string(p_up_list[name4_loc])

        p4_r::Vector{Float64} = BCI.bounds(p4_low,p4_up,p4_num,p4_grid)
        p3_r::Vector{Float64} = BCI.bounds(p3_low,p3_up,p3_num,p3_grid)
        p2_r::Vector{Float64} = BCI.bounds(p2_low,p2_up,p2_num,p2_grid)
        p1_r::Vector{Float64} = BCI.bounds(p1_low,p1_up,p1_num,p1_grid)
        u4_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u4_num,u4_grid)
        u3_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u3_num,u3_grid)
        u2_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u2_num,u2_grid)
        u1_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u1_num,u1_grid)

        filename = name1*name2*name3*name4*"#"*p1_low_st*"-"*p1_up_st*p1_grid*p1_num_st*"#"*p2_low_st*"-"*p2_up_st*p2_grid*p2_num_st*"#"*p3_low_st*"-"*p3_up_st*p3_grid*p3_num_st*"#"*p4_low_st*"-"*p4_up_st*p4_grid*p4_num_st*"#"*u1_grid*u1_num_st*"#"*u2_grid*u2_num_st*"#"*u3_grid*u3_num_st*"#"*u4_grid*u4_num_st*".jld2"

        #Parameters = (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num)

        println(filename)

        if mode=="ISO" # ISO
            Parameters = BCI.fload_Matrix_ISO(DataDirectory,filename)[1] # 1 is Parameters
            matrices = BCI.fload_Matrix_ISO(DataDirectory,filename)[2:end] # 1 is Parameters
        elseif mode=="AXI" # AXI
            Parameters = BCI.fload_Matrix(DataDirectory,filename)[1] # 1 is Parameters
            matrices = BCI.fload_Matrix(DataDirectory,filename)[2:end] # 1 is Parameters
        else
            println("Mode not recognized")
        end

        if interaction[1] == interaction[2] && interaction[3] == interaction[4]
            # print conversion statistic
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            println("Scorrected:")
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],TMatrix1=matrices[2])

            Fill_A_Binary!(BigM.A_Binary,interaction,Lists;SMatrix3=matrices[1],TMatrix1=matrices[2])

        end
    
        if interaction[1] == interaction[2] && interaction[3] != interaction[4]

            # print conversion statistic
            BCI.DoesConserve(matrices[1],matrices[2],matrices[3],zeros(size(matrices[3])),Parameters)
            #SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3])
            
            Fill_A_Binary!(BigM.A_Binary,interaction,Lists;SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3])


        end
    
        if interaction[1] != interaction[2] && interaction[3] == interaction[4]

            # print conversion statistic
            BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],matrices[3],Parameters)
            #SCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #BCI.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],TMatrix1=matrices[2],TMatrix2=matrices[3])

            Fill_A_Binary!(BigM.A_Binary,interaction,Lists;SMatrix3=matrices[1],TMatrix1=matrices[2],TMatrix2=matrices[3])

        end
    
        if interaction[1] != interaction[2] && interaction[3] != interaction[4]
            BCI.DoesConserve(matrices[1],matrices[2],matrices[3],matrices[4],Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3],TMatrix2=matrices[4])
            
            Fill_A_Binary!(BigM.A_Binary,interaction,Lists;SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3],TMatrix2=matrices[4])

        end

    end # for

end

function LoadMatrices_Emi_Struct(BigM::BigMatrices,DataDirectory::String,Lists::ListStruct;mode="AXI")

    name_list = Lists.name_list
    p_up_list = Lists.p_up_list
    p_low_list = Lists.p_low_list
    p_grid_list = Lists.p_grid_list
    p_num_list = Lists.p_num_list
    u_grid_list = Lists.u_grid_list
    u_num_list = Lists.u_num_list
    interaction_list_Emi = Lists.interaction_list_Emi

    if isempty(interaction_list_Emi) # no sync interactions to load
        return
    end

    fill!(BigM.A_Emi,0f0);

    for i in eachindex(interaction_list_Emi)

        interaction = interaction_list_Emi[i]

        name1 = "Pho"
        name2 = interaction[1]
        #= to be implimented
        name1 = interaction[1]
        name2 = interaction[2]
        name3 = interaction[3]
        =#

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        # name3_loc = findfirst(==(name3),name_list)

        p1_grid = p_grid_list[name1_loc]
        p2_grid = p_grid_list[name2_loc]
        # p3_grid = p_grid_list[name3_loc]
        u1_grid = u_grid_list[name1_loc]
        u2_grid = u_grid_list[name2_loc]
        # u3_grid = u_grid_list[name3_loc]

        p1_num = p_num_list[name1_loc]
        p2_num = p_num_list[name2_loc]
        # p3_num = p_num_list[name3_loc]
        u1_num = u_num_list[name1_loc]
        u2_num = u_num_list[name2_loc]
        # u3_num = u_num_list[name3_loc]

        p1_num_st = string(p_num_list[name1_loc])
        p2_num_st = string(p_num_list[name2_loc])
        # p3_num_st = string(p_num_list[name3_loc])
        u1_num_st = string(u_num_list[name1_loc])
        u2_num_st = string(u_num_list[name2_loc])
        # u3_num_st = string(u_num_list[name3_loc])

        p1_low = p_low_list[name1_loc]
        p2_low = p_low_list[name2_loc]
        # p3_low = p_low_list[name3_loc]

        p1_low_st = string(p_low_list[name1_loc])
        p2_low_st = string(p_low_list[name2_loc])
        # p3_low_st = string(p_low_list[name3_loc])

        p1_up = p_up_list[name1_loc]
        p2_up = p_up_list[name2_loc]
        # p3_up = p_up_list[name3_loc]

        p1_up_st = string(p_up_list[name1_loc])
        p2_up_st = string(p_up_list[name2_loc])
        # p3_up_st = string(p_up_list[name3_loc])

        # p3_r::Vector{Float64} = BCI.bounds(p3_low,p3_up,p3_num,p3_grid)
        p2_r::Vector{Float64} = BCI.bounds(p2_low,p2_up,p2_num,p2_grid)
        p1_r::Vector{Float64} = BCI.bounds(p1_low,p1_up,p1_num,p1_grid)
        # u3_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u3_num,u3_grid)
        u2_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u2_num,u2_grid)
        u1_r::Vector{Float64} = BCI.bounds(BCI.u_low,BCI.u_up,u1_num,u1_grid)

        filename = "sync"*name2*"#"*p1_low_st*"-"*p1_up_st*p1_grid*p1_num_st*"#"*p2_low_st*"-"*p2_up_st*p2_grid*p2_num_st*"#"*u1_grid*u1_num_st*"#"*u2_grid*u2_num_st*".jld2";

        println(filename)

        if mode=="ISO" # ISO
            Parameters = BCI.fload_Matrix_SyncISO(DataDirectory,filename)[1] # 1 is Parameters
            matrix = BCI.fload_Matrix_SyncISO(DataDirectory,filename)[2] # 1 is Parameters
        elseif mode=="AXI" # AXI
            Parameters = BCI.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
            matrix = BCI.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
            matrix = BCI.fload_Matrix_Sync(DataDirectory,filename) # remove later
        else
            println("Mode not recognized")
        end
            
        # some SMatrix values are greater than float32 precision!
        PhaseSpaceFactors_Sync_Undo!(matrix,p1_r,u1_r,p2_r,u2_r)
        SMatrix = Float32.(matrix) / 3f8 # NOTE THIS FACTOR OF C IS UNCONFIRMED

        Fill_A_Emi!(BigM.A_Emi,interaction,Lists;SMatrix3=SMatrix)

    end # for

end
