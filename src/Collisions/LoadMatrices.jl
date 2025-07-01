function GainCorrection!(Output::Tuple{Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, Array{Float64, 9}, Array{Float64, 9}, Array{Float64, 6}, Array{Float64, 6}})

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Output[1]

    GainMatrix3 = Output[2]
    GainMatrix4 = Output[3]
    LossMatrix1 = Output[4]
    LossMatrix2 = Output[5]

    p1_r = DC.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = DC.deltaVector(p1_r);
    p1_d_full = [p1_d; DC.deltaVector([p1_r[end]; 2*p1_r[end]])];
    u1_r = DC.bounds(DC.u_low,DC.u_up,u1_num,u1_grid);
    u1_d = DC.deltaVector(u1_r);

    p2_r = DC.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = DC.deltaVector(p2_r);
    p2_d_full = [p2_d; DC.deltaVector([p2_r[end]; 2*p2_r[end]])];
    u2_r = DC.bounds(DC.u_low,DC.u_up,u2_num,u2_grid);
    u2_d = DC.deltaVector(u2_r);

    p3_r = DC.bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = DC.deltaVector(p3_r);
    p3_d_full = [p3_d; DC.deltaVector([p3_r[end]; 2*p3_r[end]])];
    u3_r = DC.bounds(DC.u_low,DC.u_up,u3_num,u3_grid);
    u3_d = DC.deltaVector(u3_r);

    p4_r = DC.bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = DC.deltaVector(p4_r);
    p4_d_full = [p4_d; DC.deltaVector([p4_r[end]; 2*p4_r[end]])];
    u4_r = DC.bounds(DC.u_low,DC.u_up,u4_num,u4_grid);
    u4_d = DC.deltaVector(u4_r);

    #=     for k in axes(SMatrix3,3), l in axes(SMatrix3, 4), m in axes(SMatrix3,5), n in axes(SMatrix3,6)
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
    end =#


    for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        
        GainSumN3 = zero(Float64)
        GainSumN4 = zero(Float64)
        LossSumN1 = LossMatrix1[p1,u1,h1,p2,u2,h2]
        LossSumN2 = LossMatrix2[p2,u2,h2,p1,u1,h1]
        
        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
            GainSumN3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        end
        for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
        GainSumN4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        end

        if name1 == name3
            if GainSumN3 != 0e0
                Correction = LossSumN1/GainSumN3
                @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) .*= Correction
            else
                GainMatrix3[p1,u1,h1,p1,u1,h1,p2,u2,h2] += LossSumN1
            end
            
        end

        if name2 == name4
            if GainSumN4 != 0e0
                Correction = LossSumN2/GainSumN4
                @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]) .*= Correction
            else
                GainMatrix4[p2,u2,h2,p1,u1,h1,p2,u2,h2] += LossSumN2
            end
        end

    end

    return nothing 

end

function LoadMatrices_Binary(M_Bin::Array{Float32,2},DataDirectory::String,PhaseSpace::PhaseSpaceStruct)

    Binary_list = PhaseSpace.Binary_list

    if isempty(Binary_list) # no binary interactions to load
        return
    end
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    fill!(M_Bin,0f0);

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        name3_loc = findfirst(==(name3),name_list)
        name4_loc = findfirst(==(name4),name_list)

        # ele pos swap for compton
        if name1 == "Pos" && name2 == "Pho" && name3 == "Pos" && name4 == "Pho"
            name1 = "Ele"
            name2 = "Pho"
            name3 = "Ele"
            name4 = "Pho"
        end

        px1_grid = px_grid_list[name1_loc]
        px2_grid = px_grid_list[name2_loc]
        px3_grid = px_grid_list[name3_loc]
        px4_grid = px_grid_list[name4_loc]
        py1_grid = py_grid_list[name1_loc]
        py2_grid = py_grid_list[name2_loc]
        py3_grid = py_grid_list[name3_loc]
        py4_grid = py_grid_list[name4_loc]
        pz1_grid = pz_grid_list[name1_loc]
        pz2_grid = pz_grid_list[name2_loc]
        pz3_grid = pz_grid_list[name3_loc]
        pz4_grid = pz_grid_list[name4_loc]

        px1_num = px_num_list[name1_loc]
        px2_num = px_num_list[name2_loc]
        px3_num = px_num_list[name3_loc]
        px4_num = px_num_list[name4_loc]
        py1_num = py_num_list[name1_loc]
        py2_num = py_num_list[name2_loc]
        py3_num = py_num_list[name3_loc]
        py4_num = py_num_list[name4_loc]
        pz1_num = pz_num_list[name1_loc]
        pz2_num = pz_num_list[name2_loc]
        pz3_num = pz_num_list[name3_loc]
        pz4_num = pz_num_list[name4_loc]

        px1_low = px_low_list[name1_loc]
        px2_low = px_low_list[name2_loc]
        px3_low = px_low_list[name3_loc]
        px4_low = px_low_list[name4_loc]
        #=py1_low = py_low_list[name1_loc]
        py2_low = py_low_list[name2_loc]
        py3_low = py_low_list[name3_loc]
        py4_low = py_low_list[name4_loc]
        pz1_low = pz_low_list[name1_loc]
        pz2_low = pz_low_list[name2_loc]
        pz3_low = pz_low_list[name3_loc]
        pz4_low = pz_low_list[name4_loc]=#
        
        px1_up = px_up_list[name1_loc]
        px2_up = px_up_list[name2_loc]
        px3_up = px_up_list[name3_loc]
        px4_up = px_up_list[name4_loc]
        #=py1_up = py_up_list[name1_loc]
        py2_up = py_up_list[name2_loc]
        py3_up = py_up_list[name3_loc]
        py4_up = py_up_list[name4_loc]
        pz1_up = pz_up_list[name1_loc]
        pz2_up = pz_up_list[name2_loc]
        pz3_up = pz_up_list[name3_loc]
        pz4_up = pz_up_list[name4_loc]=#

        m1 = PhaseSpace.Grids.mass_list[name1_loc]
        m2 = PhaseSpace.Grids.mass_list[name2_loc]
        m3 = PhaseSpace.Grids.mass_list[name3_loc]
        m4 = PhaseSpace.Grids.mass_list[name4_loc]

        pxr1::Vector{Float64} = PhaseSpace.Grids.pxr_list[name1_loc]
        pxr2::Vector{Float64} = PhaseSpace.Grids.pxr_list[name2_loc]
        pxr3::Vector{Float64} = PhaseSpace.Grids.pxr_list[name3_loc]
        pxr4::Vector{Float64} = PhaseSpace.Grids.pxr_list[name4_loc]
        pyr1::Vector{Float64} = PhaseSpace.Grids.pyr_list[name1_loc]
        pyr2::Vector{Float64} = PhaseSpace.Grids.pyr_list[name2_loc]
        pyr3::Vector{Float64} = PhaseSpace.Grids.pyr_list[name3_loc]
        pyr4::Vector{Float64} = PhaseSpace.Grids.pyr_list[name4_loc]
        pzr1::Vector{Float64} = PhaseSpace.Grids.pzr_list[name1_loc]
        pzr2::Vector{Float64} = PhaseSpace.Grids.pzr_list[name2_loc]
        pzr3::Vector{Float64} = PhaseSpace.Grids.pzr_list[name3_loc]
        pzr4::Vector{Float64} = PhaseSpace.Grids.pzr_list[name4_loc]

        #filename = name1*name2*name3*name4*"#"*px1_low*"-"*px1_up*px1_grid*px1_num*"#"*px2_low*"-"*px2_up*px2_grid*px2_num*"#"*px3_low*"-"*px3_up*px3_grid*px3_num*"#"*px4_low*"-"*px4_up*px4_grid*px4_num*"#"*py1_grid*py1_num*"#"*py2_grid*py2_num*"#"*py3_grid*py3_num*"#"*py4_grid*py4_num*".jld2"

        filename = name1*name2*name3*name4
        filename *= "#"*string(px1_low)*"-"*string(px1_up)*px1_grid*string(px1_num)
        filename *= "#"*py1_grid*string(py1_num)
        filename *= "#"*pz1_grid*string(pz1_num)
    
        filename *= "#"*string(px2_low)*"-"*string(px2_up)*px2_grid*string(px2_num)
        filename *= "#"*py2_grid*string(py2_num)
        filename *= "#"*pz2_grid*string(pz2_num)
    
        filename *= "#"*string(px3_low)*"-"*string(px3_up)*px3_grid*string(px3_num)
        filename *= "#"*py3_grid*string(py3_num)
        filename *= "#"*pz3_grid*string(pz3_num)
    
        filename *= "#"*string(px4_low)*"-"*string(px4_up)*px4_grid*string(px4_num)
        filename *= "#"*py4_grid*string(py4_num)
        filename *= "#"*pz4_grid*string(pz4_num)
        
        filename *= ".jld2";

        #Parameters=(name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)
        #filename = DC.FileName(Parameters)

        #Parameters = (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num)

        filename = "CompTestWeighted22.jld2"

        println(filename)

        Output = DC.fload_Matrix(DataDirectory,filename,corrected=true) # UNCOMMENT LATER
        #Output = fload_Matrix(DataDirectory,filename,corrected=true)
        Parameters = Output[1]
        GainMatrix3 = Output[2]
        GainMatrix4 = Output[3]
        LossMatrix1 = Output[4]
        LossMatrix2 = Output[5]

        if name1 == name2 && name3 == name4
            # print conversion statistic
            DC.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            GainCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            println("Scorrected:")
            DC.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(pxr3,pyr3,pxr4,pyr4,pxr1,pyr1,pxr2,pyr2,SMatrix3=matrices[1],TMatrix1=matrices[2])

            Fill_M_Bin!(BigM.M_Bin,interaction,PhaseSpace;SMatrix3=matrices[1],TMatrix1=matrices[2])

        end
    
        if name1 == name2 && name3 != name4

            # print conversion statistic
            DC.DoesConserve(matrices[1],matrices[2],matrices[3],zeros(size(matrices[3])),Parameters)
            #GainCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #DC.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3])
            
            Fill_M_Bin!(BigM.M_Bin,interaction,Lists;SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3])


        end
    
        if name1 != name2 && name3 == name4

            # print conversion statistic
            DC.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],matrices[3],Parameters)
            #GainCorrection!(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)
            #println("Scorrected:")
            #DC.DoesConserve(matrices[1],zeros(size(matrices[1])),matrices[2],zeros(size(matrices[2])),Parameters)

            PhaseSpaceFactors_Binary_Undo!(p3_r,u3_r,p4_r,u4_r,p1_r,u1_r,p2_r,u2_r,SMatrix3=matrices[1],TMatrix1=matrices[2],TMatrix2=matrices[3])

            Fill_M_Bin!(BigM.M_Bin,interaction,Lists;SMatrix3=matrices[1],TMatrix1=matrices[2],TMatrix2=matrices[3])

        end
    
        if name1 != name2 && name3 != name4

            #PhaseSpaceFactors_Binary_Undo!(pxr3,pyr3,pzr3,pxr4,pyr4,pzr4,pxr1,pyr1,pzr1,pxr2,pyr2,pzr2,SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3],TMatrix2=matrices[4])

            # print conversion statistic
            DC.DoesConserve2(Output) # UNCOMMENT LATER
            GainCorrection!(Output)
            println("")
            println("GainCorrected:")
            println("")
            DC.DoesConserve2(Output) # UNCOMMENT LATER
            
            #Fill_M_Bin!(BigM.M_Bin,interaction,PhaseSpace;SMatrix3=matrices[1],SMatrix4=matrices[2],TMatrix1=matrices[3],TMatrix2=matrices[4])
            Fill_M_Bin!(M_Bin,interaction,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

        end

    end # for

end

function LoadMatrices_Emi(M_Emi::Array{Float32,2},DataDirectory::String,PhaseSpace::PhaseSpaceStruct)

    Emi_list = PhaseSpace.Emi_list

    if isempty(Emi_list) # no emission interactions to load
        return
    end
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    fill!(M_Emi,0f0);

    for i in eachindex(Emi_list)

        interaction = Emi_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        name3_loc = findfirst(==(name3),name_list)

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
    
        px1_num::String = string(px_num_list[name1_loc])
        px2_num::String = string(px_num_list[name2_loc])
        px3_num::String = string(px_num_list[name3_loc])
        py1_num::String = string(py_num_list[name1_loc])
        py2_num::String = string(py_num_list[name2_loc])
        py3_num::String = string(py_num_list[name3_loc])
        pz1_num::String = string(pz_num_list[name1_loc])
        pz2_num::String = string(pz_num_list[name2_loc])
        pz3_num::String = string(pz_num_list[name3_loc])
    
        px1_low::String = string(px_low_list[name1_loc])
        px2_low::String = string(px_low_list[name2_loc])
        px3_low::String = string(px_low_list[name3_loc])
        py1_low::String = string(py_low_list[name1_loc])
        py2_low::String = string(py_low_list[name2_loc])
        py3_low::String = string(py_low_list[name3_loc])
        pz1_low::String = string(pz_low_list[name1_loc])
        pz2_low::String = string(pz_low_list[name2_loc])
        pz3_low::String = string(pz_low_list[name3_loc])
        
        px1_up::String = string(px_up_list[name1_loc])
        px2_up::String = string(px_up_list[name2_loc])
        px3_up::String = string(px_up_list[name3_loc])
        py1_up::String = string(py_up_list[name1_loc])
        py2_up::String = string(py_up_list[name2_loc])
        py3_up::String = string(py_up_list[name3_loc])
        pz1_up::String = string(pz_up_list[name1_loc])
        pz2_up::String = string(pz_up_list[name2_loc])
        pz3_up::String = string(pz_up_list[name3_loc])
    
        pxr1::Vector{Float64} = PhaseSpace.Grids.pxr_list[name1_loc]
        pxr2::Vector{Float64} = PhaseSpace.Grids.pxr_list[name2_loc]
        pxr3::Vector{Float64} = PhaseSpace.Grids.pxr_list[name3_loc]
        pyr1::Vector{Float64} = PhaseSpace.Grids.pyr_list[name1_loc]
        pyr2::Vector{Float64} = PhaseSpace.Grids.pyr_list[name2_loc]
        pyr3::Vector{Float64} = PhaseSpace.Grids.pyr_list[name3_loc]
        pzr1::Vector{Float64} = PhaseSpace.Grids.pzr_list[name1_loc]
        pzr2::Vector{Float64} = PhaseSpace.Grids.pzr_list[name2_loc]
        pzr3::Vector{Float64} = PhaseSpace.Grids.pzr_list[name3_loc]

        #filename = "sync"*name2*"#"*px1_low*"-"*px1_up*px1_grid*px1_num*"#"*px3_low*"-"*px3_up*px3_grid*px3_num*"#"*py1_grid*py1_num*"#"*py3_grid*py3_num_st*".jld2";

        #filename = "syncEle#-14.0#4.0#72#-5.0#4.0#72#8#8.jld2";
        #filename = "syncEle#-14.0-7.0l84#0.0-7.0l56#u8#u8.jld2";
        filename = "syncTest3.jld2"

        println(filename)

        #Parameters = DC.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
        #matrix = DC.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
        matrix = fload_Matrix_Sync(DataDirectory,filename)[2] # remove later
            
        # some SMatrix values are greater than float32 precision!
        #PhaseSpaceFactors_Sync_Undo!(matrix,p2_r,u2_r,p1_r,u1_r)
        #PhaseSpaceFactors_Emi_Undo!(matrix,pxr3,pyr3,pxr1,pyr1)

        Fill_M_Emi!(M_Emi,interaction,PhaseSpace;GainMatrix3=matrix)
    
    end # for

end

##### REMOVE LATER #####

function fload_Matrix(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        GainMatrix3 = f["GainMatrix3"];
        GainMatrix4 = f["GainMatrix4"];
        LossMatrix1 = f["LossMatrix1"];
        LossMatrix2 = f["LossMatrix2"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    # old format
    #(name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num) = Parameters

    return (Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

end

function fload_Matrix_Sync(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        GainMatrix3 = f["GainMatrix3"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3)

end

function DoesConserve2(Output::Tuple{Tuple,Array{Float64,9},Array{Float64,9},Array{Float64,6},Array{Float64,6}};Tuple_Output::Bool = false)

    # Output is tuple generated by fload_Matrix

    Parameters = Output[1] # always Output[1]

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    GainMatrix3 = Output[2]
    GainMatrix4 = Output[3]
    LossMatrix1 = Output[4]
    LossMatrix2 = Output[5]

    p1_r = DC.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = DC.deltaVector(p1_r);
    p1_d_full = [p1_d; DC.deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ = DC.deltaEVector(p1_r,mu1);
    E1_Δ_full = [E1_Δ; DC.deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    E1_d_full = E1_Δ_full ./ p1_d_full;

    p2_r = DC.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = DC.deltaVector(p2_r);
    p2_d_full = [p2_d; DC.deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ = DC.deltaEVector(p2_r,mu2);
    E2_Δ_full = [E2_Δ; DC.deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    E2_d_full = E2_Δ_full ./ p2_d_full;

    p3_r = DC.bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = DC.deltaVector(p3_r);
    p3_d_full = [p3_d; DC.deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ = DC.deltaEVector(p3_r,mu3);
    E3_Δ_full = [E3_Δ; DC.deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    E3_d_full = E3_Δ_full ./ p3_d_full

    p4_r = DC.bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = DC.deltaVector(p4_r);
    p4_d_full = [p4_d; DC.deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ = DC.deltaEVector(p4_r,mu4);
    E4_Δ_full = [E4_Δ; DC.deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    E4_d_full = E4_Δ_full ./ p4_d_full

    SsumN3 = 0
    TsumN1 = 0
    SsumE3 = 0
    TsumE1 = 0

    SsumN4 = 0
    TsumN2 = 0
    SsumE4 = 0
    TsumE2 = 0

    NGainMatrix3 = zeros(Float64,size(LossMatrix1))
    NLossMatrix1 = zeros(Float64,size(LossMatrix1))
    EGainMatrix3 = zeros(Float64,size(LossMatrix1))
    ELossMatrix1 = zeros(Float64,size(LossMatrix1))

    NGainMatrix4 = zeros(Float64,size(LossMatrix1))
    NLossMatrix2 = zeros(Float64,size(LossMatrix1))
    EGainMatrix4 = zeros(Float64,size(LossMatrix1))
    ELossMatrix2 = zeros(Float64,size(LossMatrix1))

    for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
        SsumN3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        SsumE3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d_full[p3]
        NGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        EGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d_full[p3]
        end
    end

    for p1 in axes(GainMatrix4, 4), u1 in axes(GainMatrix4,5), h1 in axes(GainMatrix4,6), p2 in axes(GainMatrix4,7), u2 in axes(GainMatrix4,8), h2 in axes(GainMatrix4,9)
        for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
        SsumN4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        SsumE4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d_full[p4]
        NGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        EGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d_full[p4]
        end
    end

    for p1 in axes(LossMatrix1,1), u1 in axes(LossMatrix1, 2), h1 in axes(LossMatrix1,3), p2 in axes(LossMatrix1,4), u2 in axes(LossMatrix1,5), h2 in axes(LossMatrix1,6)
        TsumN1 += LossMatrix1[p1,u1,h1,p2,u2,h2]
        TsumE1 += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d_full[p1]
        TsumN2 += LossMatrix2[p2,u2,h2,p1,u1,h1]
        TsumE2 += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d_full[p2]
        NLossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]
        NLossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]
        ELossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d_full[p1]
        ELossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d_full[p2]
    end

    NErrMatrix1 = (NGainMatrix3 .- NLossMatrix1) ./ NLossMatrix1
    EErrMatrix1 = (EGainMatrix3 .- ELossMatrix1) ./ ELossMatrix1
    meanNErr1 = sum(abs.(NErrMatrix1)) / length(NLossMatrix1)
    meanEErr1 = sum(abs.(EErrMatrix1)) / length(NLossMatrix1)
    stdN1 = sqrt(sum((NErrMatrix1 .- meanNErr1).^2)/length(NLossMatrix1))
    stdE1 = sqrt(sum((EErrMatrix1 .- meanEErr1).^2)/length(NLossMatrix1))

    NErrMatrix2 = (NGainMatrix4 .- NLossMatrix2) ./ NLossMatrix2
    EErrMatrix2 = (EGainMatrix4 .- ELossMatrix2) ./ ELossMatrix2
    meanNErr2 = sum(abs.(NErrMatrix2)) / length(NLossMatrix2)
    meanEErr2 = sum(abs.(EErrMatrix2)) / length(NLossMatrix2)
    stdN2 = sqrt(sum((NErrMatrix2 .- meanNErr2).^2)/length(NLossMatrix2))
    stdE2 = sqrt(sum((EErrMatrix2 .- meanEErr2).^2)/length(NLossMatrix2))

    println("sumSN3 = "*string(SsumN3))
    println("sumSN4 = "*string(SsumN4))
    println("sumTN1 = "*string(TsumN1))    
    println("sumTN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumSN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumTN = "*string(TsumN))

    println("#")

    println("sumSE3 = "*string(SsumE3))
    println("sumSE4 = "*string(SsumE4))
    println("sumTE1 = "*string(TsumE1))  
    println("sumTE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumSE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumTE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string(SsumE-TsumE))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

    println("#")
    println("#")
    println("mean error in N = $meanNErr1")
    println("std of error in  N = $stdN1")
    println("mean error in E = $meanEErr1")
    println("std of error in E = $stdE1")
    println("#")
    println("#")
    println("mean error in N = $meanNErr2")
    println("std of error in  N = $stdN2")
    println("mean error in E = $meanEErr2")
    println("std of error in E = $stdE2")

    if Tuple_Output == true
        return NGainMatrix3, NLossMatrix1, NErrMatrix1,NGainMatrix4, NLossMatrix2, NErrMatrix2
    else
        return nothing
    end
    

end