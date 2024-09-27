function LoadMatricies_Binary(CollisionMatriciesBinary::Dict{Vector{String},Tuple},DataDirectory::String,Lists::Tuple;mode="AXI")

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list_Binary,interaction_list_Sync) = Lists

    if isempty(interaction_list_Binary) # no binary interactions to load
        return
    end

    for i in eachindex(interaction_list_Binary)

        interaction = interaction_list_Binary[i]

        name1_loc = findfirst(==(interaction[1]),name_list)
        name2_loc = findfirst(==(interaction[2]),name_list)
        name3_loc = findfirst(==(interaction[3]),name_list)
        name4_loc = findfirst(==(interaction[4]),name_list)

        name1 = name_list[name1_loc]
        name2 = name_list[name2_loc]
        name3 = name_list[name3_loc]
        name4 = name_list[name4_loc]

        nump1 = nump_list[name1_loc]
        nump2 = nump_list[name2_loc]
        nump3 = nump_list[name3_loc]
        nump4 = nump_list[name4_loc]
        numt1 = numt_list[name1_loc]
        numt2 = numt_list[name2_loc]
        numt3 = numt_list[name3_loc]
        numt4 = numt_list[name4_loc]

        nump1s = string(nump_list[name1_loc])
        nump2s = string(nump_list[name2_loc])
        nump3s = string(nump_list[name3_loc])
        nump4s = string(nump_list[name4_loc])
        numt1s = string(numt_list[name1_loc])
        numt2s = string(numt_list[name2_loc])
        numt3s = string(numt_list[name3_loc])
        numt4s = string(numt_list[name4_loc])

        pl1 = pl_list[name1_loc]
        pl2 = pl_list[name2_loc]
        pl3 = pl_list[name3_loc]
        pl4 = pl_list[name4_loc]

        pl1s = string(pl_list[name1_loc])
        pl2s = string(pl_list[name2_loc])
        pl3s = string(pl_list[name3_loc])
        pl4s = string(pl_list[name4_loc])

        pu1 = pu_list[name1_loc]
        pu2 = pu_list[name2_loc]
        pu3 = pu_list[name3_loc]
        pu4 = pu_list[name4_loc]

        pu1s = string(pu_list[name1_loc])
        pu2s = string(pu_list[name2_loc])
        pu3s = string(pu_list[name3_loc])
        pu4s = string(pu_list[name4_loc])

        filename = name1*name2*name3*name4*"#"*pl1s*"#"*pu1s*"#"*nump1s*"#"*pl2s*"#"*pu2s*"#"*nump2s*"#"*pl3s*"#"*pu3s*"#"*nump3s*"#"*pl4s*"#"*pu4s*"#"*nump4s*"#"*numt1s*"#"*numt2s*"#"*numt3s*"#"*numt4s*".jld2"

        Parameters = (name1,name2,name3,name4,Float64(pl3),Float64(pu3),nump3,Float64(pl4),Float64(pu4),nump4,Float64(pl1),Float64(pu1),nump1,Float64(pl2),Float64(pu2),nump2,numt3,numt4,numt1,numt2)

        println(filename)

        if mode=="ISO" # ISO
            matricies = BoltzmannCollisionIntegral.fload_Matrix_ISO(DataDirectory,filename)
        elseif mode=="AXI" # AXI
            matricies = BoltzmannCollisionIntegral.fload_Matrix(DataDirectory,filename)
        else
            println("Mode not recognized")
        end

        if interaction[1] == interaction[2] && interaction[3] == interaction[4]
            # print conversion statistic
            BoltzmannCollisionIntegral.DoesConserve(matricies[1],zeros(size(matricies[1])),matricies[2],zeros(size(matricies[2])),Parameters)
            SCorrection!(matricies[1],zeros(size(matricies[1])),matricies[2],zeros(size(matricies[2])),Parameters)
            println("Scorrected:")
            BoltzmannCollisionIntegral.DoesConserve(matricies[1],zeros(size(matricies[1])),matricies[2],zeros(size(matricies[2])),Parameters)

            pr4::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl4,pu4,nump4)
            pr3::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl3,pu3,nump3)
            pr2::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl2,pu2,nump2)
            pr1::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl1,pu1,nump1)
            tr4::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt4)
            tr3::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt3)
            tr2::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt2)
            tr1::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt1)

            #PhaseSpaceFactors_Binary_Undo!(matricies[1],zeros(Float64,size(matricies[1])),matricies[2],pr3,tr3,pr4,tr4,pr1,tr1,pr2,tr2)
            
            SMatrix = SixDtoThreeD(Float32.(matricies[1]))
            TMatrix = FourDtoTwoD(Float32.(matricies[2]))

            if mode=="ISO"
                SMatrix ./= 2*numt1
                TMatrix ./= 2*numt1
            end
            CollisionMatriciesBinary[interaction] = (SMatrix,TMatrix)
        end
    
        if interaction[1] == interaction[2] && interaction[3] != interaction[4]
            SMatrix3 = SixDtoThreeD(Float32.(matricies[1]))
            SMatrix4 = SixDtoThreeD(Float32.(matricies[2]))
            TMatrix = FourDtoTwoD(Float32.(matricies[3]))
            CollisionMatriciesBinary[interaction] = (SMatrix3,SMatrix4,TMatrix)
        end
    
        if interaction[1] != interaction[2] && interaction[3] == interaction[4]
            SMatrix = SixDtoThreeD(Float32.(matricies[1]))
            TMatrix1 = FourDtoTwoD(Float32.(matricies[2]))
            TMatrix2 = FourDtoTwoD(Float32.(matricies[3]))
            CollisionMatriciesBinary[interaction] = (SMatrix,TMatrix1,TMatrix2)
        end
    
        if interaction[1] != interaction[2] && interaction[3] != interaction[4]
            SMatrix3 = SixDtoThreeD(Float32.(matricies[1]))
            SMatrix4 = SixDtoThreeD(Float32.(matricies[2]))
            TMatrix1 = FourDtoTwoD(Float32.(matricies[3]))
            TMatrix2 = FourDtoTwoD(Float32.(matricies[4]))
            CollisionMatriciesBinary[interaction] = (SMatrix3,SMatrix4,TMatrix1,TMatrix2)
        end

    end # for

end


function SCorrection!(SMatrix3::Array{T,6},SMatrix4::Array{T,6},TMatrix1::Array{T,4},TMatrix2::Array{T,4},Parameters) where T

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays

    (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2) = Parameters

    pr3 = BoltzmannCollisionIntegral.prange(p3l,p3u,nump3);
    dp3 = BoltzmannCollisionIntegral.deltaVector(pr3);
    dp3full = [dp3; BoltzmannCollisionIntegral.deltaVector([pr3[end]; 2*pr3[end]])];
    tr3 = BoltzmannCollisionIntegral.trange(numt3);
    dμ3 = BoltzmannCollisionIntegral.deltaVector(tr3);

    pr4 = BoltzmannCollisionIntegral.prange(p4l,p4u,nump4);
    dp4 = BoltzmannCollisionIntegral.deltaVector(pr4);
    dp4full = [dp4; BoltzmannCollisionIntegral.deltaVector([pr4[end]; 2*pr4[end]])];
    tr4 = BoltzmannCollisionIntegral.trange(numt4);
    dμ4 = BoltzmannCollisionIntegral.deltaVector(tr4);

    pr1 = BoltzmannCollisionIntegral.prange(p1l,p1u,nump1);
    dp1 = BoltzmannCollisionIntegral.deltaVector(pr1);
    dp1full = [dp1; BoltzmannCollisionIntegral.deltaVector([pr1[end]; 2*pr1[end]])];
    tr1 = BoltzmannCollisionIntegral.trange(numt1);
    dμ1 = BoltzmannCollisionIntegral.deltaVector(tr1);

    pr2 = BoltzmannCollisionIntegral.prange(p2l,p2u,nump2);
    dp2 = BoltzmannCollisionIntegral.deltaVector(pr2);
    dp2full = [dp2; BoltzmannCollisionIntegral.deltaVector([pr2[end]; 2*pr2[end]])];
    tr2 = BoltzmannCollisionIntegral.trange(numt2);
    dμ2 = BoltzmannCollisionIntegral.deltaVector(tr2);

    for k in axes(SMatrix3,3), l in axes(SMatrix3, 4), m in axes(SMatrix3,5), n in axes(SMatrix3,6)
        SsumN3 = zero(T)
        SsumN4 = zero(T)
        TsumN1 = zero(T)
        TsumN2 = zero(T)
        for i in axes(SMatrix3,1), j in axes(SMatrix3,2) 
        SsumN3 += SMatrix3[i,j,k,l,m,n]*dp3full[i]*dμ3[j]
        end
        for i in axes(SMatrix4,1), j in axes(SMatrix4,2) 
            SsumN4 += SMatrix4[i,j,k,l,m,n]*dp4full[i]*dμ4[j]
            end
        TsumN1 += TMatrix1[k,l,m,n]*dp1full[k]*dμ1[l]
        TsumN2 += TMatrix2[m,n,k,l]*dp2full[m]*dμ2[n]
        SCor = (TsumN1+TsumN2)/(SsumN3+SsumN4)
        @view(SMatrix3[:,:,k,l,m,n]) .*= SCor
        @view(SMatrix4[:,:,k,l,m,n]) .*= SCor
    end


end


function LoadMatricies_Sync(CollisionMatriciesSync::Dict{Vector{String},Array{Float32,2}},DataDirectory::String,Lists::Tuple;mode="AXI")

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list_Binary,interaction_list_Sync) = Lists

    if isempty(interaction_list_Sync) # no sync interactions to load
        return
    end

    for i in eachindex(interaction_list_Sync)

        interaction = interaction_list_Sync[i]

        name1_loc = findfirst(==("Pho"),name_list)
        name2_loc = findfirst(==(interaction[1]),name_list)

        name1 = name_list[name1_loc]
        name2 = name_list[name2_loc]

        nump1 = nump_list[name1_loc]
        nump2 = nump_list[name2_loc]

        numt1 = numt_list[name1_loc]
        numt2 = numt_list[name2_loc]

        nump1s = string(nump_list[name1_loc])
        nump2s = string(nump_list[name2_loc])

        numt1s = string(numt_list[name1_loc])
        numt2s = string(numt_list[name2_loc])

        pl1 = pl_list[name1_loc]
        pl2 = pl_list[name2_loc]

        pl1s = string(pl_list[name1_loc])
        pl2s = string(pl_list[name2_loc])

        pu1 = pu_list[name1_loc]
        pu2 = pu_list[name2_loc]

        pu1s = string(pu_list[name1_loc])
        pu2s = string(pu_list[name2_loc])

        filename = "sync"*name2*"#"*pl1s*"#"*pu1s*"#"*nump1s*"#"*pl2s*"#"*pu2s*"#"*nump2s*"#"*numt1s*"#"*numt2s*".jld2"

        Parameters = (name1,name2,Float64(pl1),Float64(pu1),nump1,Float64(pl2),Float64(pu2),nump2,numt1,numt2)

        println(filename)

        if mode=="ISO" # ISO
            matrix = BoltzmannCollisionIntegral.fload_Matrix_SyncISO(DataDirectory,filename)
        elseif mode=="AXI" # AXI
            matrix = BoltzmannCollisionIntegral.fload_Matrix_Sync(DataDirectory,filename)
        else
            println("Mode not recognized")
        end
            
        # some SMatrix values are greater than float32 precision 
            pr2::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl2,pu2,nump2)
            pr1::Vector{Float64} = BoltzmannCollisionIntegral.prange(pl1,pu1,nump1)
            tr2::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt2)
            tr1::Vector{Float64} = BoltzmannCollisionIntegral.trange(numt1)

        PhaseSpaceFactors_Sync_Undo!(matrix,pr1,tr1,pr2,tr2)
        SMatrix = FourDtoTwoD(Float32.(matrix))

        if mode=="ISO"
            SMatrix ./= 2*numt1
        end

        CollisionMatriciesSync[interaction] = SMatrix

    end # for

end
