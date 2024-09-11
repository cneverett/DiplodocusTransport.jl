function LoadMatricies(CollisionMatricies::Dict{Vector{String},Tuple},DataDirectory::String,Lists::Tuple)

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    for i in eachindex(interaction_list)

        interaction = interaction_list[i]

        name1_loc = findfirst(==(interaction[1]),name_list)
        name2_loc = findfirst(==(interaction[2]),name_list)
        name3_loc = findfirst(==(interaction[3]),name_list)
        name4_loc = findfirst(==(interaction[4]),name_list)

        name1 = name_list[name1_loc]
        name2 = name_list[name2_loc]
        name3 = name_list[name3_loc]
        name4 = name_list[name4_loc]

        nump1 = string(nump_list[name1_loc])
        nump2 = string(nump_list[name2_loc])
        nump3 = string(nump_list[name3_loc])
        nump4 = string(nump_list[name4_loc])
        numt1 = string(numt_list[name1_loc])
        numt2 = string(numt_list[name2_loc])
        numt3 = string(numt_list[name3_loc])
        numt4 = string(numt_list[name4_loc])

        pl1 = string(pl_list[name1_loc])
        pl2 = string(pl_list[name2_loc])
        pl3 = string(pl_list[name3_loc])
        pl4 = string(pl_list[name4_loc])

        pu1 = string(pu_list[name1_loc])
        pu2 = string(pu_list[name2_loc])
        pu3 = string(pu_list[name3_loc])
        pu4 = string(pu_list[name4_loc])

        filename = name1*name2*name3*name4*"#"*pl1*"#"*pu1*"#"*nump1*"#"*pl2*"#"*pu2*"#"*nump2*"#"*pl3*"#"*pu3*"#"*nump3*"#"*pl4*"#"*pu4*"#"*nump4*"#"*numt1*"#"*numt2*"#"*numt3*"#"*numt4*".jld2"

        Parameters = (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2)

        println(filename)

        matricies = BoltzmannCollisionIntegral.fload_Matrix(DataDirectory,filename)

        # print conversion statistic
        BoltzmannCollisionIntegral.DoesConservere(matricies[1],matricies[2],matricies[3],matricies[4],Parameters)

        if interaction[1] == interaction[2] && interaction[3] == interaction[4]
            SMatrix = SixDtoThreeD(Float32.(matricies[1]))
            TMatrix = FourDtoTwoD(Float32.(matricies[2]))
            CollisionMatricies[interaction] = (SMatrix,TMatrix)
        end
    
        if interaction[1] == interaction[2] && interaction[3] != interaction[4]
            SMatrix3 = SixDtoThreeD(Float32.(matricies[1]))
            SMatrix4 = SixDtoThreeD(Float32.(matricies[2]))
            TMatrix = FourDtoTwoD(Float32.(matricies[3]))
            CollisionMatricies[interaction] = (SMatrix3,SMatrix4,TMatrix)
        end
    
        if interaction[1] != interaction[2] && interaction[3] == interaction[4]
            SMatrix = SixDtoThreeD(Float32.(matricies[1]))
            TMatrix1 = FourDtoTwoD(Float32.(matricies[2]))
            TMatrix2 = FourDtoTwoD(Float32.(matricies[3]))
            CollisionMatricies[interaction] = (SMatrix,TMatrix1,TMatrix2)
        end
    
        if interaction[1] != interaction[2] && interaction[3] != interaction[4]
            SMatrix3 = SixDtoThreeD(Float32.(matricies[1]))
            SMatrix4 = SixDtoThreeD(Float32.(matricies[2]))
            TMatrix1 = FourDtoTwoD(Float32.(matricies[3]))
            TMatrix2 = FourDtoTwoD(Float32.(matricies[4]))
            CollisionMatricies[interaction] = (SMatrix3,SMatrix4,TMatrix1,TMatrix2)
        end

    end

end