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

function LoadMatrices_Binary(M_Bin::AbstractMatrix{F},DataDirectory::String,PhaseSpace::PhaseSpaceStruct;mode::ModeType=Ani(),corrected::Bool=true) where F<:AbstractFloat

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

    fill!(M_Bin,zero(F));

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        # Memory optimisation by allowing Ele and Pos populations to be modelled as identical thus only requiring one to be defined
        if !isnothing(findfirst(==("Pos"),[name1,name2,name3,name4])) && isnothing(findfirst(==("Pho"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho") # "Pos" compton scattering
                return # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            elseif (name1,name2,name3,name4) == ("Ele","Pos","Pho","Pho") # "Ele" "Pos" annihilation
                GainScale = 1/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            name1_loc = findfirst(==(name1),name_list)
            name2_loc = findfirst(==(name2),name_list)
            name3_loc = findfirst(==(name3),name_list)
            name4_loc = findfirst(==(name4),name_list)
        end

        # ele pos swap for compton i.e. use ele compton matrices for pos compton
        if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho")
            name1 = "Ele"
            name3 = "Ele"
        end

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        px4_grid::String = px_grid_list[name4_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        py4_grid::String = py_grid_list[name4_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
        pz4_grid::String = pz_grid_list[name4_loc]

        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        px4_num::Int64 = px_num_list[name4_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        py4_num::Int64 = py_num_list[name4_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]
        pz4_num::Int64 = pz_num_list[name4_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]
        px4_low::Float64 = px_low_list[name4_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]
        px4_up::Float64 = px_up_list[name4_loc]

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]
        m4::Float64 = PhaseSpace.Grids.mass_list[name4_loc]

        Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)

        filename::String = BinaryFileName(Parameters)

        println(filename)

        Output = BinaryFileLoad_Matrix(DataDirectory,filename,corrected=corrected)
        Parameters = Output[1]
        GainMatrix3 = Output[2] .* GainScale
        GainMatrix4 = Output[3] .* GainScale
        LossMatrix1 = Output[4]
        LossMatrix2 = Output[5]

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        DoesConserve(Output) # print conversion statistic
        Fill_M_Bin!(M_Bin,name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,mode=mode)

    end # for

end

function LoadMatrices_Emi(M_Emi::AbstractMatrix{F},DataDirectory::String,PhaseSpace::PhaseSpaceStruct;corrected::Bool=true) where F<:AbstractFloat

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

    fill!(M_Emi,zero(F));

    for i in eachindex(Emi_list)

        interaction = Emi_list[i]

        name1::String = interaction.name1
        name2::String = interaction.name2
        name3::String = interaction.name3
        type::String = interaction.EmiName

        # ele pos swap for synchrotron
        if name1 == "Pos" && name2 == "Pos" && name3 == "Pho"
            name1 = "Ele"
            name2 = "Ele"
            name3 = "Pho"
        end

        Ext::Vector{Float64} = interaction.Ext

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        name3_loc = findfirst(==(name3),name_list)

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]

        z1::Float64 = PhaseSpace.Grids.charge_list[name1_loc]
        z2::Float64 = PhaseSpace.Grids.charge_list[name2_loc]
        z3::Float64 = PhaseSpace.Grids.charge_list[name3_loc]

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
    
        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]

        Parameters = (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,Ext)

        filename = DC.EmissionFileName(Parameters)

        println(filename)

        #Parameters = DC.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
        #matrix = DC.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
        matrix = DC.EmissionFileLoad_Matrix(DataDirectory,filename)[2] # remove later

        if corrected
            EmissionCorrection!(PhaseSpace,matrix,Parameters)
        end

        Fill_M_Emi!(M_Emi,interaction,PhaseSpace;GainMatrix3=matrix)
    
    end # for

end