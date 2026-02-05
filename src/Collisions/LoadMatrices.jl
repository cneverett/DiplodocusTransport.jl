function LoadMatrices_Binary(M_Bin::AbstractMatrix{F},Binary_list::Vector{BinaryStruct},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,mode::ModeType=Ani(),corrected::Bool=true,Bin_sparse::Bool=false) where F<:AbstractFloat

    Bin_Norm = PhaseSpace.Characteristic.Bin_Norm
    
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
                GainScale = 1.0/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                LossScale = 1.0/4.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                LossScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            LossMatrix1 = 1.0
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
        
        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainScale *= Bin_Norm 
        LossScale *= Bin_Norm

        Parameters = Output[1]
        GainMatrix3 = Output[2] .* GainScale
        GainMatrix4 = Output[3] .* GainScale
        LossMatrix1 = Output[4] .* LossScale
        LossMatrix2 = Output[5] .* LossScale

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        DoesConserve(Output) # print conversion statistic
        Fill_M_Bin!(M_Bin,name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,mode=mode,Bin_sparse=Bin_sparse)

    end # for

end

function LoadMatrices_Emi(M_Emi::AbstractMatrix{F},Emission_list::Vector{EmiStruct},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,Emi_corrected::Bool=true,Emi_sparse::Bool=false) where F<:AbstractFloat

    Emi_Norm = PhaseSpace.Characteristic.Emi_Norm
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Space = PhaseSpace.Space

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

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num

    for i in eachindex(Emi_list)

        interaction = Emi_list[i]

        name1::String = interaction.name1
        name2::String = interaction.name2
        name3::String = interaction.name3
        type::String = interaction.EmiName
        Ext_sampled::Vector{Float64} = interaction.Ext_sampled
        mode::ModeType = interaction.mode
        Force::Bool = interaction.Force
        Domain::Union{Vector{Int64},Nothing} = interaction.Domain

        if !isnothing(findfirst(==("Pos"),[name1,name2,name3])) && isnothing(findfirst(==("Pho"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3) == ("Pos","Pos","Pho") # "Pos" synchrotron 
                return # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            end
        end

        # ele pos swap for synchrotron
        if name1 == "Pos" && name2 == "Pos" && name3 == "Pho"
            name1 = "Ele"
            name2 = "Ele"
            name3 = "Pho"
        end

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

        Parameters = (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,Ext_sampled)

        filename = EmissionFileName(Parameters)

        println(filename)

        #Parameters = DC.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
        #matrix = DC.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
        GainMatrix3_All = EmissionFileLoad_Matrix(DataDirectory,filename)[2] # remove later

        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainMatrix3_All .*= Emi_Norm 

        # apply energy correction
        if Emi_corrected
            for Ext_idx in eachindex(Ext_sampled)
                EmissionCorrection!(PhaseSpace,GainMatrix3_All,Parameters,Ext_idx)
            end
        end

        name_locs = (name1_loc,name2_loc,name3_loc)

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            if isnothing(Domain) || in(off_space,Domain)

                if type=="Sync" && Force 
                    
                    B_field = PhaseSpace.Grids.B_field[x,y,z]
                    Ext_idx = findmin(abs.(Ext_sampled .- B_field))[2]
                    force = SyncRadReact(mode,Ext_sampled[Ext_idx])

                    Fill_I_Emi!(M_Emi,PhaseSpace,force,x,y,z,name1_loc)
                    Fill_J_Emi!(M_Emi,PhaseSpace,force,x,y,z,name1_loc)
                    Fill_K_Emi!(M_Emi,PhaseSpace,force,x,y,z,name1_loc)

                end

                #= Fill_M_Emi! is called for each spatial grid point as the emission correction is dependent on space through the electromagnetic fields =#
                Fill_M_Emi!(M_Emi,PhaseSpace,name_locs,x,y,z;GainMatrix3=GainMatrix3_All,mode=mode)

            else
                continue
            end

        end # space loop

    end # for

end

