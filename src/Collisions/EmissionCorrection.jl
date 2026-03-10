function EmissionCorrection!(PhaseSpace::PhaseSpaceStruct,GainMatrix3::AbstractArray{Float64,6},Parameters,Ext_idx::Int64)

    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,Ext) = Parameters

    Characteristic = PhaseSpace.Characteristic

    if type == "Sync"

        force = SyncRadReact(Ani(),Ext[Ext_idx])
        Momentum = PhaseSpace.Momentum
        Space = PhaseSpace.Space
        scheme = Momentum.scheme
        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        name_list = PhaseSpace.name_list

        Grids = PhaseSpace.Grids
        dE_list = Grids.dE_list

        name1_loc = findfirst(==(name1),name_list)
        name3_loc = findfirst(==(name3),name_list)

        dE1 = dE_list[name1_loc]
        dE3 = dE_list[name3_loc]

        p3r = Grids.pxr_list[name3_loc]
        p1m = Grids.mpx_list[name1_loc]
        

        for px in axes(GainMatrix3, 4), py in axes(GainMatrix3,5), pz in axes(GainMatrix3,6) # loop over p1 states

            # critical frequency
            ω0 = abs((z1*1.6e-19*Ext[Ext_idx]))/(p1m[pz]*9.11e-31)
            pc = log10(1.054e-34*ω0/(9.11e-31*3e8^2)*(p1m[pz])^3)
            pmin = p3r[1]

            if pc < pmin
                println("Critical frequency below minimum momentum, no correction applied")
                continue
            end
        
            GainSumE3 = zero(Float64)
            LossSumE1 = zero(Float64)

            pxp = px+1
            pxm = px-1

            px_num = px1_num

            #= Boundary Conditions:
                Flux on boundaries should be zero i.e. no particles leave/enter from the domain bounds
            =#
            if pxp > px_num
                pxp = px_num
            end
            if pxm < 1
                pxm = 1
            end

            I_plus = zero(Float64)
            I_minus = zero(Float64)

            I_plus += IFluxFunction(force,PhaseSpace,name1_loc,"plus",1,1,1,1,px,py,pz)
            I_minus -= IFluxFunction(force,PhaseSpace,name1_loc,"minus",1,1,1,1,px,py,pz)

            # scheme
            if scheme == "upwind"
                if sign(I_plus) == 1
                    i_plus_right = 0
                    i_plus_left = 1
                else
                    i_plus_right = 1
                    i_plus_left = 0
                end
                if sign(I_minus) == 1
                    i_minus_right = 1
                    i_minus_left = 0
                else
                    i_minus_right = 0
                    i_minus_left = 1
                end
            elseif scheme == "central"
                i_plus_right = 0.5
                i_plus_left = 0.5
                i_minus_right = 0.5
                i_minus_left = 0.5
            else
                error("Unknown scheme")
            end

            
            #=
            ________________________________
            a |  I_m   | I_m+I_p | I_p    |
            __|________|_________|________|_
                    b-1       b        b+1  
            =#

            Mom_norm = MomentumSpaceNorm(Grids,name1_loc,px,py,pz)

            # normalised fluxes
            if px != pxp
                #LossSumE1 += convert(Float64,(I_plus * i_plus_right) / ((pxr[pxp+1]-pxr[pxp])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) * dE1[pxp]
                LossSumE1 += convert(Float64,(I_plus * i_plus_left) / Mom_norm)  * (dE1[px]-dE1[px+1])
            end
            if px != pxm
                LossSumE1 += convert(Float64,(I_minus * i_minus_right) / Mom_norm) * (dE1[px]-dE1[px-1])
                #LossSumE1 += convert(Float64,(I_minus * i_minus_left) / ((pxr[pxm+1]-pxr[pxm])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) * dE1[pxm]
            end

            # calculate total rate of energy gain from p1 state
            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                GainSumE3 += GainMatrix3[p3,u3,h3,px,py,pz] * dE3[p3]
            end

            vol = VolFunction(PhaseSpace,1,1,1,1)

            GainSumE3 *= vol

            if GainSumE3 != 0e0
                Correction = (LossSumE1)/GainSumE3
                if Correction < 0.0 
                    println("Negative correction factor, check flux calculations, setting correction to zero")
                    Correction = 0.0
                end
                println("$Correction")
                @view(GainMatrix3[:,:,:,px,py,pz]) .= Correction * @view(GainMatrix3[:,:,:,px,py,pz])
            end

        end # loop over p1 states

    else
        return nothing

    end

end