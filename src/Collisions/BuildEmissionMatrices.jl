"""
    BuildEmissionMatrices(PhaseSpace,DataDirectory;loading_check)

Function that builds the emission matrices associated with emissive interactions and their corresponding reaction forces if applicable. If there are such interactions, first space is allocated for the arrays, then data is loaded into these arrays from the desired `DataDirectory` location and finally the big matrices are returned as an immutable `EmissionMatricesStruct`.
"""
function BuildEmissionMatrices(PhaseSpace::PhaseSpaceStruct,Emission_list::Vector{EmiStruct},DataDirectory::String;loading_check::Bool=false,Emi_corrected::Bool=true,Emi_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"
    
    if Emi_sparse == true && Emi_corrected == false
        error("Emi_sparse=true requires Emi_corrected=true")
    end

    if isempty(Emission_list)
        if Emi_sparse
            M_Emi = spzeros(Precision,0,0)
        else
            M_Emi = zeros(Precision,0,0)
        end 
    else
        M_Emi = Allocate_M_Emi(PhaseSpace,loading_check,Precision,Emi_sparse)
        LoadMatrices_Emi(M_Emi,Emission_list,DataDirectory,PhaseSpace,Emi_corrected,Emi_sparse)
    end

    return EmissionMatricesStruct{Precision}(M_Emi,Emission_list)

end


"""
    Allocate_M_Emi(Momentum::MomentumStruct)

Allocates a big matrix `M_Emi` which stores the interaction rates for emission interactions between all particles in the simulation.
"""
function Allocate_M_Emi(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType,Emi_sparse::Bool)

    Momentum = PhaseSpace.Momentum
    Space = PhaseSpace.Space

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = Space.x_num*Space.y_num*Space.z_num

    n = n_momentum*n_space

    size = (n)^2*sizeof(Precision)

    if loading_check

        if size > 1e9
            println("M_Emi will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Emi will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Emi will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Emi will be approx. $size bytes in memory if dense")
        end

        #println("Do you want to proceed? (y/n)")
        #proceed = readline()
        #proceed = "y" ## FIX LATER
        #if proceed != "y"
        #    error("Error: User aborted emission interaction matrix loading")
        #end

    else

        println("Building emission interaction matrix")
        if size > 1e9
            println("M_Emi is approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Emi is approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Emi is approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Emi is approx. $size bytes in memory if dense")
        end

    end

    if Emi_sparse
        M_Emi = spzeros(Precision,n,n)
    else 
        M_Emi = zeros(Precision,n,n)
    end

    return M_Emi

end

function Fill_M_Emi!(M_Emi::AbstractMatrix{F},PhaseSpace::PhaseSpaceStruct,names::Tuple{Int64,Int64,Int64},x::Int64,y::Int64,z::Int64;GainMatrix2=nothing,GainMatrix3=nothing,LossMatrix1=nothing,mode::ModeType=Ani()) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids

    (name1,name2,name3) = names

    # interaction is name1 -> name2 + name3
    # absorption of name1 and emission of name2 not implemented with name1 == name2

    #GainMatrix_to_M_Emi!(M_Emi,PhaseSpace,GainMatrix2,name2,name1,x,y,z,mode)
    GainMatrix_to_M_Emi!(M_Emi,PhaseSpace,GainMatrix3,name3,name1,x,y,z,mode)
    #LossMatrix_to_M_Emi!(M_Emi,PhaseSpace,LossMatrix1,name1,x,y,z,mode)

    #GainMatrix2 = nothing
    GainMatrix3 = nothing
    #LossMatrix1 = nothing

    GC.gc()

end


function GainMatrix_to_M_Emi!(M_Emi::AbstractMatrix{F},PhaseSpace::PhaseSpaceStruct,GainMatrix::Array{Float64,6},name2::Int64,name1::Int64,x::Int64,y::Int64,z::Int64,mode::ModeType) where F<:Union{Float32,Float64}

    vol = VolFunction(PhaseSpace,x,y,z)

    dpy1 = PhaseSpace.Grids.dpy_list[name1]
    dpy2 = PhaseSpace.Grids.dpy_list[name2]
    dpz1 = PhaseSpace.Grids.dpz_list[name1]
    dpz2 = PhaseSpace.Grids.dpz_list[name2]
    
    for px1 in axes(GainMatrix,4), px2 in axes(GainMatrix,1)

        val = 0.0
        w = 1.0

        if mode isa Iso

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in axes(GainMatrix,5), pz1 in axes(GainMatrix,6), py2 in axes(GainMatrix,2), pz2 in axes(GainMatrix,3)
                val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpy2[py2] * dpz2[pz2] * dpy1[py1] * dpz1[pz1]
            end

            w /= (sum(dpy1) * sum(dpz1)) * (sum(dpy2) * sum(dpz2))

        end

        for py1 in axes(GainMatrix,5), py2 in axes(GainMatrix,2)
            
            if mode isa Axi

                # average over incoming and outgoing phi angles (pz)
                for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3)
                    val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpz2[pz2] * dpz1[pz1]
                end

                w /= sum(dpz1) * sum(dpz2)

            end

            for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3)

                if mode isa Ani 

                    val += GainMatrix[px2,py2,pz2,px1,py1,pz1]

                end

                a = GlobalIndices_To_StateIndex(x,y,z,px2,py2,pz2,name2,PhaseSpace)
                b = GlobalIndices_To_StateIndex(x,y,z,px1,py1,pz1,name1,PhaseSpace)


                M_Emi[a,b] += convert(F,val*w*vol) 

            end

        end

    end

end

function LossMatrix_to_M_Emi!(M_Emi::Matrix{F},PhaseSpace::PhaseSpaceStruct,LossMatrix::Array{Float64,2},name1::Int64,x::Int64,y::Int64,z::Int64,mode::ModeType) where F<:Union{Float32,Float64}

    vol = VolFunction(PhaseSpace,x,y,z)

    dpy1 = PhaseSpace.Grids.dpy_list[name1]
    dpz1 = PhaseSpace.Grids.dpz_list[name1]

    for px1 in axes(LossMatrix,1)
        
        val = 0.0 
        w = 1.0 

        if mode isa Iso

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in axes(LossMatrix,2), pz1 in axes(LossMatrix,3)
                val += LossMatrix[px1,py1,pz1] * dpy1[py1] * dpz1[pz1]
            end

            w /= (sum(dpy1) * sum(dpz1))

        end

        for py1 in axes(LossMatrix,2)

            if mode isa Axi

                # average over incoming phi angles (pz)
                for pz1 in axes(LossMatrix,3)
                    val += LossMatrix[px1,py1,pz1] * dpz1[pz1]
                end

                w /= sum(dpz1)

            end

            for pz1 in axes(LossMatrix,3)

                if mode isa Ani 

                    val += LossMatrix[px1,py1,pz1]

                end

                a = GlobalIndices_To_StateIndex(x,y,z,px1,py1,pz1,name1,PhaseSpace)
                b = a

                M_Emi[a,b] -= convert(F,val*w*vol) 

            end

        end
        
    end

end


"""
    Fill_I_Emi!(I_Flux,PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx)

Generates `I_Flux` terms in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_I_Emi!(I_Flux::AbstractMatrix{T},PhaseSpace::PhaseSpaceStruct,Force::ForceType,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme

    BCp = momentum_coords.xp_BC
    BCm = momentum_coords.xm_BC 

    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert BCp isa Closed || BCm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx

    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
        b = a 
        pxp = px+1
        pxm = px-1
        if pxp > px_num # right boundary
            if BCp isa Closed || BCp isa Open
                pxp = px_num
            elseif BCp isa Periodic
                pxp = 1
            end
        end
        if pxm < 1 # left boundary
            if BCm isa Closed || BCm isa Open
                pxm = 1
            elseif BCm isa Periodic
                pxm = px_num
            end
        end

        bp = GlobalIndices_To_StateIndex(x,y,z,pxp,py,pz,name,PhaseSpace)
        bm = GlobalIndices_To_StateIndex(x,y,z,pxm,py,pz,name,PhaseSpace)

        # integration sign introduced here
        I_plus = IFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        I_minus = -IFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(I_plus) == 1
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(I_minus) == 1
                h_minus_right = 1
                h_minus_left = 0
            else
                h_minus_right = 0
                h_minus_left = 1
            end
        elseif scheme == "central"
            h_plus_right = 0.5
            h_plus_left = 0.5
            h_minus_right = 0.5
            h_minus_left = 0.5
        else
            error("Unknown scheme")
        end

        
        #=
        ________________________________
        a |  I_m   | I_m+I_p | I_p    |
        __|________|_________|________|_
                b-1       b        b+1  
        =#

        ## note on boundaries:
        # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
        # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

        Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
        Mom_Normp = MomentumSpaceNorm(Grids,name,pxp,py,pz)
        Mom_Normm = MomentumSpaceNorm(Grids,name,pxm,py,pz)

        # normalised fluxes
        if b != bp
            I_Flux[a,bp] += convert(T,(I_plus * h_plus_right) / Mom_Normp)
            I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
        elseif BCp isa Open # b=bp
            I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
        end
        if b != bm
            I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
            I_Flux[a,bm] += convert(T,(I_minus * h_minus_left) / Mom_Normm) 
        elseif BCm isa Open # b=bm
            I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
        end

    end # end coordinates loop

end

"""
    Fill_J_Emi!(J_Flux,PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx)

Generates `J_Flux` term in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_J_Emi!(J_Flux::AbstractMatrix{T},PhaseSpace::PhaseSpaceStruct,Force::ForceType,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.yp_BC
    BCm = momentum_coords.ym_BC 
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert (BCp isa Closed) && (BCm isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx
    
    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
        b = a 
        pyp = py+1
        pym = py-1
        if pyp > py_num # right boundary
            if BCp isa Closed || BCp isa Open
                pyp = py_num
            elseif BCp isa Periodic
                pyp = 1
            end
        end
        if pym < 1 # left boundary
            if BCm isa Closed || BCm isa Open
                pym = 1
            elseif BCm isa Periodic
                pym = py_num
            end
        end

        bp = GlobalIndices_To_StateIndex(x,y,z,px,pyp,pz,name,PhaseSpace)
        bm = GlobalIndices_To_StateIndex(x,y,z,px,pym,pz,name,PhaseSpace)

        # integration sign introduced here
        J_plus = JFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        J_minus = -JFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(J_plus) == 1
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(J_minus) == 1
                h_minus_right = 1
                h_minus_left = 0
            else
                h_minus_right = 0
                h_minus_left = 1
            end
        elseif scheme == "central"
            h_plus_right = 0.5
            h_plus_left = 0.5
            h_minus_right = 0.5
            h_minus_left = 0.5
        else
            error("Unknown scheme")
        end

        #=
        ________________________________
        a |  J_m   | J_m+J_p | J_p    |
        __|________|_________|________|_
            b-1        b        b+1  
        =#

        ## note on boundaries:
        # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
        # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

        Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
        Mom_Normp = MomentumSpaceNorm(Grids,name,px,pyp,pz)
        Mom_Normm = MomentumSpaceNorm(Grids,name,px,pym,pz)

        # normalised fluxes
        if b != bp
            J_Flux[a,bp] += convert(T,(J_plus * h_plus_right) / Mom_Normp)
            J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
        elseif BCp isa Open # b=bp
            J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
        end
        if b != bm
            J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
            J_Flux[a,bm] += convert(T,(J_minus * h_minus_left) / Mom_Normm)
        elseif BCm isa Open # b=bm
            J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
        end

    end
end

"""
    Fill_K_Flux!(K_Flux,PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx)

Generates `K_Flux` terms in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_K_Flux!(K_Flux::AbstractMatrix{T},PhaseSpace::PhaseSpaceStruct,Force::ForceType,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.zp_BC
    BCm = momentum_coords.zm_BC 
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    @assert (BCp isa Periodic) && (BCm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx

    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
        b = a 
        pzp = pz+1
        pzm = pz-1
        if pzp > pz_num # right boundary
            if BCp isa Closed || BCp isa Open
                pzp = pz_num
            elseif BCp isa Periodic
                pzp = 1
            end
        end
        if pzm < 1 # left boundary
            if BCm isa Closed || BCm isa Open
                pzm = 1
            elseif BCm isa Periodic
                pzm = pz_num
            end
        end

        bp = GlobalIndices_To_StateIndex(x,y,z,px,py,pzp,name,PhaseSpace)
        bm = GlobalIndices_To_StateIndex(x,y,z,px,py,pzm,name,PhaseSpace)

        # integration sign introduced here
        K_plus = KFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        K_minus = -KFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(K_plus) == 1
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(K_minus) == 1
                h_minus_right = 1
                h_minus_left = 0
            else
                h_minus_right = 0
                h_minus_left = 1
            end
        elseif scheme == "central"
            h_plus_right = 1/2
            h_plus_left = 1/2
            h_minus_right = 1/2
            h_minus_left = 1/2
        else
            error("Unknown scheme")
        end

        #=
        ________________________________
        a |  K_m   | K_m+K_p | K_p    |
        __|________|_________|________|_
            b-1        b        b+1  
        =#

        ## note on boundaries:
        # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
        # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

        Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
        Mom_Normp = MomentumSpaceNorm(Grids,name,px,py,pzp)
        Mom_Normm = MomentumSpaceNorm(Grids,name,px,py,pzm)

        if b != bp
            K_Flux[a,bp] += convert(T,(K_plus * h_plus_right) / Mom_Normp)
            K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
        elseif BCp isa Open # b=bp
            K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
        end
        if b != bm
            K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
            K_Flux[a,bm] += convert(T,(K_minus * h_minus_left) / Mom_Normm)
        elseif BCm isa Open # b=bm
            K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
        end
                
    end
end
