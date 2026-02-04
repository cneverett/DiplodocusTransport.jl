"""
    BuildFluxMatrices(PhaseSpace;debug_mode=false,Precision=Float32)

Function that builds the flux matrices associated with coordinate forces and regular forces. First space is allocated for the arrays, then the fluxes are filled and finally the flux matrices are returned as an immutable `FluxMatricesStruct`.

# Optional arguments
- `debug_mode::Bool=false`: If true, all flux matrices are built and returned. If false, only the essential flux matrices are built and returned to save memory.
- `Precision::DataType=Float32`: The data type for the elements of the flux matrices and volume elements. Defines the machine error of the simulation. Can be either `Float32` or `Float64`.
"""
function BuildFluxMatrices(PhaseSpace::PhaseSpaceStruct;debug_mode::Bool=false,Precision::DataType=Float32)

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    (Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux) = Allocate_Flux(PhaseSpace,Precision,debug_mode)

    if debug_mode

        # Fill flux matrices
        Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)
        Fill_B_Flux!(B_Flux,PhaseSpace)
        Fill_C_Flux!(C_Flux,PhaseSpace)
        Fill_D_Flux!(D_Flux,PhaseSpace)
        Fill_I_Flux!(I_Flux,PhaseSpace) 
        Fill_J_Flux!(J_Flux,PhaseSpace)
        Fill_K_Flux!(K_Flux,PhaseSpace) 
        Fill_Vol!(Vol,PhaseSpace)
        @. F_Flux = I_Flux + J_Flux + K_Flux + B_Flux + C_Flux + D_Flux  

    else

        (Ap_Flux,Am_Flux,F_Flux,Vol) = Allocate_Flux(PhaseSpace,Precision,debug_mode)
        # Fill flux matrices
        Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)
        Fill_B_Flux!(F_Flux,PhaseSpace)
        Fill_C_Flux!(F_Flux,PhaseSpace)
        Fill_D_Flux!(F_Flux,PhaseSpace)
        Fill_I_Flux!(F_Flux,PhaseSpace) 
        Fill_J_Flux!(F_Flux,PhaseSpace)
        Fill_K_Flux!(F_Flux,PhaseSpace) 
        Fill_Vol!(Vol,PhaseSpace)

        if isdiag(Ap_Flux)
            Ap_Flux = convert(Vector{Precision},diag(Ap_Flux))::Vector{Precision}
        else
            error("Aplus flux is not diagonal, which should be the case for stationary spacetimes")
        end
        if isdiag(Am_Flux)
            Am_Flux = convert(Vector{Precision},diag(Am_Flux))::Vector{Precision}
        else
            error("Aminus flux is not diagonal, which should be the case for stationary spacetimes")
        end

    end

    return FluxMatricesStruct{Precision}(Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

end

"""
    Allocate_Flux(PhaseSpace,Precision,debug_mode)

Allocates arrays for fluxes and volume elements.
"""
function Allocate_Flux(PhaseSpace::PhaseSpaceStruct,Precision::DataType,debug_mode::Bool)

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    n_space = x_num*y_num*z_num

    n = n_momentum*n_space

    # boundary terms included in arrays
    # fluxes allocated with all zeros
    # time fluxes
    Ap_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n) 
    Am_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
    # volume element 
    Vol::Vector{Precision} = zeros(Precision,n_space)
    # sum of space and momentum fluxes 
    F_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
    if debug_mode
        # space fluxes
        B_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        C_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        D_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        # momentum fluxes
        I_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        J_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        K_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
    else
        # space fluxes
        B_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
        C_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
        D_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
        # momentum fluxes
        I_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
        J_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
        K_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64} 
    end

    return (Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

end

"""
    Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)

Generates `Ap_Flux` and `Am_Flux` terms.
"""
function Fill_A_Flux!(Ap_Flux::SparseMatrixCSC{T,Int64},Am_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)


        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a

                # integration sign introduced here
                A_plus = AFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                A_minus = -AFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                # normalisation
                norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # fill
                Ap_Flux[a,b] += convert(T,A_plus / norm)
                Am_Flux[a,b] += convert(T,A_minus / norm)

            end
        end
    end
end

"""
    Fill_Vol!(Vol,PhaseSpace)

Generates volume elements for each space sub-domain.
"""
function Fill_Vol!(Vol::Vector{T},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        space = (x-1)*y_num*z_num+(y-1)*z_num+z

        Vol[space] = convert(T,VolFunction(PhaseSpace,1,x,y,z))

    end

end

"""
    Fill_I_Flux!(I_Flux,PhaseSpace)

Generates `I_Flux` terms.
"""
function Fill_I_Flux!(I_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

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

    Forces = PhaseSpace.Forces
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

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

            for f in 1:length(Forces)
                # integration sign introduced here
                I_plus = IFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"plus",1,x,y,z,px,py,pz)
                I_minus = -IFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"minus",1,x,y,z,px,py,pz)
        
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

            end # Forces loop
        end # end coordinates loop
    end
end

"""
    Fill_J_Flux!(J_Flux,PhaseSpace)

Generates `J_Flux` terms.
"""
function Fill_J_Flux!(J_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

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

    Forces = PhaseSpace.Forces
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

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

            for f in 1:length(Forces)
                # integration sign introduced here
                J_plus = JFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"plus",1,x,y,z,px,py,pz)
                J_minus = -JFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"minus",1,x,y,z,px,py,pz)

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

            end # Force loop
        end
    end
end

"""
    Fill_K_Flux!(K_Flux,PhaseSpace)

Generates `K_Flux` terms.
"""
function Fill_K_Flux!(K_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

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

    Forces = PhaseSpace.Forces
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

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

            for f in 1:length(Forces)
                # integration sign introduced here
                K_plus = KFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"plus",1,x,y,z,px,py,pz)
                K_minus = -KFluxFunction(Forces[f],space_coords,momentum_coords,Grids,Characteristic,name,"minus",1,x,y,z,px,py,pz)

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
                
            end # force loop
        end
    end
end


"""
    Fill_B_Flux!(B_Flux,PhaseSpace)

Generates `B_Flux` terms.
"""
function Fill_B_Flux!(B_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme

    BCp = space_coords.xp_BC
    BCm = space_coords.xm_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
            b = a 
            xp = x+1
            xm = x-1
            if xp > x_num # right boundary
                if BCp isa Closed || BCp isa Open
                    xp = x_num
                elseif BCp isa Periodic
                    xp = 1
                end
            end
            if xm < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    xm = 1
                elseif BCm isa Periodic
                    xm = x_num
                end
            end

            bp = GlobalIndices_To_StateIndex(xp,y,z,px,py,pz,name,PhaseSpace)
            bm = GlobalIndices_To_StateIndex(xm,y,z,px,py,pz,name,PhaseSpace)

            # integration sign introduced here
            B_plus = BFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            B_minus = -BFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(B_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(B_minus) == 1
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
            a |  B_m   | B_m+B_p | B_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                B_Flux[a,bp] += convert(T,(B_plus * h_plus_right) / Mon_Norm) 
                B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
            elseif BCp isa Open # b=bp
                B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
            end
            if b != bm
                B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
                B_Flux[a,bm] += convert(T,(B_minus * h_minus_left) / Mon_Norm) 
            elseif BCm isa Open # b=bm
                B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    Fill_C_Flux!(C_Flux,PhaseSpace)

Generates `C_Flux` terms.
"""
function Fill_C_Flux!(C_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = space_coords.yp_BC
    BCm = space_coords.ym_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
            b = a 
            yp = y+1
            ym = y-1
            if yp > y_num # right boundary
                if BCp isa Closed || BCp isa Open
                    yp = y_num
                elseif BCp isa Periodic
                    yp = 1
                end
            end
            if ym < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    ym = 1
                elseif BCm isa Periodic
                    ym = y_num
                end
            end

            bp = GlobalIndices_To_StateIndex(x,yp,z,px,py,pz,name,PhaseSpace)
            bm = GlobalIndices_To_StateIndex(x,ym,z,px,py,pz,name,PhaseSpace)

            # integration sign introduced here
            C_plus = CFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            C_minus = -CFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(C_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(C_minus) == 1
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
            a |  C_m   | C_m+C_p | C_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                C_Flux[a,bp] += convert(T,(C_plus * h_plus_right) / Mom_Norm) 
                C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
            elseif BCp isa Open # b=bp
                C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
            end
            if b != bm
                C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
                C_Flux[a,bm] += convert(T,(C_minus * h_minus_left) / Mom_Norm) 
            elseif BCm isa Open # b=bm
                C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    Fill_D_Flux!(D_Flux,PhaseSpace)

Generates `D_Flux` terms.
"""
function Fill_D_Flux!(D_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme

    BCp = space_coords.zp_BC
    BCm = space_coords.zm_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
            b = a 
            zp = z+1
            zm = z-1
            if zp > z_num # right boundary
                if BCp isa Closed || BCp isa Open || BCp isa Reflective
                    zp = z_num
                elseif BCp isa Periodic
                    zp = 1
                end
            end
            if zm < 1 # left boundary
                if BCm isa Closed || BCm isa Open || BCm isa Reflective
                    zm = 1
                elseif BCm isa Periodic
                    zm = z_num
                end
            end

            if BCp isa Reflective && (z+1 > z_num)
                @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                # mirror u momentum at boundary
                py_ref = py_num - py + 1
                bp = GlobalIndices_To_StateIndex(x,y,zp,px,py_ref,pz,name,PhaseSpace)
            else
                bp = GlobalIndices_To_StateIndex(x,y,zp,px,py,pz,name,PhaseSpace)
            end
            
            if BCm isa Reflective && (z-1 < 1)
                @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                # mirror u momentum at boundary
                py_ref = py_num - py + 1
                bm = GlobalIndices_To_StateIndex(x,y,zm,px,py_ref,pz,name,PhaseSpace)
            else
                bm = GlobalIndices_To_StateIndex(x,y,zm,px,py,pz,name,PhaseSpace)
            end

            # integration sign introduced here
            D_plus = DFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            D_minus = -DFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(D_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(D_minus) == 1
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
            a |  D_m   | D_m+D_p | D_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                D_Flux[a,bp] += convert(T,(D_plus * h_plus_right) / Mon_Norm) 
                D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
            elseif BCp isa Open # b=bp
                D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
            end
            if b != bm
                D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
                D_Flux[a,bm] += convert(T,(D_minus * h_minus_left) / Mon_Norm) 
            elseif BCm isa Open # b=bm
                D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    MomentumSpaceNorm(Grids,species_index,px_idx,py_idx,pz_idx)

Returns the normalization factor for a species distribution function in given momentum space sub-domain.
"""
function MomentumSpaceNorm(Grids::GridsStruct,species_index,px_idx,py_idx,pz_idx)

    dpx = Grids.dpx_list[species_index][px_idx]
    dpy = Grids.dpy_list[species_index][py_idx]
    dpz = Grids.dpz_list[species_index][pz_idx]

    return dpx*dpy*dpz

end