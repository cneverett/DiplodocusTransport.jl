"""
    BuildFluxMatrices(PhaseSpace;debug_mode=false,Precision=Float32)

Function that builds the flux matrices associated with coordinate forces and regular forces. First space is allocated for the arrays, then the fluxes are filled and finally the flux matrices are returned as an immutable `FluxMatricesStruct`.

# Optional arguments
- `debug_mode::Bool=false`: If true, all flux matrices are built and returned. If false, only the essential flux matrices are built and returned to save memory.
- `Precision::DataType=Float32`: The data type for the elements of the flux matrices and volume elements. Defines the machine error of the simulation. Can be either `Float32` or `Float64`.
"""
function BuildFluxMatrices(PhaseSpace::PhaseSpaceStruct;debug_mode::Bool=false,Precision::DataType=Float32)

    if Precision != Float32 && Precision != Float64
        error("Precision must be either Float32 or Float64")
    end

    (Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux) = Allocate_Flux(PhaseSpace,Precision,debug_mode)

    if debug_mode

        # Fill flux matrices
        Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)
        # to be implemented
        #Fill_B_Flux!(B_Flux,Lists,PhaseSpace)
        #Fill_C_Flux!(C_Flux,Lists,PhaseSpace)
        #Fill_D_Flux!(D_Flux,Lists,PhaseSpace)
        Fill_I_Flux!(I_Flux,PhaseSpace) 
        Fill_J_Flux!(J_Flux,PhaseSpace)
        Fill_K_Flux!(K_Flux,PhaseSpace) 
        Fill_Vol!(Vol,PhaseSpace)
        @. F_Flux = I_Flux + J_Flux + K_Flux #+ B_Flux + C_Flux + D_Flux  

    else

        (Ap_Flux,Am_Flux,F_Flux,Vol) = Allocate_Flux(PhaseSpace,Precision,debug_mode)
        # Fill flux matrices
        Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)
        # to be implemented
        #Fill_B_Flux!(F_Flux,Lists,PhaseSpace)
        #Fill_C_Flux!(F_Flux,Lists,PhaseSpace)
        #Fill_D_Flux!(F_Flux,Lists,PhaseSpace)
        Fill_I_Flux!(F_Flux,PhaseSpace) 
        Fill_J_Flux!(F_Flux,PhaseSpace)
        Fill_K_Flux!(F_Flux,PhaseSpace) 
        Fill_Vol!(Vol,PhaseSpace)

        if isdiag(Ap_Flux)
            Ap_Flux = diag(Ap_Flux)::Vector{Precision}
        else
            error("Aplus flux is not diagonal, which should be the case for stationary spacetimes")
        end
        if isdiag(Am_Flux)
            Am_Flux = diag(Am_Flux)::Vector{Precision}
        else
            error("Aminus flux is not diagonal, which should be the case for stationary spacetimes")
        end

    end

    return FluxMatricesStruct{Precision}(Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

end

"""
    Allocate_Flux(PhaseSpace,MatrixType,VectorType)

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
        B_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0) # change 0 to n when implemented
        C_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0) # change 0 to n when implemented
        D_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0) # change 0 to n when implemented
        # momentum fluxes
        I_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
        J_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
        K_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
    else
        # space fluxes
        B_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0) 
        C_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0)
        D_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0)
        # momentum fluxes
        I_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0)
        J_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0)
        K_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,0,0) 
    end

    return (Ap_Flux,Am_Flux,F_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

end

function Fill_A_Flux!(Ap_Flux::SparseMatrixCSC{T,Int64},Am_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    offset = PhaseSpace.Grids.momentum_species_offset

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = x_num*y_num*z_num

    tr = Grids.tr # A_flux uses only first time step!
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr
    
    for name in 1:length(name_list)
        off_name = offset[name]
        pxr = Grids.pxr_list[name]
        pyr = Grids.pyr_list[name]
        pzr = Grids.pzr_list[name]

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = (pz-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                b = a

                # integration sign introduced here
                A_plus = AFluxFunction(space_coords,momentum_coords,tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                A_minus = -AFluxFunction(space_coords,momentum_coords,tr[1],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])

                # normalisation
                norm = (pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])

                # fill
                Ap_Flux[a,b] += convert(T,A_plus / norm)
                Am_Flux[a,b] += convert(T,A_minus / norm)

            end
        end
    end
end

function Fill_Vol!(Vol::Vector{T},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num

    tr = Grids.tr # time step assumed to be constant!
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        space = (x-1)*y_num*z_num+(y-1)*z_num+z

        Vol[space] = convert(T,VolFunction(space_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1]))

    end

end

function Fill_I_Flux!(I_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.xp_BC
    BCm = momentum_coords.xm_BC 

    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    if typeof(BCp) != Closed || typeof(BCm) != Closed
        error("I flux boundaries incorrectly defined, i.e. not closed")
    end

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = x_num*y_num*z_num

    tr = Grids.tr # time step assumed to be constant!
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr
    
    for name in 1:length(name_list)
        off_name = offset[name]
        pxr = Grids.pxr_list[name]
        pyr = Grids.pyr_list[name]
        pzr = Grids.pzr_list[name]

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,species_index,PhaseSpace)
            b = a 
            pxp = px+1
            pxm = px-1
            if pxp > px_num # right boundary
                if typeof(BCp) == Closed || typeof(BCp) == Open
                    pxp = px_num
                elseif typeof(BCp) == Periodic
                    pxp = 1
                end
            end
            if pxm < 1 # left boundary
                if typeof(BCm) == Closed || typeof(BCm) == Open
                    pxm = 1
                elseif typeof(BCm) == Periodic
                    pxm = px_num
                end
            end

            bp = GlobalIndices_To_StateIndex(x,y,z,pxp,py,pz,species_index,PhaseSpace)
            bm = GlobalIndices_To_StateIndex(x,y,z,pxm,py,pz,species_index,PhaseSpace)

            I_plus = zero(type)
            I_minus = zero(type)

            for f in 1:length(Forces)
                # integration sign introduced here
                I_plus = IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                I_minus = -IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
        
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

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                # normalised fluxes
                if b != bp
                    I_Flux[a,bp] += convert(T,(I_plus * i_plus_right) / ((pxr[pxp+1]-pxr[pxp])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                    I_Flux[a,b] += convert(T,(I_plus * i_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                elseif typeof(BCp) == Open # b=bp
                    I_Flux[a,b] += convert(T,(I_plus * i_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                end
                if b != bm
                    I_Flux[a,b] += convert(T,(I_minus * i_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                    I_Flux[a,bm] += convert(T,(I_minus * i_minus_left) / ((pxr[pxm+1]-pxr[pxm])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                elseif typeof(BCm) == Open # b=bm
                    I_Flux[a,b] += convert(T,(I_minus * i_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) 
                end

            end # Forces loop
        end # end coordinates loop
    end
end

function Fill_J_Flux!(J_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.yp_BC
    BCm = momentum_coords.ym_BC 
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    if typeof(BCp) != Closed || typeof(BCm) != Closed
        error("J flux boundaries incorrectly defined, i.e. not closed")
    end

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = x_num*y_num*z_num

    tr = Grids.tr # time step assumed to be constant!
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr
    
    for name in 1:length(name_list)

        off_name = offset[name]
        pxr = Grids.pxr_list[name]
        pyr = Grids.pyr_list[name]
        pzr = Grids.pzr_list[name]

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,species_index,PhaseSpace)
            b = a 
            pyp = py+1
            pym = py-1
            if pyp > py_num # right boundary
                if typeof(BCp) == Closed || typeof(BCp) == Open
                    pyp = py_num
                elseif typeof(BCp) == Periodic
                    pyp = 1
                end
            end
            if pym < 1 # left boundary
                if typeof(BCm) == Closed || typeof(BCm) == Open
                    pym = 1
                elseif typeof(BCm) == Periodic
                    pym = py_num
                end
            end

            bp = GlobalIndices_To_StateIndex(x,y,z,px,pyp,pz,species_index,PhaseSpace)
            bm = GlobalIndices_To_StateIndex(x,y,z,px,pym,pz,species_index,PhaseSpace)

            J_plus = zero(type)
            J_minus = zero(type)

            for f in 1:length(Forces)
                # integration sign introduced here
                J_plus = JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                J_minus = -JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pzr[pz],pzr[pz+1],name_list[name])

                # scheme
                if scheme == "upwind"
                    if sign(J_plus) == 1
                        j_plus_right = 0
                        j_plus_left = 1
                    else
                        j_plus_right = 1
                        j_plus_left = 0
                    end
                    if sign(J_minus) == 1
                        j_minus_right = 1
                        j_minus_left = 0
                    else
                        j_minus_right = 0
                        j_minus_left = 1
                    end
                elseif scheme == "central"
                    j_plus_right = 0.5
                    j_plus_left = 0.5
                    j_minus_right = 0.5
                    j_minus_left = 0.5
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

                # normalised fluxes
                if b != bp
                    J_Flux[a,bp] += convert(T,(J_plus * j_plus_right) / ((pxr[px+1]-pxr[px])*(pyr[pyp+1]-pyr[pyp])*(pzr[pz+1]-pzr[pz])))
                    J_Flux[a,b] += convert(T,(J_plus * j_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                elseif typeof(BCp) == Open # b=bp
                    J_Flux[a,b] += convert(T,(J_plus * j_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                end
                if b != bm
                    J_Flux[a,b] += convert(T,(J_minus * j_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                    J_Flux[a,bm] += convert(T,(J_minus * j_minus_left) / ((pxr[px+1]-pxr[px])*(pyr[pym+1]-pyr[pym])*(pzr[pz+1]-pzr[pz])))
                elseif typeof(BCm) == Open # b=bm
                    J_Flux[a,b] += convert(T,(J_minus * j_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                end

            end # Force loop
        end
    end
end

function Fill_K_Flux!(K_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.zp_BC
    BCm = momentum_coords.zm_BC 
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    if typeof(BCp) != Periodic || typeof(BCm) != Periodic
        error("K flux boundaries incorrectly defined, i.e. not periodic")
    end

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = x_num*y_num*z_num

    tr = Grids.tr # time step assumed to be constant!
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr
    
    for name in 1:length(name_list)

        off_name = offset[name]
        pxr = Grids.pxr_list[name]
        pyr = Grids.pyr_list[name]
        pzr = Grids.pzr_list[name]

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,species_index,PhaseSpace)
            b = a 
            pzp = pz+1
            pzm = pz-1
            if pzp > pz_num # right boundary
                if typeof(BCp) == Closed || typeof(BCp) == Open
                    pzp = pz_num
                elseif typeof(BCp) == Periodic
                    pzp = 1
                end
            end
            if pzm < 1 # left boundary
                if typeof(BCm) == Closed || typeof(BCm) == Open
                    pzm = 1
                elseif typeof(BCm) == Periodic
                    pzm = pz_num
                end
            end

            bp = GlobalIndices_To_StateIndex(x,y,z,px,py,pzp,species_index,PhaseSpace)
            bm = GlobalIndices_To_StateIndex(x,y,z,px,py,pzm,species_index,PhaseSpace)

            K_plus = zero(type)
            K_minus = zero(type)

            for f in 1:length(Forces)
                # integration sign introduced here
                K_plus = KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz+1],name_list[name])
                K_minus = -KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],name_list[name])

                # scheme
                if scheme == "upwind"
                    if sign(K_plus) == 1
                        k_plus_right = 0
                        k_plus_left = 1
                    else
                        k_plus_right = 1
                        k_plus_left = 0
                    end
                    if sign(K_minus) == 1
                        k_minus_right = 1
                        k_minus_left = 0
                    else
                        k_minus_right = 0
                        k_minus_left = 1
                    end
                elseif scheme == "central"
                    k_plus_right = 1/2
                    k_plus_left = 1/2
                    k_minus_right = 1/2
                    k_minus_left = 1/2
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

                if b != bp
                    K_Flux[a,bp] += convert(T,(K_plus * k_plus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pzp+1]-pzr[pzp])))
                    K_Flux[a,b] += convert(T,(K_plus * k_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                elseif typeof(BCp) == Open # b=bp
                    K_Flux[a,b] += convert(T,(K_plus * k_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                end
                if b != bm
                    K_Flux[a,b] += convert(T,(K_minus * k_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                    K_Flux[a,bm] += convert(T,(K_minus * k_minus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pzm+1]-pzr[pzm])))
                elseif typeof(BCm) == Open # b=bm
                    K_Flux[a,b] += convert(T,(K_minus * k_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                end
                
            end # force loop
        end
    end
end
