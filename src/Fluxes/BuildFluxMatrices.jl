"""
    BuildFluxMatrices(PhaseSpace)

Function that builds the flux matrices associated with coordinate forces and regular forces. First space is allocated for the arrays, then the fluxes are filled and finally the flux matrices are returned as an immutable `FluxMatricesStruct`.
"""
function BuildFluxMatrices(PhaseSpace::PhaseSpaceStruct;MatrixType::DataType=Matrix{Float32},VectorType::DataType=Vector{Float32})

    if eltype(MatrixType) != eltype(VectorType)
        error("elements of MatrixType and VectorType must be the same type")
    end

    (Ap_Flux,Am_Flux,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux,Vol) = Allocate_Flux(PhaseSpace,MatrixType,VectorType)

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

    return FluxMatricesStruct{MatrixType,VectorType}(Ap_Flux,Am_Flux,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux,Vol)

end

"""
    Allocate_Flux(PhaseSpace,MatrixType,VectorType)

Allocates arrays for fluxes and volume elements.
"""
function Allocate_Flux(PhaseSpace::PhaseSpaceStruct,MatrixType::DataType,VectorType::DataType)

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

    size = n*n*sizeof(eltype(MatrixType))

    println("Building flux matrices")
        if size > 1e9
            println("Flux is approx. $(size/1e9) GB in memory")
        elseif size > 1e6
            println("Flux is approx. $(size/1e6) MB in memory")
        elseif size > 1e3
            println("Flux is approx. $(size/1e3) KB in memory")
        else
            println("Flux is approx. $size bytes in memory")
        end

    # boundary terms included in arrays
    # fluxes allocated with all zeros
    # time fluxes
    Ap_Flux::MatrixType = zeros(eltype(MatrixType),n,n) 
    Am_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    # space fluxes
    B_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    C_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    D_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    # momentum fluxes
    I_Flux::MatrixType = zeros(eltype(MatrixType),n,n) 
    J_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    K_Flux::MatrixType = zeros(eltype(MatrixType),n,n)
    # volume element 
    Vol::VectorType = zeros(eltype(VectorType),n_space)

    return (Ap_Flux,Am_Flux,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux,Vol)

end

function Fill_A_Flux!(Ap_Flux::AbstractMatrix{<:AbstractFloat},Am_Flux::AbstractMatrix{<:AbstractFloat},PhaseSpace::PhaseSpaceStruct)

    type = eltype(Ap_Flux)

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

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

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

                # fill including normalisation
                Ap_Flux[a,b] += convert(type,A_plus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))
                Am_Flux[a,b] += convert(type,A_minus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz])))

            end
        end
    end
end

function Fill_Vol!(Vol::AbstractVector{<:AbstractFloat},PhaseSpace::PhaseSpaceStruct)

    type = eltype(Vol)

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

        Vol[space] = convert(type,VolFunction(space_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1]))

    end

end

function Fill_I_Flux!(I_Flux::AbstractMatrix{<:AbstractFloat},PhaseSpace::PhaseSpaceStruct)

    type = eltype(I_Flux)

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

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

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = (pz-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space

                a = (pz-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                b = a
                pxp = px+1
                pxm = px-1
                bp = (pz-1)*px_num*py_num+(py-1)*px_num+pxp+off_name+off_space
                bm = (pz-1)*px_num*py_num+(py-1)*px_num+pxm+off_name+off_space

                #= Boundary Conditions:
                    Flux on boundaries should be zero i.e. no particles leave/enter from the domain bounds
                =#
                if pxp > px_num
                    pxp = px_num
                    bp = (pz-1)*px_num*py_num+(py-1)*px_num+pxp+off_name+off_space
                    # b = bp therefore no I_plus flux only I_minus i.e. particles only leave/enter from left boundary (p<p_max)
                end
                if pxm < 1
                    pxm = 1
                    bm = (pz-1)*px_num*py_num+(py-1)*px_num+pxm+off_name+off_space
                    # b = bm therefore no I_minus flux only I_plus i.e. particles only leave/enter from right boundary (p_min<p)
                end

                I_plus = 0f0
                I_minus = 0f0

                for f in 1:length(Forces)
                    # integration sign introduced here
                    I_plus += IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                    I_minus -= IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
            
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

                    # normalised fluxes
                    if b != bp
                        I_Flux[a,bp] += convert(type,(I_plus * i_plus_right) / ((pxr[pxp+1]-pxr[pxp])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz])) 
                        I_Flux[a,b] += convert(type,(I_plus * i_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz])) 
                    end
                    if b != bm
                        I_Flux[a,b] += convert(type,(I_minus * i_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz])) 
                        I_Flux[a,bm] += convert(type,(I_minus * i_minus_left) / ((pxr[pxm+1]-pxr[pxm])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz])) 
                    end

                end # Forces loop
            end
        end
    end
end

function Fill_J_Flux!(J_Flux::AbstractMatrix{<:AbstractFloat},PhaseSpace::PhaseSpaceStruct)

    type = eltype(J_Flux)

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

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

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = (pz-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                b = a
                pyp = py+1
                pym = py-1
                bp = (pz-1)*px_num*py_num+(pyp-1)*px_num+px+off_name+off_space
                bm = (pz-1)*px_num*py_num+(pym-1)*px_num+px+off_name+off_space

                #= Boundary Conditions:
                    Flux on boundaries should be zero i.e. no particles leave/enter from the domain bounds
                =#
                if pyp > py_num
                    pyp = py_num
                    bp = (pz-1)*px_num*py_num+(pyp-1)*px_num+px+off_name+off_space
                    # b = bp therefore no J_plus flux only J_minus i.e. particles only leave/enter from left boundary (u<1)
                end
                if pym < 1
                    pym = 1
                    bm = (pz-1)*px_num*py_num+(pym-1)*px_num+px+off_name+off_space
                    # b = bm therefore no J_minus flux only J_plus i.e. particles only leave/enter from right boundary (u>-1)
                end

                J_plus = 0f0
                J_minus = 0f0

                for f in 1:length(Forces)
                    # integration sign introduced here
                    J_plus += JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                    J_minus -= JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pzr[pz],pzr[pz+1],name_list[name])

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

                    # normalised fluxes
                    if b != bp
                        J_Flux[a,bp] += convert(type,(J_plus * j_plus_right) / ((pxr[px+1]-pxr[px])*(pyr[pyp+1]-pyr[pyp]))*(pzr[pz+1]-pzr[pz]))
                        J_Flux[a,b] += convert(type,(J_plus * j_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz]))
                    end
                    if b != bm
                        J_Flux[a,b] += convert(type,(J_minus * j_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz]))
                        J_Flux[a,bm] += convert(type,(J_minus * j_minus_left) / ((pxr[px+1]-pxr[px])*(pyr[pym+1]-pyr[pym]))*(pzr[pz+1]-pzr[pz]))
                    end

                end # Force loop

            end
        end
    end
end

function Fill_K_Flux!(K_Flux::AbstractMatrix{<:AbstractFloat},PhaseSpace::PhaseSpaceStruct)

    type = eltype(K_Flux)

    Space = PhaseSpace.Space
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    Forces = PhaseSpace.Forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

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

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = (pz-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                b = a
                pzp = pz+1
                pzm = pz-1
                bp = (pzp-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                bm = (pzm-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space

                #= Boundary Conditions:
                    Flux on boundaries are periodic i.e. particles leave/enter from the opposite bound
                =#
                if pzp > pz_num
                    pzp = 1
                    bp = (pzp-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                end
                if pzm < 1
                    pzm = pz_num
                    bm = (pzm-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                end

                K_plus = 0f0
                K_minus = 0f0

                for f in 1:length(Forces)
                    # integration sign introduced here
                    K_plus += KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz+1],name_list[name])
                    K_minus -= KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],name_list[name])

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

                    if b != bp
                        K_Flux[a,bp] += convert(type,(K_plus * k_plus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pzp+1]-pzr[pzp]))
                        K_Flux[a,b] += convert(type,(K_plus * k_plus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz]))
                    end
                    if b != bm
                        K_Flux[a,b] += convert(type,(K_minus * k_minus_right) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pz+1]-pzr[pz]))
                        K_Flux[a,bm] += convert(type,(K_minus * k_minus_left) / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))*(pzr[pzm+1]-pzr[pzm]))
                    end
            end # force loop
            end
        end
    end
end
