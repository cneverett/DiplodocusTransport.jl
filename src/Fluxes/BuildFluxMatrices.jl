function Allocate_Flux(PhaseSpace::PhaseSpaceStruct)

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
    # time fluxes
    Ap_Flux::Array{Float32,2} = zeros(Float32,n,n) 
    Am_Flux::Array{Float32,2} = zeros(Float32,n,n)
    # space fluxes
    B_Flux::Array{Float32,2} = zeros(Float32,n,n)
    C_Flux::Array{Float32,2} = zeros(Float32,n,n)
    D_Flux::Array{Float32,2} = zeros(Float32,n,n)
    # momentum fluxes
    I_Flux::Array{Float32,2} = zeros(Float32,n,n) 
    J_Flux::Array{Float32,2} = zeros(Float32,n,n)
    K_Flux::Array{Float32,2} = zeros(Float32,n,n)
    # volume element 
    Vol::Vector{Float32} = zeros(Float32,n_space)

    return (Ap_Flux,Am_Flux,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux,Vol)

end

function Build_Flux(FluxM::FluxMatricesStruct,PhaseSpace::PhaseSpaceStruct)

    # ensure empty to begin
    fill!(FluxM.Ap_Flux,0f0)
    fill!(FluxM.Am_Flux,0f0)
    #fill!(FluxM.B_Flux,0f0)
    #fill!(FluxM.C_Flux,0f0)
    #fill!(FluxM.D_Flux,0f0)
    fill!(FluxM.I_Flux,0f0)
    fill!(FluxM.J_Flux,0f0)
    fill!(FluxM.K_Flux,0f0)
    fill!(FluxM.Vol,0f0)

    Fill_A_Flux!(FluxM.Ap_Flux,FluxM.Am_Flux,PhaseSpace)

    # to be implemented
    #Fill_B_Flux!(FluxM.B_Flux,Lists,PhaseSpace)
    #Fill_C_Flux!(FluxM.C_Flux,Lists,PhaseSpace)
    #Fill_D_Flux!(FluxM.D_Flux,Lists,PhaseSpace)

    Fill_I_Flux!(FluxM.I_Flux,PhaseSpace) 
    Fill_J_Flux!(FluxM.J_Flux,PhaseSpace)
    Fill_K_Flux!(FluxM.K_Flux,PhaseSpace) 
    
    Fill_Vol!(FluxM.Vol,PhaseSpace)

end

function Fill_A_Flux!(Ap_Flux::Array{Float32},Am_Flux::Array{Float32},PhaseSpace::PhaseSpaceStruct)

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

                A_plus = AFluxFunction(space_coords,momentum_coords,tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                A_minus = AFluxFunction(space_coords,momentum_coords,tr[1],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])

                Ap_Flux[a,b] += A_plus
                Am_Flux[a,b] -= A_minus

                # normalisation
                Ap_Flux[a,b] /= (pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])#*(phi1-phi0)
                Am_Flux[a,b] /= (pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py])#*(phi1-phi0)

            end
        end
    end
end

function Fill_Vol!(Vol::Vector{Float32},PhaseSpace::PhaseSpaceStruct)

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

        Vol[space] = VolFunction(space_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1])

    end

end

function Fill_I_Flux!(I_Flux::Array{Float32},PhaseSpace::PhaseSpaceStruct)

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

    # TOFIX: change coordinate grids later
    #t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    #x_r = [0,1f0] # cylindrical radial
    #y_r = [0,Float32(2*pi)] # cylindrical theta
    #z_r = [0,1f0] # cylindrical z

    #phi0 = 0f0
    #phi1 = Float32(2*pi)
    
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

                I_plus = 0f0
                I_minus = 0f0

                for f in 1:length(Forces)
                    I_plus += IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                    I_minus += IFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pyr[py],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                end

                # outflow momentum boundaries These may not be correct
                if pxp > px_num
                    pxp = px_num
                    bp = (pz-1)*px_num*py_num+(py-1)*px_num+pxp+off_name+off_space
                end
                if pxm < 1
                    pxm = 1
                    bm = (pz-1)*px_num*py_num+(py-1)*px_num+pxm+off_name+off_space
                end
                
                #=
                ________________________________
                a |  -I_m  | I_m-I_p | I_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#
                if scheme == "central"
                    if b != bp
                        I_Flux[a,bp] += I_plus / ((pxr[px+2]-pxr[px+1])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                        I_Flux[a,b] -= I_plus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                    end
                    if b != bm
                        I_Flux[a,bm] -= I_minus / ((pxr[px]-pxr[px-1])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                        I_Flux[a,b] += I_minus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                    end
                else
                    error("Scheme note recognised")
                end
            end
        end
    end
end

function Fill_J_Flux!(J_Flux::Array{Float32},PhaseSpace::PhaseSpaceStruct)

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

    # TOFIX: change coordinate grids later
    #t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    #x_r = [0,1f0] # cylindrical radial
    #y_r = [0,Float32(2*pi)] # cylindrical theta
    #z_r = [0,1f0] # cylindrical z

    #phi0 = 0f0
    #phi1 = Float32(2*pi)
    
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

                J_plus = 0f0
                J_minus = 0f0

                for f in 1:length(Forces)
                    J_plus += JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py+1],pzr[pz],pzr[pz+1],name_list[name])
                    J_minus += JFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pzr[pz],pzr[pz+1],name_list[name])
                end

                # Flux on boundaries should be zero
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

                #=
                ________________________________
                a |  -J_m  | J_m-J_p | J_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#

                if b != bp
                    J_Flux[a,bp] += J_plus / ((pxr[px+1]-pxr[px])*(pyr[pyp+1]-pyr[pyp]))#*(phi1-phi0)
                    J_Flux[a,b] -= J_plus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                end
                if b != bm
                    J_Flux[a,bm] -= J_minus / ((pxr[px+1]-pxr[px])*(pyr[pym+1]-pyr[pym]))#*(phi1-phi0)
                    J_Flux[a,b] += J_minus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                end

            end
        end
    end
end

function Fill_K_Flux!(K_Flux::Array{Float32},PhaseSpace::PhaseSpaceStruct)

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

    # TOFIX: change coordinate grids later
    #t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    #x_r = [0,1f0] # cylindrical radial
    #y_r = [0,Float32(2*pi)] # cylindrical theta
    #z_r = [0,1f0] # cylindrical z

    #phi0 = 0f0
    #phi1 = Float32(2*pi)
    
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

                K_plus = 0f0
                K_minus = 0f0

                for f in 1:length(Forces)
                    K_plus += KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz+1],name_list[name])
                    K_minus += KFluxFunction(Forces[f],space_coords,momentum_coords,tr[1],tr[2],xr[x],xr[x+1],yr[y],yr[y+1],zr[z],zr[z+1],pxr[px],pxr[px+1],pyr[py],pyr[py+1],pzr[pz],name_list[name])
                end

                # Boundary conditions periodic these may not be correct!
                if pzp > pz_num
                    pzp = 1
                    bp = (pzp-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                end
                if pzm < 1
                    pzm = pz_num
                    bm = (pzm-1)*px_num*py_num+(py-1)*px_num+px+off_name+off_space
                end

                #=
                ________________________________
                a |  -K_m  | K_m-K_p | K_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#

                if b != bp
                    K_Flux[a,bp] += K_plus / ((pxr[px+2]-pxr[px+1])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                    K_Flux[a,b] -= K_plus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                end
                if b != bm
                    K_Flux[a,bm] -= K_minus / ((pxr[px]-pxr[px-1])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                    K_Flux[a,b] += K_minus / ((pxr[px+1]-pxr[px])*(pyr[py+1]-pyr[py]))#*(phi1-phi0)
                end
            end
        end
    end
end
