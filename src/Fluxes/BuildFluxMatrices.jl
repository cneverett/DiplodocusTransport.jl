function Allocate_Flux(Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with more coordinates

    Ap_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom) 
    Am_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom) 

    I_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom) # boundary terms included in arrays
    J_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom)
    K_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom)

    return (Ap_Flux,Am_Flux,I_Flux,J_Flux,K_Flux)

end

#= function Allocate_Flux_Force(Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with more coordinates

    I_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom) # boundary terms included in arrays
    J_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom)
    K_Flux::AbstractArray{Float32} = zeros(Float32,n_st,n_mom,n_mom)

    return (I_Flux,J_Flux,K_Flux)

end =#

function Build_Flux(FluxM::FluxMatrices,Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    # ensure empty to begin
    fill!(FluxM.Ap_Flux,0f0)
    fill!(FluxM.Am_Flux,0f0)
    fill!(FluxM.I_Flux,0f0)
    fill!(FluxM.J_Flux,0f0)
    fill!(FluxM.K_Flux,0f0)

    Fill_A_Flux!(FluxM.Ap_Flux,FluxM.Am_Flux,Lists,SpaceTime)

    Fill_I_Flux!(FluxM.I_Flux,Lists,SpaceTime) 
    Fill_J_Flux!(FluxM.J_Flux,Lists,SpaceTime)
    Fill_K_Flux!(FluxM.K_Flux,Lists,SpaceTime)  

end

function Fill_A_Flux!(Ap_Flux::Array{Float32},Am_Flux::Array{Float32},Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    spacetime_coords = SpaceTime.spacetime_coordinates
    momentum_coords = SpaceTime.momentum_coordinates

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list
    p_low_list = Lists.p_low_list
    p_up_list = Lists.p_up_list
    p_grid_list = Lists.p_grid_list
    u_grid_list = Lists.u_grid_list

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with coordinates

    # TOFIX: change coordinate grids later
    t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    x_r = [0,1f0] # cylindrical radial
    y_r = [0,Float32(2*pi)] # cylindrical theta
    z_r = [0,1f0] # cylindrical z

    phi0 = 0f0
    phi1 = Float32(2*pi)
    
    for i in 1:length(name_list)
        off = offset[i]
        name = name_list[i]
        p_r = BCI.bounds(p_low_list[i],p_up_list[i],p_num_list[i],p_grid_list[i])
        u_r = BCI.bounds(BCI.u_low,BCI.u_up,u_num_list[i],u_grid_list[i])

        for j in 1:n_st
            for k in 1:p_num_list[i], l in 1:u_num_list[i]

                a = (l-1)*p_num_list[i]+k+off
                b = a

                A_plus = AFluxFunction(spacetime_coords,momentum_coords,t_r[j+1],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l],u_r[l+1],phi0,phi1,name)
                A_minus = AFluxFunction(spacetime_coords,momentum_coords,t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l],u_r[l+1],phi0,phi1,name)
                Ap_Flux[j,a,b] += A_plus
                Am_Flux[j,a,b] -= A_minus

                # normalisation
                Ap_Flux[j,a,b] /= (p_r[k+1]-p_r[k])*(u_r[l+1]-u_r[l])#*(phi1-phi0)
                Am_Flux[j,a,b] /= (p_r[k+1]-p_r[k])*(u_r[l+1]-u_r[l])#*(phi1-phi0)
            end
        end
    end
end

function Fill_I_Flux!(I_Flux::Array{Float32},Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    spacetime_coords = SpaceTime.spacetime_coordinates
    momentum_coords = SpaceTime.momentum_coordinates

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list
    p_low_list = Lists.p_low_list
    p_up_list = Lists.p_up_list
    p_grid_list = Lists.p_grid_list
    u_grid_list = Lists.u_grid_list

    forces = SpaceTime.forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with coordinates

    # TOFIX: change coordinate grids later
    t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    x_r = [0,1f0] # cylindrical radial
    y_r = [0,Float32(2*pi)] # cylindrical theta
    z_r = [0,1f0] # cylindrical z

    phi0 = 0f0
    phi1 = Float32(2*pi)

    tmp = 0f0
    
    for i in 1:length(name_list)
        off = offset[i]
        name = name_list[i]
        p_r = BCI.bounds(p_low_list[i],p_up_list[i],p_num_list[i],p_grid_list[i])
        u_r = BCI.bounds(BCI.u_low,BCI.u_up,u_num_list[i],u_grid_list[i])

        for j in 1:n_st
            for k in 1:p_num_list[i], l in 1:u_num_list[i]

                a = (l-1)*p_num_list[i]+k+off
                b = a
                bp = b+1
                bm = b-1

                I_plus = 0f0
                I_minus = 0f0

                for f in 1:length(forces)
                    I_plus += IFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k+1],u_r[l],u_r[l+1],phi0,phi1,name)
                    I_minus += IFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],u_r[l],u_r[l+1],phi0,phi1,name)
                end

                # outflow momentum boundaries These may not be correct
                if bp > (l)*p_num_list[i]+off
                    bp = (l)*p_num_list[i]+off
                end
                if bm < (l-1)*p_num_list[i]+1+off
                    bm = (l-1)*p_num_list[i]+1+off
                end
                
                #=
                ________________________________
                a |  -I_m  | I_m-I_p | I_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#

                if b != bp
                    I_Flux[j,a,bp] += I_plus
                    I_Flux[j,a,b] -= I_plus
                    tmp += I_Flux[j,a,bp]
                end
                if b != bm
                    I_Flux[j,a,bm] -= I_minus
                    I_Flux[j,a,b] += I_minus
                    tmp += I_Flux[j,a,bm]
                end
                tmp += I_Flux[j,a,b]

                if tmp != 0f0
                    error("I_Flux not conservative")
                end
            end
        end
    end
end

function Fill_J_Flux!(J_Flux::Array{Float32},Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    spacetime_coords = SpaceTime.spacetime_coordinates
    momentum_coords = SpaceTime.momentum_coordinates

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list
    p_low_list = Lists.p_low_list
    p_up_list = Lists.p_up_list
    p_grid_list = Lists.p_grid_list
    u_grid_list = Lists.u_grid_list

    forces = SpaceTime.forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with coordinates

    # TOFIX: change coordinate grids later
    t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    x_r = [0,1f0] # cylindrical radial
    y_r = [0,Float32(2*pi)] # cylindrical theta
    z_r = [0,1f0] # cylindrical z

    phi0 = 0f0
    phi1 = Float32(2*pi)

    tmp = 0f0
    
    for i in 1:length(name_list)
        off = offset[i]
        name = name_list[i]
        p_r = BCI.bounds(p_low_list[i],p_up_list[i],p_num_list[i],p_grid_list[i])
        u_r = BCI.bounds(BCI.u_low,BCI.u_up,u_num_list[i],u_grid_list[i])

        for j in 1:n_st
            for k in 1:p_num_list[i], l in 1:u_num_list[i]

                a = (l-1)*p_num_list[i]+k+off
                b = a
                bp = b+1*p_num_list[i]
                bm = b-1*p_num_list[i]

                J_plus = 0f0
                J_minus = 0f0

                for f in 1:length(forces)
                    J_plus += JFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l+1],phi0,phi1,name)
                    J_minus += JFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l],phi0,phi1,name)
                end

                # Boundary conditions
                if bp > k+off+p_num_list[i]*(u_num_list[i]-1)
                    bp = k+off
                end
                if bm < k+off
                    bm = (u_num_list[i]-1)*p_num_list[i]+k+off
                end

                if b != bp
                    J_Flux[j,a,bp] += J_plus
                    J_Flux[j,a,b] -= J_plus
                    tmp += J_Flux[j,a,bp]
                end
                if b != bm
                    J_Flux[j,a,bm] -= J_minus
                    J_Flux[j,a,b] += J_minus
                    tmp += J_Flux[j,a,bm]
                end
                tmp += J_Flux[j,a,b]

                if tmp != 0f0
                    error("J_Flux not conservative")
                end
            end
        end
    end
end

function Fill_K_Flux!(K_Flux::Array{Float32},Lists::ListStruct,SpaceTime::SpaceTimeStruct)

    spacetime_coords = SpaceTime.spacetime_coordinates
    momentum_coords = SpaceTime.momentum_coordinates

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list
    p_low_list = Lists.p_low_list
    p_up_list = Lists.p_up_list
    p_grid_list = Lists.p_grid_list
    u_grid_list = Lists.u_grid_list

    forces = SpaceTime.forces

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    n_mom = sum(p_num_list.*u_num_list)

    n_st = 1 # will change with coordinates

    # TOFIX: change coordinate grids later
    t_r = SpaceTime.t_low:SpaceTime.dt:SpaceTime.t_up
    x_r = [0,1f0] # cylindrical radial
    y_r = [0,Float32(2*pi)] # cylindrical theta
    z_r = [0,1f0] # cylindrical z

    phi0 = 0f0
    phi1 = Float32(2*pi)
    
    for i in 1:length(name_list)
        off = offset[i]
        name = name_list[i]
        p_r = BCI.bounds(p_low_list[i],p_up_list[i],p_num_list[i],p_grid_list[i])
        u_r = BCI.bounds(BCI.u_low,BCI.u_up,u_num_list[i],u_grid_list[i])

        for j in 1:n_st
            for k in 1:p_num_list[i], l in 1:u_num_list[i]

                a = (l-1)*p_num_list[i]+k+off
                b = a

                K_plus = 0f0
                K_minus = 0f0

                for f in 1:length(forces)
                    K_plus += KFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l],u_r[l+1],phi1,name)
                    K_minus += KFluxFunction(forces[f],spacetime_coords,momentum_coords,t_r[j+1],t_r[j],x_r[1],x_r[2],y_r[1],y_r[2],z_r[1],z_r[2],p_r[k],p_r[k+1],u_r[l],u_r[l+1],phi0,name)
                end

                K_Flux[j,a,b] += K_plus
                K_Flux[j,a,b] -= K_minus
            end
        end
    end
end
