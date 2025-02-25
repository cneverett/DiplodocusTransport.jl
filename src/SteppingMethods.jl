function (g::Euler)(df::fType,f::fType,t,dt)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)
    fill!(df,Float32(0))

    # reset arrays
    fill!(g.df,Float32(0))
    
    # update changes due to Binary S and T interactions
    #if isempty(g.Lists.interaction_list_Binary) == false 
    #    update_Big_Binary!(g,f)
    #end

    if g.Implicit
        # TOFIX add fluxes inside mul!
        mul!(g.df,g.A_Binary_Mul,f)
        @. g.Jac = 2*g.A_Binary_Mul
        g.df .= (I-dt*g.Jac)\g.df
        @. df = dt*g.df
        
    else
        # this version is allocating, probably temp. Maybe put a temp array in Euler for this.
        fill!(g.temp,Float32(0))
        if isempty(g.Lists.interaction_list_Binary) == false
            update_Big_Binary!(g,f)
            @. g.temp += ((dt/(1f0))*2*pi/2) * g.A_Binary_Mul 
        end
        if isempty(g.Lists.interaction_list_Emi) == false
            @. g.temp += ((dt/(1f0))*2*pi/2) * g.BigM.A_Emi
        end
        @. g.temp -= @view(g.FluxM.Ap_Flux[1,:,:])
        @. g.temp -= @view(g.FluxM.Am_Flux[1,:,:])
        @. g.temp += @view(g.FluxM.K_Flux[1,:,:]) 
        @. g.temp += @view(g.FluxM.J_Flux[1,:,:]) 
        @. g.temp += @view(g.FluxM.I_Flux[1,:,:])

        if isinf(sum(g.temp))
            error("overflow in arrays")
            #@. g.temp = g.temp*(g.temp!=Inf)
        end

        # improve allocations with LU decomp??
        g.df_temp .= @view(g.FluxM.Ap_Flux[1,:,:]) * ones(Float32,size(g.temp,1)) # make vector as diagonal
        g.temp ./= g.df_temp
        #mul!(g.df,@view(g.FluxM.Ap_Flux[1,:,:])\g.temp,f)
        mul!(g.df,g.temp,f)
        @. df = g.df

    end

end

function update_Big_Binary!(g::SteppingMethod,f)

    # Thanks to Emma Godden for solving this bug
    temp = reshape(g.A_Binary_Mul,size(g.A_Binary_Mul,1)*size(g.A_Binary_Mul,2))
    mul!(temp,g.BigM.A_Binary,f) # temp is linked to A_Binary_Mul so it gets edited while maintaining is 2D shape
    #mul!(reshape(g.A_Binary_Mul,size(g.A_Binary_Mul,1)*size(g.A_Binary_Mul,2)),g.BigM.A_Binary,f)

    # below done in stepping method
    #if g.Implicit
    #   @. g.Jac += 2*g.A_Binary_Mul # assign jacobian elements
    #end
    #mul!(g.df_temp,temp,f)
    #@. g.df += g.df_temp 

    return nothing

end