function (g::Euler)(df::fType,f::fType,t,dt)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)
    fill!(df,Float32(0))

    # reset arrays
    fill!(g.df,Float32(0))
    
    # update changes due to Binary S and T interactions
    if isempty(g.Lists.interaction_list_Binary) == false 
        update_Big_Binary!(g,f)
    end

    #if isempty(g.Lists.interaction_list_Emi) == false
    #    update_ΔS_Sync!(g,Matrices_Synchrotron,f)
    #    @. df += g.ΔfS_list
    #end

    if g.Implicit
        # TOFIX add fluxes inside mul!
        mul!(g.df,g.A_Binary_Mul,f)
        @. g.Jac = 2*g.A_Binary_Mul
        g.df .= (I-dt*g.Jac)\g.df
        @. df = dt*g.df
        
    else
        # TOFIX add fluxes inside mul!
        #mul!(g.df,g.A_Binary_Mul,f)
        #@. df = dt*g.df

        # this version is allocating, probably temp. Maybe put a temp array in Euler for this.
        temp = (dt*2*pi/2) .* g.A_Binary_Mul .- (@view(g.FluxM.Ap_Flux[1,:,:]) .+ @view(g.FluxM.Am_Flux[1,:,:])) .+ @view(g.FluxM.K_Flux[1,:,:]) .+ @view(g.FluxM.J_Flux[1,:,:]) .+  @view(g.FluxM.I_Flux[1,:,:])
        mul!(g.df,(@view(g.FluxM.Ap_Flux[1,:,:])\temp),f)
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

end