function (g::Euler)(df::fType,f::fType,t,dt)

    # limit u to be positive
    @. f = f*(f>=0f0)
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

        fill!(g.Jac,0f0)
        @. g.Jac = 2*g.A_Binary_Mul

        df .= (I-dt*g.J)\g.df
    else
        # TOFIX add fluxes inside mul!
        mul!(g.df_temp,g.A_Binary_Mul,f)
        @. df = dt*g.df_temp
    end

end

function update_Big_Binary!(g::SteppingMethod,f)

    # Thanks to Emma Godden for solving this bug
    #temp = reshape(g.A_Binary_Mul,size(g.A_Binary_Mul,1)*size(g.A_Binary_Mul,2))
    #mul!(temp,g.BigM.A_Binary,f) # temp is linked to A_Binary_Mul so it gets edited while maintaining is 2D shape
    mul!(reshape(g.A_Binary_Mul,size(g.A_Binary_Mul,1)*size(g.A_Binary_Mul,2)),g.BigM.A_Binary,f)

    # below done in stepping method
    #if g.Implicit
    #   @. g.Jac += 2*g.A_Binary_Mul # assign jacobian elements
    #end
    #mul!(g.df_temp,temp,f)
    #@. g.df += g.df_temp 

end