"""
    Euler(dg,g,t,dt)

Explicit and Implicit Euler time-stepping method for the Boltzmann equation. Explicit/Implicit is defined by the boolian `Implicit` flag when defining the EulerStruct.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````


# Implicit method:


"""
function (Euler::EulerStruct)(df::fType,f::fType,dt0,dt)

    # limit u to be positive, now done in solver
    #@. f = f*(f>=0f0)
    fill!(df,Float32(0))

    # reset arrays
    fill!(Euler.df,Float32(0))


    if Euler.Implicit
        fill!(Euler.Jac,Float32(0))
        fill!(Euler.temp,Float32(0))
        # TOFIX add fluxes inside mul!
        #mul!(g.df,g.M_Bin_Mul_Step,f)
        #@. g.Jac = 2*g.M_Bin_Mul_Step
        #g.df .= (I-dt*g.Jac)\g.df
        #@. df = dt*g.df
        
    else
        # this version is allocating, probably temp. Maybe put a temp array in Euler for this.
        fill!(Euler.temp,Float32(0))
        if isempty(Euler.PhaseSpace.Binary_list) == false
            update_Big_Bin!(Euler,f)
            @. Euler.temp += Euler.M_Bin_Mul_Step 
        end
        if isempty(Euler.PhaseSpace.Emi_list) == false
            update_Big_Emi!(Euler,f)
            @. Euler.temp += Euler.M_Emi_Step
        end
        # Signs here perhaps incorrect?!
        #  @. Euler.temp -= Euler.FluxM.B_Flux
        #  @. Euler.temp -= Euler.FluxM.C_Flux
        #  @. Euler.temp -= Euler.FluxM.D_Flux
        @. Euler.temp -= Euler.FluxM.K_Flux 
        @. Euler.temp -= Euler.FluxM.J_Flux 
        @. Euler.temp -= Euler.FluxM.I_Flux

        if isinf(sum(Euler.temp))
            error("overflow in arrays")
            #@. g.temp = g.temp*(g.temp!=Inf)
        end

        if Euler.PhaseSpace.Time.t_grid != "u" # non-uniform time stepping
            Euler.temp .*= dt / dt0
        end

        @. Euler.temp -= Euler.FluxM.Ap_Flux
        @. Euler.temp -= Euler.FluxM.Am_Flux

        #if Euler.Implicit

        #else
        #    Euler.LU = copy(Euler.FluxM.Ap_Flux)
        #    ldiv!(Euler.df,lu!(Euler.LU),Euler.temp)
        #end
        # improve allocations with LU decomp??
        Euler.df_temp .= Euler.FluxM.Ap_Flux * ones(Float32,size(Euler.temp,1)) # make vector as diagonal
        Euler.temp ./= Euler.df_temp
        #mul!(g.df,@view(g.FluxM.Ap_Flux[1,:,:])\g.temp,f)
        mul!(Euler.df,Euler.temp,f)
        @. df = Euler.df

        #println("$df")
        #error("")

    end

end

function update_Big_Bin!(method::SteppingMethodType,f)

    if size(method.BigM.M_Bin) != (length(f)^2,length(f))
        error("M_Bin is not the correct size")
    end

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    FluxM = method.FluxM 

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = FluxM.Vol

    # Thanks to Emma Godden for fixing a bug here
    temp = reshape(method.M_Bin_Mul_Step,length(f)*length(f))

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

        tempView = @view temp[n_momentum^2*off_space+1:n_momentum^2*(off_space+1)]
        fView = @view f[n_momentum*off_space+1:n_momentum*(off_space+1)]

        mul!(tempView,method.BigM.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        #println(maximum(method.BigM.M_Bin[:,168+144]))
        #println(method.M_Bin_Mul_Step[169:168+408,168+144])
        #println(method.M_Bin_Mul_Step[169:168+408,24])
        #println(Vol[off_space+1])
        #println(sum(fView))
        #println("")


        # multiply by volume element
        tempView .*= Vol[off_space+1]

        # assign jacobian elements
        if method.Implicit
            JacView = @view method.Jac[n_momentum*off_space+1:n_momentum*(off_space+1),n_momentum*off_space+1:n_momentum*(off_space+1)]
            @. JacView += 2*method.M_Bin_Mul_Step 
        end 

    end

    return nothing

end

function update_Big_Emi!(method::SteppingMethodType,f)

    if size(method.BigM.M_Emi) != (length(f),length(f))
        error("M_Bin is not the correct size")
    end

    PhaseSpace = method.PhaseSpace
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Momentum = PhaseSpace.Momentum
    FluxM = method.FluxM 

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_space = x_num+y_num+z_num
    n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

    Vol = FluxM.Vol

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

        tempView = @view method.M_Emi_Step[n_momentum*off_space+1:n_momentum*(off_space+1),n_momentum*off_space+1:n_momentum*(off_space+1)]

        tempView .= method.BigM.M_Emi

        # multiply by volume element
        tempView .*= Vol[off_space+1]

    end

    # assign jacobian elements
    if method.Implicit
        @. method.Jac += method.M_Emi_Step 
    end

    return nothing

end