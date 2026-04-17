mutable struct OutputStruct
    
    f::Vector{Vector{AbstractFloat}}
    t::Vector{AbstractFloat}

    function OutputStruct(f0::Vector{T},n_save::Int64) where T<:Union{Float32,Float64}

        self = new()
        self.f = [similar(f0) for _ in 1:n_save]
        self.t = Vector{T}(undef,n_save)

        return self
    
    end

end


function Solve(method::AbstractSteppingMethod,dt_initial::AbstractFloat,t_save::Vector{<:AbstractFloat};progress::Bool=false,fileName::String=nothing,fileLocation::String=pwd(),Verbose::Int64=1)

    if isdir(fileLocation) == false
        mkpath(fileLocation)
    end

    PhaseSpace = method.PhaseSpace

    f=zeros(eltype(method.f_init),size(method.f_init))
    copyto!(f,method.f_init)
    n_save = length(t_save)
    output = OutputStruct(f,n_save)

    # save initial state (step 1)
    @. output.f[1] = f
    output.t[1] = t_save[1]

    # progress bar
    if progress
        p = Progress(n_save)
    end

    @. method.f = method.f_init # reset f to initial condition at start of each solve (important for multiple solves in same session)
    #println(sum(method.f))
    
    dt = dt_initial

    for i in 2:n_save # start at 2 since initial state already saved

        t_start = t_save[i-1]
        t_stop = t_save[i]

        save = false
        # perform timestep
        while !save
            dt,save = method(t_start,t_stop,dt,Verbose)
            t_start += dt
        end

        # update t and dt and save state after timestep
        copyto!(f,method.f)
        #println("$(method.f)")
        @. output.f[i] = f
        output.t[i] = t_stop
        
        if progress
            next!(p)
        end
    end
    if progress
        finish!(p)
    end

    if !isnothing(fileName)
        # save output to file
        println("Saving Simulation output")

        SolutionFileSave(output,PhaseSpace,fileLocation,fileName);

    end


    return output

end

