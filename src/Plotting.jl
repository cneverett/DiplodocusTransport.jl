"""
    NumPlot_AllSpecies(sol,num_species,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the relative number density (compaired to initial value) of each species as a function of time.
"""
function NumPlot_AllSpecies(sol,num_species,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="Relative Number Density",xlabel="Time",ylabel="Frac. Change")
    else
        ax = Axis(fig,title="Relative Number Density",xlabel="Time",ylabel="Frac. Change")
    end

    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in eachindex(f_list)
                num = NumberDensity(reshape(f_list[j],(nump_list[j],numt_list[j])),Δp_list[j],Δμ_list[j])
                scatter!(ax,sol.t[i],num/numInit_list[j]-1,marker = :circle,colormap = :viridis, colorrange = [1, length(f_list)+1],color = j)
            end
        end
    end

    return fig
end

"""
    TotMomentumPlot_AllSpecies(sol,num_species,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the total momentum of each species as a function of time.
"""
function TotMomentumPlot_AllSpecies(sol,num_species,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="Total Momentum",xlabel="Time",ylabel="Log10 Momentum")
    else
        ax = Axis(fig,title="Total Momentum",xlabel="Time",ylabel="Log10 Momentum")
    end
    

    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in eachindex(f_list)
                f_reshape = reshape(f_list[j],(nump_list[j],numt_list[j]))
                num = NumberDensity(f_reshape,Δp_list[j],Δμ_list[j])
                mom = Momentum(f_reshape,Δp_list[j],meanp_list[j],Δμ_list[j],num)
            scatter!(ax,sol.t[i],log10(mom),marker=:circle,colormap=:viridis,colorrange=[1,length(f_list)+1],color=j)
            end
        end
    end
    return fig
end

"""
    EnergyPlot_AllSpecies(sol,num_species,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the relative energy (compaired to initial value) of each species as a function of time.
"""
function EnergyPlot_AllSpecies(sol,num_species,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="Relative Energy",xlabel="Time",ylabel="Frac. Change")
    else
        ax = Axis(fig,title="Relatvie Energy",xlabel="Time",ylabel="Frac. Change")
    end

    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in eachindex(f_list)
                f_reshape = reshape(f_list[j],(nump_list[j],numt_list[j]))
                num = NumberDensity(f_reshape,Δp_list[j],Δμ_list[j])
                eng = Energy(f_reshape,ΔE_list[j],Δμ_list[j],num)
                scatter!(ax,sol.t[i],eng/engInit_list[j]-1,marker=:circle,colormap=:viridis,colorrange=[1,length(f_list)+1],color=j)
            end
        end
    end
    return fig
end

#=


function TempPlot(sol,dμ,dEkin,outfreq)
    T=0f0
    num = 0f0
    p1=scatter(legend=false,title="Temperature",xlabel="Time",ylabel="Log10 Temperature")
    for i in axes(sol.t)[1]
        if i%outfreq == 0
        num = numberDensity(reshape(sol.u[i],(nump,numt)),dpdμ)
        T = temperature(reshape(sol.u[i],(nump,numt)),dEkin,dμ,num)
        scatter!(p1,[sol.t[i]],[log10(T)],marker=(:circle,2,:black))
        end
    end
    return p1
end

=#


"""
    PDistributionPlot_AllSpecies(sol,num_species,meanp_list,Δμ_list,tempInit_list,mass_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the momentum distribution (distribution integrated over angles) of each species as a function of time.
"""
function PDistributionPlot_AllSpecies(sol,num_species,meanp_list,Δμ_list,tempInit_list,mass_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    else
        ax = Axis(fig,title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    end

    # time dependent
    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in 1:num_species
                my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
            # integrate distribution over μ
            distp = reshape(f_list[j],(nump_list[j],numt_list[j]))*Δμ_list[j]
            stairs!(ax,log10.(meanp_list[j]),log10.(distp),color=my_colors[i],step=:center)
            end
        end
    end

    # expected distribution
    Utof_list!(f_list,sol.u[end]) # final distribution for normalisation
    for j in 1:num_species
        # MJ distribution for initial Temperature
        MJPoints = MaxwellJuttner(meanp_list[j],tempInit_list[j],mass_list[j])
        #integrate distribution over μ
        distp = reshape(f_list[j],(nump_list[j],numt_list[j]))*Δμ_list[j]
        # plot expected distribution (normlised by last output)
        scatterlines!(ax,log10.(meanp_list[j]),log10.((MJPoints/maximum(MJPoints)*maximum(distp))),color=:red,markersize = 0.0)
    end

    return fig
end

"""
    μDistributionPlot_AllSpecies(sol,num_species,meanμ_list,Δp_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the μ distribution (distribution integrated over momenta) of each species as a function of time.
"""
function μDistributionPlot_AllSpecies(sol,num_species,meanμ_list,Δp_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in eachindex(f_list)
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="μ-Distribution",xlabel="μ",ylabel="μ-Distribution",limits=((-1,1),(nothing,nothing)))
    else
        ax = Axis(fig[1,1],title="μ-Distribution",xlabel="μ",ylabel="μ-Distribution",limits=((-1,1),(nothing,nothing)))
    end

    # time dependent
    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in 1:num_species
                my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
                # integrate distribution over p
                distμ = (Δp_list[j]' * reshape(f_list[j],(nump_list[j],numt_list[j])))'
                stairs!(ax,meanμ_list[j],distμ,color=my_colors[i],step=:center)
            end
        end
    end

    #initial distribution
    Utof_list!(f_list,sol.u[1])
    for j in 1:num_species
        distμ = (Δp_list[j]' * reshape(f_list[j],(nump_list[j],numt_list[j])))'
        stairs!(ax,meanμ_list[j],distμ,color=:red,step=:center,linestyle=:dash)
    end
 
    return fig
end


"""
    pμDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,nump_list,numt_list,outfreq;fig=nothing)

Returns a plot of the μ distribution as a function of p of each species as a function of time.
"""
function pμDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,nump_list,numt_list,outfreq;fig=nothing)

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in eachindex(f_list)
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="μ-Distribution as fun. of p",xlabel="p",ylabel="μ-Distribution",limits=((-5,4),(-1,1)))
    else
        ax = Axis(fig,title="μ-Distribution as fun. of p",xlabel="p",ylabel="μ-Distribution",limits=((-5,4),(-1,1)))
    end

    # time dependent
    for i in eachindex(sol.t)
        if i%outfreq == 0
            my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
            Utof_list!(f_list,sol.u[i])
            for j in 1:num_species
                for k in 1:size(meanp_list[j])[1]
                    distμ = reshape(f_list[j],(nump_list[j],numt_list[j]))[k,:]
                    minval = minimum(distμ)
                    maxval = maximum(distμ)
                    meanval = mean(distμ)
                    var = (maxval-minval)/meanval
                    scatterlines!(ax,distμ./mean(distμ).-1 .+ log10(meanp_list[j][k]),meanμ_list[j],color=my_colors[i],markersize=0.0)
                end
            end
        end
    end

    return fig
end


"""
    AllPlots(sol,num_species,Δp_list,Δμ_list,ΔE_list,meanp_list,meanμ_list,numInit_list,engInit_list,tempInit_list,mass_list,nump_list,numt_list,outfreq)

Returns a collation of figures for analysing the output of the simulation.
"""
function AllPlots(sol,num_species,Δp_list,Δμ_list,ΔE_list,meanp_list,meanμ_list,numInit_list,engInit_list,tempInit_list,mass_list,nump_list,numt_list,outfreq)

    figure = Figure(size=(1000,1000))

    NumPlot_AllSpecies(sol,num_species,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq,fig=figure[1,1])
    #TotMomentumPlot_AllSpecies(sol,num_species,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq,fig[1,2])
    EnergyPlot_AllSpecies(sol,num_species,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq,fig=figure[1,2])
    PDistributionPlot_AllSpecies(sol,num_species,meanp_list,Δμ_list,tempInit_list,mass_list,nump_list,numt_list,outfreq,fig=figure[2:3,1:2])
    μDistributionPlot_AllSpecies(sol,num_species,meanμ_list,Δp_list,nump_list,numt_list,outfreq,fig=figure[4:5,1:2])
    pμDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,nump_list,numt_list,outfreq,fig=figure[6:7,1:2])

    return figure

end

