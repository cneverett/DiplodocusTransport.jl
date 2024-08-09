function NumPlot_AllSpecies(sol,num_species,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig,title="Relative Number Density",xlabel="Time",ylabel="Frac. Change Num. Den.")

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

function TotMomentumPlot_AllSpecies(sol,num_species,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig,title="Total Momentum",xlabel="Time",ylabel="Log10 Momentum")

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

function TotEnergyPlot_AllSpecies(sol,num_species,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig,title="Relative Energy",xlabel="Time",ylabel="Frac. Change Energy")

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



function PDistributionPlot_AllSpecies(sol,num_species,meanp_list,Δμ_list,tempInit_list,mass_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig,title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))

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


function μDistributionPlot_AllSpecies(sol,num_species,meanμ_list,Δp_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in eachindex(f_list)
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig[1,1],title="μ-Distribution",xlabel="μ",ylabel="μ-Distribution",limits=((-1,1),(nothing,nothing)))

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



function pμDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,nump_list,numt_list,outfreq,fig=Figure())

    f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in eachindex(f_list)
        f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    ax = Axis(fig[1,1],title="μ-Distribution as fun. of p",xlabel="p",ylabel="μ-Distribution",limits=((-5,4),(-1,1)))

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



function AllPlots(sol,num_species,Δp_list,Δμ_list,ΔE_list,meanp_list,meanμ_list,numInit_list,engInit_list,tempInit_list,mass_list,nump_list,numt_list,outfreq)

    fig = Figure(size=(1000,1000))

    NumPlot_AllSpecies(sol,num_species,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq,fig[1,1])
    TotMomentumPlot_AllSpecies(sol,num_species,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq,fig[1,2])
    TotEnergyPlot_AllSpecies(sol,num_species,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq,fig[2,1])
    PDistributionPlot_AllSpecies(sol,num_species,meanp_list,Δμ_list,tempInit_list,mass_list,nump_list,numt_list,outfreq,fig[3:4,1:2])
    μDistributionPlot_AllSpecies(sol,num_species,meanμ_list,Δp_list,nump_list,numt_list,outfreq,fig[5:6,1:2])
    pμDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,nump_list,numt_list,outfreq,fig[7:8,1:2])

    #p3 = TempPlot(sol,dμ,dEkin,outfreq)
    #p4 = PDistributionPlot(sol,meanp,dμ,TInit,outfreq,f0)
    #p5 = TotEnergyPlot(sol,dμ,ΔE,outfreq,f0)
    #p6 = μDistributionPlot(sol,meanμ,dp,outfreq,f0)
    #p7 = pμDistributionPlot(sol,meanμ,meanp,outfreq)


    #l = @layout [grid(2,2); a{1.0w,0.2h}; a{1.0w,0.2h}; a{1.0w,0.2h}]

    #plot(p1, p2, p3, p5, p4,p6,p7, layout=l,size=(1000,1000), legend=false)

    return fig

end

