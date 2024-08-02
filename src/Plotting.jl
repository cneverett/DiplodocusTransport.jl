function NumPlot_AllSpecies(sol,f_list,Δp_list,Δμ_list,numInit_list,nump_list,numt_list,outfreq)

    fig = Figure()
    ax = Axis(fig[1,1],title="Relative Number Density",xlabel="Time",ylabel="Frac. Change Num. Den.")

    for i in eachindex(sol.t)
        if i%outfreq == 0
            Utof_list!(f_list,sol.u[i])
            for j in eachindex(f_list)
                num = NumberDensity(reshape(f_list[j],(nump_list[j],numt_list[j])),Δp_list[j],Δμ_list[j])
                scatter!(ax,sol.t[i],num/numInit_list[j]-1,marker = :circle,colormap = :viridis, colorrange = [1, length(f_list)+1],color = j)
                println(num)
            end
        end
    end

    return fig
end

function TotMomentumPlot_AllSpecies(sol,f_list,Δp_list,meanp_list,Δμ_list,nump_list,numt_list,outfreq)

    fig = Figure()
    ax = Axis(fig[1,1],title="Total Momentum",xlabel="Time",ylabel="Log10 Momentum")

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

function TotEnergyPlot_AllSpecies(sol,f_list,ΔE_list,Δp_list,Δμ_list,engInit_list,nump_list,numt_list,outfreq)

    fig = Figure()
    ax = Axis(fig[1,1],title="Relative Energy",xlabel="Time",ylabel="Frac. Change Energy")

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



function PDistributionPlot(sol,f_list,meanp_list,Δμ_list,TInit,outfreq,f0)
    p1=plot(legend=false,title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",ylims=(-30,5))
    # MJ distribution for initial Temperature
    MJPoints = MaxwellJuttnerPoints(meanp,TInit)
    distp = reshape(sol.u[end],(nump,numt))*dμ / sum(dμ)
    plot!(p1,log10.(meanp),log10.((MJPoints/maximum(MJPoints)*maximum(distp))),linewidth = 2, ls=:dot, lc=:red)
    # time dependent
    my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
    for i in axes(sol.t)[1]
        if i%outfreq == 0
        distp = reshape(sol.u[i],(nump,numt))*dμ
        plot!(p1,[log10.(meanp)],[log10.(distp)],linecolor=my_colors[i],seriestype=:stepmid)
        end
    end
    return p1
end

#=

function μDistributionPlot(sol,meanμ,dp,outfreq,f0)
    p1=plot(legend=false,title="μ-Distribution",xlabel="μ",ylabel="μ-Distribution",xlims=(-1,1))
    #initial distribution
    distμ = (dp'*reshape(f0,(nump,numt)))' / sum(dp)
    plot!(p1,[meanμ],[distμ],linewidth = 2, ls=:dot, lc=:red)
    # time dependent
    my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
    for i in axes(sol.t)[1]
        if i%outfreq == 0
        gamma = i/size(sol.t)[1]
        distμ = (dp'*reshape(sol.u[i],(nump,numt)))'
        plot!(p1,[meanμ],[distμ],linecolor=my_colors[i],seriestype=:stepmid)
        end
    end  
    return p1
end


function pμDistributionPlot(sol,meanμ,meanp,outfreq)
    p1=plot(legend=false,title="μ-Distribution as fun. of p",xlabel="p",ylabel="μ-Distribution",xlims=(-5,4),color=:gray,ylims=(-1,1))
    for j in axes(sol.t)[1]
        if j%outfreq == 0
            my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
            for i in 1:size(meanp)[1]
                distμ = reshape(sol.u[j],(nump,numt))[i,:]
                minval = minimum(distμ)
                maxval = maximum(distμ)
                meanval = mean(distμ)
                var = (maxval-minval)/meanval
                plot!(p1,[distμ./mean(distμ).-1 .+ log10(meanp[i])],[meanμ],linecolor=my_colors[j])
            end
        end
    end
    return p1
end


function AllPlots(sol,f0,dpdμ,meanp,dμ,dEkin,ΔE,outfreq)
    numInit = numberDensity(reshape(f0,(nump,numt)),dpdμ)
    momInit = momentum(reshape(f0,(nump,numt)),dpdμ,meanp,numInit)
    TInit = temperature(reshape(f0,(nump,numt)),dEkin,dμ,numInit)

    p1 = NumPlot(sol,dpdμ,outfreq,f0)
    p2 = TotMomentumPlot(sol,dpdμ,meanp,outfreq)
    p3 = TempPlot(sol,dμ,dEkin,outfreq)
    p4 = PDistributionPlot(sol,meanp,dμ,TInit,outfreq,f0)
    p5 = TotEnergyPlot(sol,dμ,ΔE,outfreq,f0)
    p6 = μDistributionPlot(sol,meanμ,dp,outfreq,f0)
    p7 = pμDistributionPlot(sol,meanμ,meanp,outfreq)


    l = @layout [grid(2,2); a{1.0w,0.2h}; a{1.0w,0.2h}; a{1.0w,0.2h}]

    plot(p1, p2, p3, p5, p4,p6,p7, layout=l,size=(1000,1000), legend=false)

end

=#