
figure_theme = Theme(backgroundcolor=(:black,0),
    Figure = (backgroundcolor=(:black,0),
        size=(600,600),
        fontsize = 20.0f0),
    Axis = (backgroundcolor=(:black,0),
        titlesize = 20.0f0,
        xlabelsize = 20.0f0,
        ylabelsize = 20.0f0,
        titlecolor= RGBAf(250, 250, 250, 1.0),
        ylabelcolor=RGBAf(250, 250, 250, 1.0),
        xlabelcolor=RGBAf(250, 250, 250, 1.0),
        bottomspinecolor=RGBAf(200, 200, 200, 1.0),
        leftspinecolor=RGBAf(200, 200, 200, 1.0),
        rightspinecolor=RGBAf(200, 200, 200, 1.0),
        topspinecolor=RGBAf(200, 200, 200, 1.0),
        spinewidth = 2.0,
        xgridcolor=RGBAf(100, 100, 100, 0.5),
        ygridcolor=RGBAf(100, 100, 100, 0.5),
        xgridwidth=1.0,
        ygridwidth=1.0,
        xminorgridcolor=RGBAf(100, 100, 100, 0.5),
        yminorgridcolor=RGBAf(100, 100, 100, 0.5),
        xminorgridwidth=1.0,
        yminorgridwidth=1.0,
        xtickcolor=RGBAf(200, 200, 200, 1.0),
        ytickcolor=RGBAf(200, 200, 200, 1.0),
        xticklabelcolor=RGBAf(200, 200, 200, 1.0),
        yticklabelcolor=RGBAf(200, 200, 200, 1.0),
        ),
    )
set_theme!(merge(theme_latexfonts(),figure_theme))
set_theme!()

"""
    FracNumPlot_AllSpecies(sol,num_species,dp_list,du_list,numInit_list,p_num_list,u_num_list;fig=nothing)

Returns a plot of the relative number density (compaired to initial values) of each species as a function of time.
"""
function FracNumPlot_AllSpecies(sol,num_species,pr_list,ur_list,numInit_list,p_num_list,u_num_list,mass_list;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Number Density",xlabel="Time",ylabel="Frac. Change")
    else
        ax = Axis(fig,title="Number Density",xlabel="Time",ylabel="Frac. Change")
    end

    for i in eachindex(sol.t)
        for j in 1:num_species

            Na = FourFlow(sol.u[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Ua = HydroFourVelocity(Na)
            num = ScalarNumberDensity(Na,Ua)

            scatter!(ax,sol.t[i],num/numInit_list[j]-1,marker = :circle,colormap = :viridis, colorrange = [1, num_species+1],color = j)
        end
    end

    return fig
end

"""
    NumPlot_AllSpecies(sol,num_species,dp_list,du_list,p_num_list,u_num_list;fig=nothing)

Returns a plot of the number density of each species as a function of time.
"""
function NumPlot_AllSpecies(sol,num_species,pr_list,ur_list,p_num_list,u_num_list,mass_list;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Number Density",xlabel="Time",ylabel=L"$\log_{10}$ Number Density $[\mathrm{m}^3]$")
    else
        ax = Axis(fig,title="Number Density",xlabel="Time",ylabel=L"$\log_{10}$ Number Density $[\mathrm{m}^3]$")
    end

    num = zeros(Float64,num_species)

    for i in eachindex(sol.t)
        for j in 1:num_species

            Na = FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Ua = HydroFourVelocity(Na)
            num[j] = ScalarNumberDensity(Na,Ua)

            if num[j] != 0
                scatter!(ax,sol.t[i],log10(num[j]),marker = :circle,colormap = :viridis, colorrange = [1, num_species+1],color = j)
            end 
    
        end

        num_tot = sum(num)
        if num_tot != 0
            scatter!(ax,sol.t[i],log10(num_tot),marker = :circle,color = :white)
        end

    end

    return fig
end

"""
    FracEnergyPlot_AllSpecies(sol,num_species,ΔE_list,dp_list,du_list,engInit_list,p_num_list,u_num_list;fig=nothing)

Returns a plot of the relative energy (compaired to initial value) of each species as a function of time.
"""
function FracEnergyPlot_AllSpecies(sol,num_species,pr_list,ur_list,engInit_list,p_num_list,u_num_list,mass_list;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Relative Energy",xlabel="Time",ylabel="Frac. Change")
    else
        ax = Axis(fig,title="Relatvie Energy",xlabel="Time",ylabel="Frac. Change")
    end

    for i in eachindex(sol.t)
        for j in 1:num_species

            Na = FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Ua = HydroFourVelocity(Na)
            num = ScalarNumberDensity(Na,Ua)

            Tab = StressEnergyTensor(sol.u[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])

            eng = ScalarEnergyDensity(Tab,Ua,num)

            scatter!(ax,sol.t[i],eng/engInit_list[j]-1,marker=:circle,colormap=:viridis,colorrange=[1,num_species+1],color=j)
        end
    end
    return fig
end

"""
    EnergyPlot_AllSpecies(sol,num_species,ΔE_list,dp_list,du_list,u_num_list;fig=nothing)

Returns a plot of the energy of each species as a function of time.
"""
function EnergyPlot_AllSpecies(sol,num_species,pr_list,ur_list,p_num_list,u_num_list,mass_list;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure(size=(600,300))
        ax = Axis(fig[1,1],title="Energy",xlabel="Time",ylabel=L"$\log_{10}$ Energy $[m_\text{Ele}c]$")
    else
        ax = Axis(fig,title="Relatvie Energy",xlabel="Time",ylabel=L"$\log_{10}$ Energy $[m_\text{Ele}c]$")
    end

    eng = zeros(Float64,num_species)

    for i in eachindex(sol.t)
        for j in 1:num_species

            Na = FourFlow(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
            Ua = HydroFourVelocity(Na)
            num = ScalarNumberDensity(Na,Ua)

            Tab = StressEnergyTensor(sol.f[i].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])

            eng[j] = ScalarEnergyDensity(Tab,Ua,num)

            if eng[j] != 0
                scatter!(ax,sol.t[i],log10(eng[j]),marker=:circle,colormap=:viridis,colorrange=[1,num_species+1],color=j)
            end
        end

        eng_tot = sum(eng)
        if eng_tot != 0
            scatter!(ax,sol.t[i],log10(eng_tot),marker = :circle,color = :white)
        end
    end
    return fig
end

"""
    PDistributionPlot_AllSpecies(sol,num_species,meanp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,outfreq;fig=nothing)

Returns a plot of the momentum distribution (distribution integrated over angles) of each species as a function of time.

# Optional Arguments
- `Flux`: Default false, if true and a photon population is present then that population is multiplied by the mean momentum, converting it from the momentum distribution to a flux i.e. ``F_\\nu``.
- `MaxwellJuttner`: Default true, if true the expected thermal distribtuion is plotted.
"""
function PDistributionPlot_AllSpecies(sol,num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt;mode="AXI",fig=nothing,animation=false,animation_end=0,time=0,MaxwellJuttner=true,Flux=false)

    if isnothing(fig)
        fig = Figure()
        time = sol.t[end]
        ax = Axis(fig[1,1],title="p-Distribution,   t=$time",xlabel=L"$\log_{10}$ Momentum $[m_\text{Ele}c]$",ylabel=L"$\log_{10}$ p-Distribution",limits=((nothing,nothing),(-20,5)))
    else
        if animation == true && length(sol.t) == 1
            ax = Axis(fig[1,1],title="p-Distribution,   t=$(round(time,sigdigits=5))",xlabel=L"$\log_{10}$ Momentum",ylabel=L"$\log_{10}$ p-Distribution",limits=((nothing,nothing),(-40,5)))
        elseif animation == true && length(sol.t) > 1
            ax = current_axis(fig)
            empty!(ax)
            ax.title = "P-Distribution,   t=$(round(time,sigdigits=5))"
        else
            ax = Axis(fig,title="p-Distribution",xlabel=L"$\log_{10}$ Momentum",ylabel=L"$\log_{10}$ p-Distribution",limits=((nothing,nothing),(-40,5)))
        end
    end

    max = -19
    min = -20
    # time dependent, not including initial state
    count = 1
    for i in 2:size(sol.t)[1]
        if isapprox(sol.t[i],count*out_dt) # first output in that range
            count += 1
            #println(sol.t[i])
            for j in 1:num_species
                my_colors = [cgrad(:grayC,rev=false)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
                # integrate distribution over μ
                if mode=="AXI"
                    dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
                    dist .= reshape(sol.f[i].x[j],(p_num_list[j],u_num_list[j]))
                    # unscale by dp*du 
                    for k in axes(dist,1), l in axes(dist,2)
                        dist[k,l] /= dp_list[j][k] * du_list[j][l]
                    end
                    if Flux == true && name_list[j] == "Pho"
                        for k in axes(dist,1)
                            dist[k,:] .*= meanp_list[j][k]
                        end
                    end
                    dist_p = dist*du_list[j]
                elseif mode=="ISO"
                    dist_p = sol.u[i].x[j] * 2
                end
                #stairs!(ax,log10.(meanp_list[j]),log10.(dist_p),color=my_colors[i],step=:center)
                stairs!(ax,log10.(meanp_list[j]),sign.(dist_p).*(log10.(abs.(dist_p)).+45),color=my_colors[i],step=:center)
                #println(dist_p)
                #stairs!(ax,log10.(meanp_list[j]),log10.(abs.(dist_p)),colormap = (:viridis,sol.t[i]/sol.t[end]),colorrange = [1,num_species+1], color = j,step=:center)
                
                maxlocal = maximum(replace!(log10.(abs.(dist_p)),NaN=>-19))
                if maxlocal == NaN
                    maxlocal = -19
                end
                max = maximum([max,maxlocal])
            end
        end
    end
    #println(max)

    # initial distribution
    for j in 1:num_species
        # integrate distribution over u
        if mode=="AXI"
            dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
                dist .= reshape(sol.f[1].x[j],(p_num_list[j],u_num_list[j]))
                # unscale by dp*du 
                for k in axes(dist,1), l in axes(dist,2)
                    dist[k,l] /= dp_list[j][k] * du_list[j][l]
                end
                dist_p = dist*du_list[j]
        elseif mode=="ISO"
            dist_p = sol.u[i].x[j] * 2
        end
        stairs!(ax,log10.(meanp_list[j]),log10.(dist_p),color=:blue,step=:center)
        maxlocal = maximum(log10.(dist_p))
        if maxlocal == NaN
            maxlocal = -19
        end
        max = maximum([max,maxlocal])
    end

    # expected distribution
    if animation == false || (animation == true && length(sol.t) == animation_end)
    if MaxwellJuttner == true
        for j in 1:num_species
            # MJ distribution for initial Temperature
            MJPoints = MaxwellJuttner_Distribution(meanp_list[j],tempInit_list[j],mass_list[j])
            #integrate distribution over μ
            if mode=="AXI"
                dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
                dist .= reshape(sol.f[end].x[j],(p_num_list[j],u_num_list[j]))
                # unscale by dp*du 
                for k in axes(dist,1), l in axes(dist,2)
                    dist[k,l] /= dp_list[j][k] * du_list[j][l]
                end
                dist_p = dist*du_list[j]
            elseif mode=="ISO"
                dist_p = sol.u[end].x[j] * 2
            end
            # plot expected distribution (normlised by last output)
            points = log10.((MJPoints/maximum(MJPoints)*maximum(dist_p)))
            maxlocal = maximum(points)
            if maxlocal == NaN || maxlocal == -Inf 
                maxlocal = -19
            end
            max = maximum([max,maxlocal])
            replace!(points,-Inf=>NaN) 
            scatterlines!(ax,log10.(meanp_list[j]),points,color=:red,markersize = 0.0)
        end
    end
    end

    # some straight lines
    #line_points = log10.(meanp_list[2].^(-2/3))
    #lines!(ax,log10.(meanp_list[2]),line_points,color=:red)
    #line_points = log10.(meanp_list[2].^(1/3))
    #lines!(ax,log10.(meanp_list[2]),line_points,color=:blue)
    #line_points = log10.(meanp_list[2].^(-2)).+7
    #lines!(ax,log10.(meanp_list[2]),line_points,color=:red)
    #line_points = log10.(meanp_list[2].^(-1/2)).-3
    #lines!(ax,log10.(meanp_list[2]),line_points,color=:green)
    #line_points = log10.(meanp_list[2].^(-3/2))
    #lines!(ax,log10.(meanp_list[2]),line_points,color=:green)

    ax.limits =((nothing,nothing),(-55,55)#= (min,max) =#)

    return fig
end

"""
    uDistributionPlot_AllSpecies(sol,num_species,meanμ_list,dp_list,p_num_list,u_num_list,out_dt;fig=nothing,animation=false,animation_end=0,time=0)

Returns a plot of the u distribution (distribution integrated over momenta) of each species as a function of time.
"""
function uDistributionPlot_AllSpecies(sol,num_species,meanu_list,dp_list,du_list,p_num_list,u_num_list,out_dt;fig=nothing,animation=false,animation_end=0,time=0)

    if isnothing(fig)
        fig = Figure()
        time = sol.t[end]
        ax = Axis(fig[1,1],title="u-Distribution,   t=$time",xlabel="u",ylabel="u-Distribution",limits=((-1,1),(nothing,nothing)))
    else
        if animation == true && length(sol.t) == 1
            ax = Axis(fig[1,1],title="u-Distribution,   t=$(round(time,sigdigits=5))",xlabel="u",ylabel="u-Distribution",limits=((-1,1),(-1e5,2.1e6)))
        elseif animation == true && length(sol.t) > 1
            ax = current_axis(fig)
            empty!(ax)
            ax.title = "u-Distribution,   t=$(round(time,sigdigits=5))"
            ax.limits = ((-1,1),(nothing,nothing))
        else
            ax = Axis(fig,title="u-Distribution",xlabel="u",ylabel="u-Distribution",limits=((-1,1),(nothing,nothing)))
        end
    end

    # time dependent
    count =1
    for i in 2:size(sol.t)[1]
        if isapprox(sol.t[i],count*out_dt)
            count +=1
            for j in 1:num_species
                my_colors = [cgrad(:grayC,rev=false)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
                # integrate distribution over p
                dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
                dist .= reshape(sol.f[i].x[j],(p_num_list[j],u_num_list[j]))
                # unscale by dp*du 
                for k in axes(dist,1), l in axes(dist,2)
                    dist[k,l] /= dp_list[j][k] * du_list[j][l]
                end
                dist_u = (dp_list[j]' * dist)'
                stairs!(ax,meanu_list[j],dist_u,color=my_colors[i],step=:center)
            end
        end
    end

    #initial distribution
    for j in 1:num_species
        dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
        dist .= reshape(sol.f[1].x[j],(p_num_list[j],u_num_list[j]))
        # unscale by dp*du 
        for k in axes(dist,1), l in axes(dist,2)
            dist[k,l] /= dp_list[j][k] * du_list[j][l]
        end
        dist_u = (dp_list[j]' * dist)'
        stairs!(ax,meanu_list[j],dist_u,color=:blue,step=:center)
    end
 
    return fig
end


"""
    puDistributionPlot_AllSpecies(sol,num_species,meanμ_list,meanp_list,p_num_list,u_num_list,outfreq;fig=nothing)

Returns a plot of the u distribution as a function of p of each species as a function of time.
"""
function puDistributionPlot_AllSpecies(sol,num_species,meanu_list,meanp_list,dp_list,du_list,p_num_list,u_num_list,out_dt;fig=nothing)

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="μ-Distribution as fun. of p",xlabel="p",ylabel="μ-Distribution",limits=((-14,4),(-1,1)))
    else
        ax = Axis(fig,title="u-Distribution as fun. of p",xlabel="p",ylabel="u-Distribution",limits=((-14,4),(-1,1)))
    end

    # time dependent
    count = 1
    for i in 2:size(sol.t)[1]
        if isapprox(sol.t[i],count*out_dt) # first output in that range
            count += 1
            my_colors = [cgrad(:grayC,rev=false)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
            for j in 1:num_species
                dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
                dist .= reshape(sol.f[i].x[j],(p_num_list[j],u_num_list[j]))
                # unscale by dp*du 
                for k in axes(dist,1), l in axes(dist,2)
                    dist[k,l] /= dp_list[j][k] * du_list[j][l]
                end
                dist_sec = zeros(Float32,u_num_list[j])
                for k in axes(dist,1)
                    dist_sec .= dist[k,:]                
                    minval = minimum(dist_sec)
                    maxval = maximum(dist_sec)
                    meanval = mean(dist_sec)
                    var = (maxval-minval)/meanval
                    scatterlines!(ax,dist_sec/meanval.-1 .+ log10(meanp_list[j][k]),meanu_list[j],color=my_colors[i],markersize=0.0)
                end
            end
        end
    end

    return fig
end


"""
    AllPlots(sol,num_species,dp_list,du_list,ΔE_list,meanp_list,meanμ_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt)

Returns a collation of figures for analysing the output of the simulation.
"""
function AllPlots(sol,num_species,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt;mode="AXI")

    if mode=="AXI"
        figure = Figure(size=(1000,1000))
    elseif mode=="ISO"
        figure = Figure(size=(1000,455))
    end

    NumPlot_AllSpecies(sol,num_species,pr_list,ur_list,numInit_list,p_num_list,u_num_list,mass_list,fig=figure[1,1],mode=mode)
    EnergyPlot_AllSpecies(sol,num_species,pr_list,ur_list,engInit_list,p_num_list,u_num_list,mass_list,fig=figure[1,2],mode=mode)
    PDistributionPlot_AllSpecies(sol,num_species,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt,fig=figure[2:3,1:2],mode=mode)
    if mode=="AXI"
        uDistributionPlot_AllSpecies(sol,num_species,meanu_list,dp_list,du_list,p_num_list,u_num_list,out_dt,fig=figure[4:5,1:2])
        puDistributionPlot_AllSpecies(sol,num_species,meanu_list,meanp_list,dp_list,du_list,p_num_list,u_num_list,out_dt,fig=figure[6:7,1:2])
    end

    return figure

end

"""
    PDistributionPlot_AllSpecies2(sol,num_species,meanp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,outfreq;fig=nothing)

Returns a plot of the momentum distribution (distribution integrated over angles) of each species as a function of time.
"""
function PDistributionPlot_AllSpecies2(sol,num_species,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    else
        ax = Axis(fig,title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    end

    max = 40
    min = -20
    # time dependent, not inlcuding initial state
    count = 1
    for i in 2:size(sol.t)[1]
        if isapprox(sol.t[i],count*out_dt) # first output in that range
            count += 1
            #println(sol.t[i])
            for j in 1:num_species
                my_colors = [cgrad(:grayC,rev=true)[z] for z ∈ range(0.0, 1.0, length = size(sol.t)[1])]
                # integrate distribution over μ
                if mode=="AXI"
                    distp = reshape(sol.u[i].x[j],(p_num_list[j],u_num_list[j]))
                    distp = dropdims(sum(distp,dims=2),dims=2)
                    distp ./= dp_list[j]
                elseif mode=="ISO"
                    distp = sol.u[i].x[j] * 2
                end
                stairs!(ax,log10.(meanp_list[j]),log10.(distp),color=my_colors[i],step=:center)
                maxlocal = maximum(log10.(distp))
                #max = maximum([max,maxlocal])
            end
        end
    end

    # initial distribution
    for j in 1:num_species
        # integrate distribution over μ
        if mode=="AXI"
            distp = reshape(sol.u[1].x[j],(p_num_list[j],u_num_list[j]))
            distp = dropdims(sum(distp,dims=2),dims=2)
            distp ./= dp_list[j]
        elseif mode=="ISO"
            distp = sol.u[1].x[j] /dp_list[i]
        end
        stairs!(ax,log10.(meanp_list[j]),log10.(distp),color=:blue,step=:center)
        maxlocal = maximum(log10.(distp))
        #max = maximum([max,maxlocal])
    end

    # expected distribution
    for j in 1:num_species
        # MJ distribution for initial Temperature
        MJPoints = MaxwellJuttner(meanp_list[j],tempInit_list[j],mass_list[j])
        #integrate distribution over μ
        if mode=="AXI"
            distp = reshape(sol.u[end].x[j],(p_num_list[j],u_num_list[j]))
            distp = dropdims(sum(distp,dims=2),dims=2)
            distp ./= dp_list[j]
        elseif mode=="ISO"
            distp = sol.u[end].x[j]/dp_list[i]
        end
        # plot expected distribution (normlised by last output)
        scatterlines!(ax,log10.(meanp_list[j]),log10.((MJPoints/maximum(MJPoints)*maximum(distp))),color=:red,markersize = 0.0)
        maxlocal = maximum(log10.((MJPoints/maximum(MJPoints)*maximum(distp))))
        #max = maximum([max,maxlocal])
    end
    ax.limits =((nothing,nothing),(min,max+1))

    return fig
end


function PDistributionPlot_AllSpecies_Ani(filename,sol,num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure()
        #ax = Axis(fig[1,1],title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    else
        ax = Axis(fig,title="p-Distribution",xlabel="Log10 Momentum",ylabel="Log10 p-Distribution",limits=((nothing,nothing),(-20,5)))
    end

    #time = Observable(0.0)

    # animation settings
    nframes = length(sol.t)
    framerate = 5
    itterator = range(1,nframes,step=1)
    record(fig,filename,itterator; framerate=framerate) do itter
        time = sol.t[itter]
        PDistributionPlot_AllSpecies(sol[1:itter],num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt;mode="AXI",fig=fig,animation=true,animation_end=nframes,time=time,MaxwellJuttner=false)
    end
    
end


function uDistributionPlot_AllSpecies_Ani(filename,sol,num_species,meanμ_list,dp_list,p_num_list,u_num_list,out_dt;mode="AXI",fig=nothing)

    if isnothing(fig)
        fig = Figure()
        #ax = Axis(fig[1,1],title="P-Distribution",xlabel="Log10 Momentum",ylabel="Log10 P-Distribution",limits=((nothing,nothing),(-20,5)))
    else
        ax = Axis(fig,title="u-Distribution",xlabel="u",ylabel="u-Distribution",limits=((-1,1),(nothing,nothing)))
    end

    #time = Observable(0.0)

    # animation settings
    nframes = length(sol.t)
    framerate = 5
    itterator = range(1,nframes,step=1)
    record(fig,filename,itterator; framerate=framerate) do itter
        time = sol.t[itter]
        uDistributionPlot_AllSpecies(sol[1:itter],num_species,meanμ_list,dp_list,du_list,p_num_list,u_num_list,out_dt;fig=fig,animation=true,animation_end=nframes,time=time)
    end
    
end


"""
    puDistributionPlot_AllSpecies_heatmap(sol,num_species,meanμ_list,meanp_list,p_num_list,u_num_list,outfreq;fig=nothing)

Returns a plot of the u distribution as a function of p of each species as a function of time.
"""
function puDistributionPlot_AllSpecies_heatmap(sol,t_int,num_species,dp_list,du_list,p_num_list,u_num_list,pr_list,ur_list;fig=nothing)

    if isnothing(fig)
        fig = Figure()
        #ax = Axis(fig[1,1],title="p-Distribution,   t=$time",xlabel=L"$\log_{10}$ Momentum $[m_\text{Ele}c]$",ylabel=L"Angle Cosine",limits=((-5,4),(-1,1)))
    else
        ax = Axis(fig,title="p-Distribution,   t=$time",xlabel=L"$\log_{10}$ Momentum $[m_\text{Ele}c]$",ylabel=L"$ $Angle Cosine $ $",limits=((-5,4),(-1,1)))
    end


    for j in 1:num_species

        ax  = Axis(fig[j,1])
        if j == 1
            ax.title = "p-Distribution,   t=$time"
        end
        dist = zeros(Float32,(p_num_list[j],u_num_list[j]))
        dist .= reshape(sol.f[t_int].x[j],(p_num_list[j],u_num_list[j]))
        # unscale by dp*du 
        for k in axes(dist,1), l in axes(dist,2)
            dist[k,l] /= dp_list[j][k] * du_list[j][l]
        end

        dist = log10.(dist)
        replace!(dist,-Inf32=>NaN)

    cmap = cgrad(:viridis,scale=:log10)
    hm = heatmap!(ax,log10.(pr_list[j]),ur_list[j],dist, colormap = cmap)

    Colorbar(fig[1:num_species,2],hm)

    end

    return fig
end


"""
    AllPlots(sol,num_species,dp_list,du_list,ΔE_list,meanp_list,meanμ_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,out_dt)

Returns a collation of figures for analysing the output of the simulation.
"""
function AllPlots_Ani(sol,num_species,name_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,filename,fps::Int;istart=nothing,istop=nothing,iframe=nothing,thermal=false,initial=false)

    fig = Figure(size=(1000,700,))
    #fig = Figure()

    # layout setup
    grid = fig[1,1] = GridLayout()
    gp = grid[1,2:3] = GridLayout()
    gu = grid[2,1] = GridLayout()
    gpu1 = grid[2,2:3] = GridLayout()
    #gpu = grid[3,2:3] = GridLayout()
    gc = grid[1,1] = GridLayout()
    gcn = grid[1,1][1,1:2] = GridLayout()
    gce = grid[1,1][2,1:2] = GridLayout()
    gl = grid[num_species+2,2:3] = GridLayout()

    # colors for particles
    my_colors = [cgrad(:roma)[z] for z ∈ range(0.0, 1.0, length = num_species+1)]
    # colors for heatmaps
    cmap = cgrad(:inferno,scale=:log10)

    # set up axes
        # u averaged distribution plot
        axp = Axis(gp[1,1],ylabel=L"$\log_{10}f(p)$",yaxisposition = :right,xticklabelsvisible = false,xticksvisible=false)
        # p averaged distribution plot
        axu = Axis(gu[1,1],xlabel=L"$\left(f(u)- \overline{f(u)}\right)/ \overline{f(u)}$",ylabel = L"$u$")
        # number density plot
        axn = Axis(gcn[1,1],ylabel=L"$\Delta n /n$",xticklabelsvisible = false,xticksvisible=false)
        # energy density plot
        axe = Axis(gce[1,1],ylabel=L"$\log_{10}(e)$",xticklabelsvisible = false,xticksvisible=false)

        ax1 = Axis(gpu1[1,1])
        ax1.ylabel= L"$u$"
        ax1.limits = ((-5,4),(-1,1))
        if num_species == 1
            ax1.xlabel=L"$\log_{10}(p)$  $[m_\text{Ele}c]$"
        else
            ax1.xticklabelsvisible = false
            ax1.xticksvisible=false
        end
        if 1 < num_species < 3
            ax2 = Axis(grid[3,2:3])
        end
        if 2 < num_species < 4
            ax3 = Axis(grid[3,2:3])
            ax2 = Axis(grid[4,2:3])
        end


        # nice labels on heatmaps
        for (i, label) in enumerate(name_list)
            if i == 1
            Box(grid[2,4], color = (my_colors[i],0.5),strokecolor = RGBAf(200, 200, 200, 1.0),strokewidth=2.0,cornerradius = (10,10,0,0),tellheight=false,width=56)
            Label(grid[2,4], label, rotation = pi/2, tellheight = false, color=RGBAf(200, 200, 200, 1.0),fontsize=24)
            else
            Box(grid[i+1,4], color = (my_colors[i],0.5),strokecolor = RGBAf(200, 200, 200, 1.0),strokewidth=2.0,cornerradius = (10,10,0,0),tellheight=false,width=56)
            Label(grid[i+1,4], label, rotation = pi/2, tellheight = false, color=RGBAf(200, 200, 200, 1.0),fontsize=20)
            end
        end

        # col and row gaps
        colgap!(grid,1,20)
        rowgap!(gc,1,10)
        if 1 < num_species
        rowgap!(grid,2,-30)
        end
        colgap!(grid,3,-53)

        #rowsize!(grid,1,Relative(1/(4)))
        #rowsize!(grid,2,Relative(3/5 / (num_species)))
        #rowsize!(grid,3,Relative(3/5 / (num_species)))
        #rowsize!(grid,4,Relative(1/5))
        
        # set up non mutating axes 
        for j in 1:num_species
            if j == 1
                d_null = zeros(Float32,(p_num_list[j]+1,u_num_list[j]))
                fill!(d_null,NaN)
                hm = heatmap!(ax1,[-14.0; log10.(pr_list[j])],ur_list[j],d_null, colormap = cmap) 
                hideydecorations!(ax1,grid=false)

            end
            if j == num_species && num_species != 1
                d_null = zeros(Float32,(p_num_list[j],u_num_list[j]))
                fill!(d_null,NaN)
                hm = heatmap!(ax2,log10.(pr_list[j]),ur_list[j],d_null, colormap = cmap)
                ax2.xlabel=L"$\log_{10}$ Momentum $[m_\text{Ele}c]$"
                ax2.yticklabelsvisible = false
                ax2.yticksvisible=false
                linkxaxes!(ax2,ax1) 
            end
            if j != 1 && j != num_species
                d_null = zeros(Float32,(p_num_list[j],u_num_list[j]))
                fill!(d_null,NaN)
                hm = heatmap!(ax3,log10.(pr_list[j]),ur_list[j],d_null, colormap = cmap)
                ax3.yticklabelsvisible = false
                ax3.yticksvisible=false
                ax3.xticklabelsvisible = false
                ax3.xticksvisible=false
                linkxaxes!(ax2,ax3) 
            end
        end

        linkxaxes!(ax1,axp) 
        linkyaxes!(ax1,axu)

        # temporary arrays
        num = zeros(Float32,num_species)
        num0 = zeros(Float32,num_species)
        eng = zeros(Float32,num_species)

        # time 
        time = Observable(0.0) 
        Label(grid[3,1],text=@lift("   t=$(round($time,sigdigits = 4))"),fontsize=24,tellheight=false,tellwidth=false,color=RGBAf(200,200,200,1.0),valign=:center,halign=:left)
        Label(grid[3,1],text=L"               $[1/n\sigma_{T} c]$",fontsize=24,tellheight=false,tellwidth=false,color=RGBAf(200,200,200,1.0),valign=:center,halign=:center)

        if isnothing(istart)
            istart = 1
        end
        if isnothing(istop)
            istop = length(sol.t)-1
        end

        d = Vector{Matrix{Float32}}(undef,num_species)
        dinitial = Vector{Matrix{Float32}}(undef,num_species)
        dthermal = Vector{Matrix{Float32}}(undef,num_species)
        for j in 1:num_species
            d[j] = zeros(Float32,(p_num_list[j],u_num_list[j]))
            dinitial[j] = zeros(Float32,(p_num_list[j],u_num_list[j]))
            dthermal[j] = zeros(Float32,(p_num_list[j],u_num_list[j]))
            dinitial[j] .= reshape(sol.f[1].x[j],(p_num_list[j],u_num_list[j]))
            println(istop)
            dthermal[j] .= reshape(sol.f[istop].x[j],(p_num_list[j],u_num_list[j]))
        end

        if thermal 
            tempFinal = zeros(Float32,num_species);
            for i in eachindex(name_list)
                Na = FourFlow(sol.f[istop].x[i],p_num_list[i],u_num_list[i],pr_list[i],ur_list[i],mass_list[i])
                #println(Na)
                Ua = HydroFourVelocity(Na)
                #println(Ua)
                Δab = ProjectionTensor(Ua)
                #println(Δab)
                Tab = StressEnergyTensor(sol.f[istop].x[i],p_num_list[i],u_num_list[i],pr_list[i],ur_list[i],mass_list[i])
                #println(Tab)
        
                n = ScalarNumberDensity(Na,Ua)
                #println(n)
                e = ScalarEnergyDensity(Tab,Ua,n)
                #println(e)
                p = ScalarPressure(Tab,Δab)
                #println(p)
                T = ScalarTemperature(p,n)
                #println(T)
        
                #numInit_list[i] = n
                #engInit_list[i] = e
                tempFinal[i] = T
            end

        end

    # animation settings
    itterator = range(istart,istop,step=1)
    framerate = fps
    record(fig,filename,itterator; framerate=framerate) do itter

        time[] = sol.t[itter]

        empty!(axp)
        empty!(axu)

        for j in 1:num_species
            dj = d[j]
            dj .= reshape(sol.f[itter].x[j],(p_num_list[j],u_num_list[j]))
            # unscale by dp*du 
            for k in axes(dj,1), l in axes(dj,2)
                #dj[k,l] /= (dp_list[j][k] * du_list[j][l])
            end
            @. dj = dj*(dj!=Inf)

            # u averaged distribution plot

                if initial
                    dp = log10.(dinitial[j]*du_list[j]) 
                    replace!(dp,-Inf32=>NaN)
                    replace!(dp,Inf32=>NaN)
                    stairs!(axp,log10.(meanp_list[j]),dp,color=:white,step=:center,linewidth=3.0)
                end

                if thermal
                    dpthermal = dthermal[j]*du_list[j] 
                    MJPoints = MaxwellJuttner_Distribution(meanp_list[j],tempFinal[j],mass_list[j])
                    # rescale as changed f definiton to number of particles
                    MJPoints .= (MJPoints .* dp_list[j]) # need to fix du_list later
                    # plot expected distribution (normlised by last output)
                    MJPoints .= log10.((MJPoints/maximum(MJPoints)*maximum(dpthermal)))
                    replace!(MJPoints,-Inf=>NaN) 
                    scatterlines!(axp,log10.(meanp_list[j]),MJPoints,color=:blue,markersize = 0.0)
                end

                dp = log10.(dj*du_list[j]) 
                replace!(dp,-Inf32=>NaN)
                replace!(dp,Inf32=>NaN)
                stairs!(axp,log10.(meanp_list[j]),dp,color=my_colors[j],step=:center,linewidth=3.0)

            # p averaged distribution plot

                if initial
                    du = (dp_list[j]' * dinitial[j])' 
                    stairs!(axu,([du[1]; du].-mean(du))./mean(du),ur_list[j][1:end],color=:white,step=:post,linewidth=3.0)
                end

                du = (dp_list[j]' * dj)'
                stairs!(axu,([du[1]; du].-mean(du))./mean(du),ur_list[j][1:end],color=my_colors[j],step=:post,linewidth=3.0)

            # distribution heat map
                dpu = log10.(dj)
                replace!(dpu,-Inf32=>NaN)

                if j == 1
                    empty!(ax1)
                    hm = heatmap!(gpu1[1,1],log10.(pr_list[j]),ur_list[j],dpu, colormap = cmap,colorrange=(-30,30),overdraw=true)
                    Colorbar(gl[1,1],hm,vertical=false,flipaxis = false,labelsize=20,label = L"$\log_{10}f(p,u)$",labelcolor=RGBAf(200,200,200,1.0),bottomspinecolor = RGBAf(200,200,200,1.0),topspinecolor = RGBAf(200,200,200,1.0),leftspinecolor = RGBAf(200,200,200,1.0),rightspinecolor = RGBAf(200,200,200,1.0),minortickcolor=RGBAf(100,100,100,1.0),tickcolor=RGBAf(200,200,200,1.0),ticklabelcolor=RGBAf(200,200,200,1.0),spinewidth=1.5)
                elseif j == num_species && num_species != 1
                    empty!(ax2)
                    hm = heatmap!(ax2,log10.(pr_list[j]),ur_list[j],dpu, colormap = cmap,colorrange=(-30,30),overdraw=true)
                else
                    hm = heatmap!(ax3,log10.(pr_list[j]),ur_list[j],dpu, colormap = cmap,colorrange=(-30,30),overdraw=true)
                end

            # frac number density 
                if itter == 1
                    num0[j] = numInit_list[j]
                else
                    #if isnothing(iframe)
                        num0[j] = num[j]
                    #else
                    #    num0[j] = numInit_list[j]+1f0 # to account for is 0
                    #end
                end
                Na = FourFlow(sol.f[itter].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])
                Ua = HydroFourVelocity(Na)
                num[j] = ScalarNumberDensity(Na,Ua)
#=                 if itter != 1
                    val = (num[j]-num0[j])/num[j]
                    #scatter!(axn,sol.t[itter],(num[j]-num0[j])/num0[j],marker = :circle,color = my_colors[j])
                    scatter!(axn,sol.t[itter],val,marker = :circle,color = my_colors[j])
                end =#

            # energy plot
                Tab = StressEnergyTensor(sol.f[itter].x[j],p_num_list[j],u_num_list[j],pr_list[j],ur_list[j],mass_list[j])

                eng[j] = ScalarEnergyDensity(Tab,Ua,num[j])

#=                 if eng[j] != 0
                    scatter!(axe,sol.t[itter],log10(eng[j]),marker=:circle,colormap=:viridis,color = my_colors[j])
                end =#  

        end # j loop 

        if itter != 1
            val = (sum(num)-sum(num0))/sum(num)
            #scatter!(axn,sol.t[itter],(num[j]-num0[j])/num0[j],marker = :circle,color = my_colors[j])
            scatter!(axn,sol.t[itter],val,marker = :circle,color = :white)

            if sum(eng) != 0
                scatter!(axe,sol.t[itter],log10(sum(eng)),marker=:circle,colormap=:viridis,color = :white)
            end
        end

        autolimits!(axp)
        autolimits!(axu)
        autolimits!(axe)
        autolimits!(axn)
        
    end
    if isnothing(iframe) == false
        return fig
    end

end
