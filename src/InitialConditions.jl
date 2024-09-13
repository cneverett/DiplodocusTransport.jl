function InitialConditions(Lists)

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    num_species = length(name_list)
    f0_list2D = Vector{Array{Float32,2}}(undef,num_species)
    for i in 1:num_species
        f0_list2D[i] = fill(Float32(0),nump_list[i],numt_list[i])
    end

    for i in eachindex(name_list)
        f0_list2D[i][40:41,:] .= 1f3 
    end

    f0_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        f0_list[i] = reshape(f0_list2D[i],nump_list[i]*numt_list[i])
    end

    u0 = zeros(Float32,sum(nump_list.*numt_list))

    f_list_to_u!(u0,f0_list)

    #test = reshape(u0,(nump_list[1],numt_list[1]))
    
    #test[40:41,:] .= 1f3

    #u0 .= reshape(test,prod(size(test)))

    return u0

end
