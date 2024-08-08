function InitialConditions(Lists)

    (name_list,nump_list,numt_list,pu_list,pl_list,interaction_list) = Lists

    u0 = zeros(Float32,sum(nump_list.*numt_list))

    test = reshape(u0,(nump_list[1],numt_list[1]))
    
    test[40,:] .= 1f0

    u0 .= reshape(test,prod(size(test)))

    return u0

end
