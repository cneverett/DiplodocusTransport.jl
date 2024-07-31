function InitialConditions()

    u0 = zeros(Float32,sum(nump_list.*numt_list))

    u0 .= 1f0

    return u0

end
