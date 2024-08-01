CollisionMatricies = Dict{Vector{String},Tuple}() # creates an empty dictionary for storing collision matricies by inteaction name

mutable struct BoltzmannEquation <: Function

    t::Float32                  # the last timestep time to calculate Δt
    #diff_coeff::Float32         # the diffusion coefficient for the system

    f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle

    ΔfS_list::Vector{Vector{Float32}}       # change in distribution function due to SMatrix
    ΔfT_list::Vector{Vector{Float32}}       # change in distribution function due to TMatrix

    function BoltzmannEquation(u0)

        self = new()

        self.t = Float32(0)
        #self.diff_coeff = diff_coeff

        # initialize distribution function vectors for indvidual species
        num_species = length(name_list)
        self.f_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
        end

        # initialize vectors for SMatrix and TMatrix changed so distribution functions for  indvidual species
        self.ΔfS_list = Vector{Vector{Float32}}(undef,num_species)
        self.ΔfT_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.ΔfS_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
            self.ΔfT_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
        end

        return self
    end

end


#=
num_species = length(name_list)

mutable struct TempArrays

    f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle

    ΔfS_list::Vector{Vector{Float32}}       # change in distribution function due to SMatrix
    ΔfT_list::Vector{Vector{Float32}}       # change in distribution function due to TMatrix

end

function TempArraysInitialise(TempArrays)
    # initialize distribution function vectors for indvidual species
    num_species = length(name_list)
    TempArrays.f_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        TempArrays.f_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end

    # initialize vectors for SMatrix and TMatrix changed so distribution functions for  indvidual species
    TempArrays.ΔfS_list = Vector{Vector{Float32}}(undef,num_species)
    TempArrays.ΔfT_list = Vector{Vector{Float32}}(undef,num_species)
    for i in 1:num_species
        TempArrays.ΔfS_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
        TempArrays.ΔfT_list[i] = fill(Float32(0),nump_list[i]*numt_list[i])
    end
end




TempArraysInitialise(TempArrays)

=#