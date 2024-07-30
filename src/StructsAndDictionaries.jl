CollisionMatricies = Dict({Vector{String},Tuple}) # creates an empty dictionary for storing collision matricies by inteaction name

mutable struct BoltzmannEquation <: Function

    t::Float32                  # the last timestep time to calculate Δt

    f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle

    ΔfS_list::Vector{Vector{Float32}}       # change in distribution function due to SMatrix
    ΔfT_list::Vector{Vector{Float32}}       # change in distribution function due to TMatrix

    function BoltzmannEquation(name_list)

        self = new()

        self.t = Float32(0)

        # initialize distribution function vectors for indvidual species
        num_species = length(name_list)
        self.f_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.f_list[i] = zeros(Float32,nump_list[i]*numt_list[i])
        end

        # initialize vectors for SMatrix and TMatrix changed so distribution functions for  indvidual species
        self.ΔfS_list = Vector{Vector{Float32}}(undef,num_species)
        self.ΔfT_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.ΔfS_list[i] = zeros(Float32,nump_list[i]*numt_list[i])
            self.ΔfT_list[i] = zeros(Float32,nump_list[i]*numt_list[i])
        end

        return self
    end

end