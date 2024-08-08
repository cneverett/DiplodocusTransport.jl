# Empty dictionary for storing collision matricies by inteaction name
CollisionMatricies = Dict{Vector{String},Tuple}()

# Struct for storing the Boltzmann equation and its solution
mutable struct BoltzmannEquation <: Function

    t::Float32                  # the last timestep time to calculate Δt
    #diff_coeff::Float32         # the diffusion coefficient for the system

    f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle

    ΔfS_list::Vector{Vector{Float32}}       # change in distribution function due to SMatrix
    ΔfT_list::Vector{Vector{Float32}}       # change in distribution function due to TMatrix

    name_list::Vector{String}   # list of particle names
    nump_list::Vector{Int64}    # list of momentum bins for each particle
    numt_list::Vector{Int64}    # list of angular bins for each particle
    pu_list::Vector{Float32}    # list of upper momentum limits for each particle
    pl_list::Vector{Float32}    # list of lower momentum limits for each particle
    interaction_list::Vector{Vector{String}} # list of interactions

    function BoltzmannEquation(Lists)

        self = new()

        (self.name_list,self.nump_list,self.numt_list,self.pu_list,self.pl_list,self.interaction_list) = Lists

        self.t = Float32(0)
        #self.diff_coeff = diff_coeff

        # initialize distribution function vectors for indvidual species
        num_species = length(self.name_list)
        self.f_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.f_list[i] = fill(Float32(0),self.nump_list[i]*self.numt_list[i])
        end

        # initialize vectors for SMatrix and TMatrix changed so distribution functions for  indvidual species
        self.ΔfS_list = Vector{Vector{Float32}}(undef,num_species)
        self.ΔfT_list = Vector{Vector{Float32}}(undef,num_species)
        for i in 1:num_species
            self.ΔfS_list[i] = fill(Float32(0),self.nump_list[i]*self.numt_list[i])
            self.ΔfT_list[i] = fill(Float32(0),self.nump_list[i]*self.numt_list[i])
        end

        return self
    end

end