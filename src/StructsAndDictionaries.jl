# Empty dictionary for storing binary collision matricies by inteaction name
CollisionMatriciesBinary = Dict{Vector{String},Tuple}()
CollisionMatriciesSync = Dict{Vector{String},Array{Float32,2}}()

# Struct for storing the Boltzmann equation and its solution
mutable struct BoltzmannEquation <: Function

    t::Float32                  # the last timestep time to calculate Δt

    #f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle
    f1DA::ArrayPartition  # advanced distribution function 
    f1DR::ArrayPartition  # retarded distributin function
    state::Bool

    ΔfS_list::ArrayPartition       # change in distribution function due to SMatrix
    ΔfT_list::ArrayPartition       # change in distribution function due to TMatrix
    ΔfS_list_temp::ArrayPartition       # temporary array for change in distribution function due to SMatrix
    ΔfS_mul_step::ArrayPartition       # temporary array the matrix multiplication step for ΔfS
    ΔfT_list_temp::ArrayPartition       # temporary array forchange in distribution function due to TMatrix

    name_list::Vector{String}   # list of particle names
    nump_list::Vector{Int64}    # list of momentum bins for each particle
    numt_list::Vector{Int64}    # list of angular bins for each particle
    pu_list::Vector{Float32}    # list of upper momentum limits for each particle
    pl_list::Vector{Float32}    # list of lower momentum limits for each particle
    interaction_list_Binary::Vector{Vector{String}} # list of Binary interactions
    interaction_list_Sync::Vector{Vector{String}} # list of Sync interactions

    function BoltzmannEquation(f1D0,Lists)

        self = new()

        (self.name_list,self.nump_list,self.numt_list,self.pu_list,self.pl_list,self.interaction_list_Binary,self.interaction_list_Sync) = Lists

        self.t = Float32(0)
        #self.diff_coeff = diff_coeff

        # initialize distribution function vectors for indvidual species
        num_species = length(self.name_list)
        #self.f_list = Vector{Vector{Float32}}(undef,num_species)
        #for i in 1:num_species
            #self.f_list[i] = fill(Float32(0),self.nump_list[i]*self.numt_list[i])
        #end
        self.f1DA = fill!(similar(f1D0),Float32(0))
        self.f1DR = fill!(similar(f1D0),Float32(0))

        # initialize vectors for SMatrix and TMatrix changed so distribution functions for  indvidual species
        self.ΔfS_list = fill!(similar(f1D0),Float32(0))
        self.ΔfT_list = fill!(similar(f1D0),Float32(0))
        self.ΔfS_list_temp = fill!(similar(f1D0),Float32(0))
        self.ΔfT_list_temp = fill!(similar(f1D0),Float32(0))

        # setup temportray arrays for multiplication step in ΔfS
        temp = ();
        for i in eachindex(self.interaction_list_Binary)
            interaction = self.interaction_list_Binary[i]
            matricies = CollisionMatriciesBinary[interaction]
            name1 = interaction[1]
            name2 = interaction[2]
            name3 = interaction[3]
            name4 = interaction[4]

            if (name1 == name2) && (name3 == name4)
                temp = temp...,zeros(Float32,size(matricies[1],1))
            end
            if (name1 == name2) && (name3 != name4)
                temp = temp...,(zeros(Float32,size(matricies[1],1),size(matricies[2],1)))
            end
            if (name1 != name2) && (name3 == name4)
                temp = temp...,zeros(Float32,size(matricies[1],1))
            end
            if (name1 != name2) && (name3 != name4)
                temp = temp...,(zeros(Float32,size(matricies[1],1),size(matricies[2],1)))
            end

        end
        self.ΔfS_mul_step = ArrayPartition(temp)


        return self
    end

end