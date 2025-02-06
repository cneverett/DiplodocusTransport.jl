# Empty dictionary for storing binary collision matrices by interaction name
Matrices_BinaryInteraction = Dict{Vector{String},Tuple}()
Matrices_Synchrotron = Dict{Vector{String},Array{Float32,2}}()
Matrices_Force = Dict{Vector{String},Array{Float32,2}}()

mutable struct BigMatrices <: Function
    
    A_Binary::Array{Float32,2}    # big matrix for binary interactions

    A_Emi::Array{Float32,2}  # big matrix for emission interactions
    J_Emi::Array{Float32,2}  # big matrix for Jacobian of emission interactions

    A_Abs::Array{Float32,2}  # big matrix for emission interactions
    J_Abs::Array{Float32,2}  # big matrix for Jacobian of emission interactions

    function BigMatrices(Lists::Tuple)
        self = new()

        self.A_Binary = Allocate_A_Binary(Lists)
        #self.A_Emi = Allocate_A_Emi(Lists)
        #self.J_Emi = Allocate_J_Emi(Lists)
        #self.A_Abs = Allocate_A_Abs(Lists)
        #self.J_Abs = Allocate_J_Abs(Lists)

        return self
    end

end

# Struct for storing the Boltzmann equation and its solution
mutable struct BoltzmannEquation <: Function

    t::Float32                  # the last time step time to calculate Δt
    dt::Float32                 # time step

    #f_list::Vector{Vector{Float32}} # vector of distribution functions for each particle
    f1DA::ArrayPartition  # advanced distribution function 
    f1DR::ArrayPartition  # retarded distribution function
    state::Bool

    ΔfS_list::ArrayPartition       # change in distribution function due to SMatrix
    ΔfT_list::ArrayPartition       # change in distribution function due to TMatrix
    ΔfS_list_temp::ArrayPartition       # temporary array for change in distribution function due to SMatrix
    ΔfS_mul_step::ArrayPartition       # temporary array the matrix multiplication step for ΔfS
    ΔfT_list_temp::ArrayPartition       # temporary array for change in distribution function due to TMatrix

    name_list::Vector{String}   # list of particle names
    p_grid_list::Vector{String}    # list of momentum grid types for each particle
    p_num_list::Vector{Int64}    # list of momentum bins for each particle
    u_grid_list::Vector{String}    # list of angular grid types for each particle
    u_num_list::Vector{Int64}    # list of angular bins for each particle
    p_up_list::Vector{Float32}    # list of upper momentum limits for each particle
    p_low_list::Vector{Float32}    # list of lower momentum limits for each particle
    interaction_list_Binary::Vector{Vector{String}} # list of Binary interactions
    interaction_list_Emi::Vector{Vector{String}} # list of Emission interactions

    A_Binary_Reshape::Array{Float32,1}    # reshaped big matrix for binary interactions
    Δf::ArrayPartition              # change in distribution function
    Δf_temp::ArrayPartition         # change in distribution function
    J::Array{Float32,2}             # Jacobian matrix

    BigM::BigMatrices

    function BoltzmannEquation(f1D0,Lists,Big_Matrices::BigMatrices,dt)

        self = new()

        (self.name_list,self.p_up_list,self.p_low_list,self.p_grid_list,self.p_num_list,self.u_grid_list,self.u_num_list,self.interaction_list_Binary,self.interaction_list_Emi) = Lists

        self.t = Float32(0)
        self.dt = dt

        # initialize distribution function vectors for individual species
        num_species = length(self.name_list)
        #self.f_list = Vector{Vector{Float32}}(undef,num_species)
        #for i in 1:num_species
            #self.f_list[i] = fill(Float32(0),self.nump_list[i]*self.numt_list[i])
        #end
        self.f1DA = fill!(similar(f1D0),Float32(0))
        self.f1DR = fill!(similar(f1D0),Float32(0))

        # initialize vectors for SMatrix and TMatrix changed so distribution functions for  individual species
        self.ΔfS_list = fill!(similar(f1D0),Float32(0))
        self.ΔfT_list = fill!(similar(f1D0),Float32(0))
        self.ΔfS_list_temp = fill!(similar(f1D0),Float32(0))
        self.ΔfT_list_temp = fill!(similar(f1D0),Float32(0))

        self.BigM = Big_Matrices

        self.A_Binary_Reshape = zeros(Float32,size(Big_Matrices.A_Binary,1))
        self.Δf = fill!(similar(f1D0),Float32(0))
        self.Δf_temp = fill!(similar(f1D0),Float32(0))
        self.J = zeros(Float32,size(Big_Matrices.A_Binary,2),size(Big_Matrices.A_Binary,2))

        # setup temporary arrays for multiplication step in ΔfS
#=         temp = ();
        for i in eachindex(self.interaction_list_Binary)
            interaction = self.interaction_list_Binary[i]
            matrices = Matrices_BinaryInteraction[interaction]
            name1 = interaction[1]
            name2 = interaction[2]
            name3 = interaction[3]
            name4 = interaction[4]

            if (name1 == name2) && (name3 == name4)
                temp = temp...,zeros(Float32,size(matrices[1],1))
            end
            if (name1 == name2) && (name3 != name4)
                temp = temp...,(zeros(Float32,size(matrices[1],1)),zeros(Float32,size(matrices[2],1)))
            end
            if (name1 != name2) && (name3 == name4)
                temp = temp...,zeros(Float32,size(matrices[1],1))
            end
            if (name1 != name2) && (name3 != name4)
                temp = temp...,(zeros(Float32,size(matrices[1],1)),zeros(Float32,size(matrices[2],1)))
            end

        end
        self.ΔfS_mul_step = ArrayPartition(temp)
 =#
        return self
    end

end

