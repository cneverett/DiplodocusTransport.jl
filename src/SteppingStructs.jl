# Struct for storing the Boltzmann equation and its solution
mutable struct Euler <: SteppingMethod

    t::Float32                  # the last time step time to calculate Δt
    dt::Float32                 # time step

    Lists::ListStruct
    BigM::BigMatrices
    FluxMC::FluxMatricesCoordinate
    FluxMF::FluxMatricesForce
    SpaceTime::SpaceTimeStruct

    Implicit::Bool

    A_Binary_Mul::Array{Float32,2}  # temporary array for matrix multiplication of binary terms
    Jac::Array{Float32,2}           # jacobian for if Implicit==True
    df::fType                       # change in distribution function
    df_temp::fType                  # change in distribution function

    function Euler(f0,Lists::ListStruct,SpaceTime::SpaceTimeStruct,Big_Matrices::BigMatrices,Flux_Matrices_Coord::FluxMatricesCoordinate,Flux_Matrices_Force::FluxMatricesForce,Implicit::Bool)

        self = new()

        self.Lists = Lists
        self.BigM = Big_Matrices
        self.FluxMC = Flux_Matrices_Coord
        self.FluxMF = Flux_Matrices_Force
        self.SpaceTime = SpaceTime

        self.Implicit = Implicit  

        self.A_Binary_Mul = zeros(Float32,length(f0),length(f0))
        if Implicit
            self.Jac = zeros(Float32,length(f0),length(f0))
        end
        self.df = fill!(similar(f0),Float32(0))
        self.df_temp = fill!(similar(f0),Float32(0))

        return self
    end

end