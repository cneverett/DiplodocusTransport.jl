mutable struct FluxMatricesCoordinate <: Function

    Ap_Flux::Array{Float32}
    Am_Flux::Array{Float32}
    I_Flux::Array{Float32}
    J_Flux::Array{Float32}
    K_Flux::Array{Float32}

    function FluxMatricesCoordinate(Lists::ListStruct,SpaceTime::SpaceTimeStruct)
        self = new()

        (self.Ap_Flux,self.Am_Flux,self.I_Flux,self.J_Flux,self.K_Flux) = Allocate_Flux_Coordinate(Lists,SpaceTime)

        Build_Flux_Coordinate(self,Lists,SpaceTime)

        return self
    end
end

mutable struct FluxMatricesForce <: Function

    I_Flux::Array{Float32}
    J_Flux::Array{Float32}
    K_Flux::Array{Float32}

    function FluxMatricesForce(Lists::ListStruct,SpaceTime::SpaceTimeStruct)
        self = new()

        (self.I_Flux,self.J_Flux,self.K_Flux) = Allocate_Flux_Force(Lists,SpaceTime)

        Build_Flux_Force(self,Lists,SpaceTime)

        return self
    end
end



