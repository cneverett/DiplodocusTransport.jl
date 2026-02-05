"""
    CharacteristicStruct()

A struct for storing the characteristic scales (in SI units) for the simulation, from which all other values are normalised to.
"""
@kwdef struct CharacteristicStruct

    # fundamental scales (SI units)
    CHAR_mass::Float64 = CONST_mEle
    CHAR_charge::Float64 = CONST_q
    CHAR_speed::Float64 = CONST_c
    CHAR_magnetic_field::Float64 = 1.0 # 1 Tesla

    # user defined scales (SI units)
    CHAR_length::Float64 = try 
            getfield(Main,Symbol("CHAR_length"))
        catch
            CONST_pcs
        end
    CHAR_number_density::Float64 = try 
            getfield(Main,Symbol("CHAR_number_density"))
        catch
            1.0e6 # 1 particle per cm^3
        end

    # derived scales (SI units)
    CHAR_time::Float64 = CHAR_length/CHAR_speed
    CHAR_momentum::Float64 = CHAR_mass*CHAR_speed
    CHAR_energy::Float64 = CHAR_mass*CHAR_speed^2
    CHAR_mass_density::Float64 = CHAR_mass*CHAR_number_density
    CHAR_energy_density::Float64 = CHAR_mass_density*CHAR_speed^2
    CHAR_pressure::Float64 = CHAR_energy_density
    CHAR_electric_field::Float64 = CHAR_magnetic_field*CONST_c

    # scales used in DiplodocusCollisions and for normalisation of Big and Flux Matrices
    Bin_Norm::Float64 = CONST_σT * CONST_c * CHAR_time * CHAR_number_density
    Emi_Norm::Float64 = CONST_σT * CONST_c * CHAR_time
    Flux_Norm::Float64 = CONST_c * CHAR_time / CHAR_length
    
end


"""
    TimeStruct()

A struct for storing the time domain of the simulation.
"""
@kwdef struct TimeStruct

    t_up::Float64    = getfield(Main,Symbol("t_up"))
    t_low::Float64   = getfield(Main,Symbol("t_low"))
    t_num::Int64     = getfield(Main,Symbol("t_num"))
    t_grid::String   = getfield(Main,Symbol("t_grid"))

end

"""
    SpaceStruct()

A struct for storing the space domain of the simulation.
"""
@kwdef struct SpaceStruct

    space_coordinates::CoordinateType = getfield(Main,Symbol("space_coords"))
    scheme::String = getfield(Main,Symbol("space_scheme"))

    x_up::Float64   = getfield(Main,Symbol("x_up"))     
    x_low::Float64  = getfield(Main,Symbol("x_low"))
    x_grid::String  = getfield(Main,Symbol("x_grid"))
    x_num::Int64    = getfield(Main,Symbol("x_num"))

    y_up::Float64   = getfield(Main,Symbol("y_up"))
    y_low::Float64  = getfield(Main,Symbol("y_low"))
    y_grid::String  = getfield(Main,Symbol("y_grid"))
    y_num::Int64    = getfield(Main,Symbol("y_num"))

    z_up::Float64   = getfield(Main,Symbol("z_up")) 
    z_low::Float64  = getfield(Main,Symbol("z_low"))
    z_grid::String  = getfield(Main,Symbol("z_grid"))
    z_num::Int64    = getfield(Main,Symbol("z_num"))

end

"""
    MomentumStruct()

A struct for storing the momentum domain of the simulation.
"""
@kwdef struct MomentumStruct 

    momentum_coordinates::CoordinateType = getfield(Main,Symbol("momentum_coords"))
    scheme::String = getfield(Main,Symbol("momentum_scheme"))

    px_up_list::Vector{Float64}     = getfield(Main,Symbol("px_up_list"))
    px_low_list::Vector{Float64}    = getfield(Main,Symbol("px_low_list"))
    px_grid_list::Vector{String}    = getfield(Main,Symbol("px_grid_list"))
    px_num_list::Vector{Int64}      = getfield(Main,Symbol("px_num_list"))

    py_up_list::Vector{Float64}     = getfield(Main,Symbol("py_up_list"))
    py_low_list::Vector{Float64}    = getfield(Main,Symbol("py_low_list"))
    py_grid_list::Vector{String}    = getfield(Main,Symbol("py_grid_list"))
    py_num_list::Vector{Int64}      = getfield(Main,Symbol("py_num_list"))

    pz_up_list::Vector{Float64}     = getfield(Main,Symbol("pz_up_list"))
    pz_low_list::Vector{Float64}    = getfield(Main,Symbol("pz_low_list"))
    pz_grid_list::Vector{String}    = getfield(Main,Symbol("pz_grid_list"))
    pz_num_list::Vector{Int64}      = getfield(Main,Symbol("pz_num_list"))

end

"""
    GridsStruct()

A struct for storing the grid values for each particle in the simulation.
"""
@kwdef mutable struct GridsStruct <: Function

    mass_list::Vector{Float64}
    charge_list::Vector{Float64}

    tr::Vector{Float64}
    dt::Vector{Float64}

    xr::Vector{Float64}
    yr::Vector{Float64}
    zr::Vector{Float64}

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    mx::Vector{Float64}
    my::Vector{Float64}
    mz::Vector{Float64}

    pxr_list::Vector{Vector{Float64}}
    pyr_list::Vector{Vector{Float64}}
    pzr_list::Vector{Vector{Float64}}

    dpx_list::Vector{Vector{Float64}}
    dpy_list::Vector{Vector{Float64}}
    dpz_list::Vector{Vector{Float64}}

    mpx_list::Vector{Vector{Float64}}
    mpy_list::Vector{Vector{Float64}}
    mpz_list::Vector{Vector{Float64}}

    dE_list::Vector{Vector{Float64}}

    momentum_species_offset::Vector{Int64}

    B_field::Array{Float64,3} # magnetic field values on the x,y,z grid
    E_field::Array{Float64,3} # electric field values on the x,y,z grid

    function GridsStruct(name_list,time::TimeStruct,space::SpaceStruct,momentum::MomentumStruct,characteristic::CharacteristicStruct,ElectromagneticField::Union{ElectroMagneticFieldType,Nothing})

        self = new()

        # time domain grids

            t_up    = time.t_up
            t_low   = time.t_low
            t_num   = time.t_num
            t_grid  = time.t_grid

            self.tr = bounds(t_low,t_up,t_num,t_grid)
            self.dt = deltaVector(self.tr)

        # space domain grids

            x_up    = space.x_up
            x_low   = space.x_low
            x_grid  = space.x_grid
            x_num   = space.x_num

            y_up    = space.y_up
            y_low   = space.y_low
            y_grid  = space.y_grid
            y_num   = space.y_num

            z_up    = space.z_up
            z_low   = space.z_low
            z_grid  = space.z_grid
            z_num   = space.z_num

            self.xr = bounds(x_low,x_up,x_num,x_grid)
            self.yr = bounds(y_low,y_up,y_num,y_grid)
            self.zr = bounds(z_low,z_up,z_num,z_grid)

            self.dx = deltaVector(self.xr)
            self.dy = deltaVector(self.yr)
            self.dz = deltaVector(self.zr)

            self.mx = meanVector(self.xr)
            self.my = meanVector(self.yr)
            self.mz = meanVector(self.zr)

        # momentum domain grids

            px_up_list      = momentum.px_up_list
            px_low_list     = momentum.px_low_list
            px_grid_list    = momentum.px_grid_list
            px_num_list     = momentum.px_num_list

            py_up_list      = momentum.py_up_list
            py_low_list     = momentum.py_low_list
            py_grid_list    = momentum.py_grid_list
            py_num_list     = momentum.py_num_list

            pz_up_list      = momentum.pz_up_list
            pz_low_list     = momentum.pz_low_list
            pz_grid_list    = momentum.pz_grid_list
            pz_num_list     = momentum.pz_num_list

            num_species     = length(name_list);

            self.mass_list  = Vector{Float64}(undef,num_species)
            self.charge_list=Vector{Float64}(undef,num_species)

            self.pxr_list   = Vector{Vector{Float64}}(undef,num_species);
            self.pyr_list   = Vector{Vector{Float64}}(undef,num_species);
            self.pzr_list   = Vector{Vector{Float64}}(undef,num_species);

            self.dpx_list   = Vector{Vector{Float64}}(undef,num_species);
            self.dpy_list   = Vector{Vector{Float64}}(undef,num_species);
            self.dpz_list   = Vector{Vector{Float64}}(undef,num_species);

            self.mpx_list   = Vector{Vector{Float64}}(undef,num_species);
            self.mpy_list   = Vector{Vector{Float64}}(undef,num_species);
            self.mpz_list   = Vector{Vector{Float64}}(undef,num_species);

            self.dE_list    = Vector{Vector{Float64}}(undef,num_species);

            self.momentum_species_offset = Vector{Int64}(undef,num_species);

            for i in eachindex(name_list)

                self.mass_list[i]   = eval(Symbol("CONST_mu"*name_list[i]));
                self.charge_list[i] = eval(Symbol("CONST_z"*name_list[i]));

                self.pxr_list[i]    = bounds(px_low_list[i],px_up_list[i],px_num_list[i],px_grid_list[i]);
                self.pyr_list[i]    = bounds(py_low_list[i],py_up_list[i],py_num_list[i],py_grid_list[i]);
                self.pzr_list[i]    = bounds(pz_low_list[i],pz_up_list[i],pz_num_list[i],pz_grid_list[i]);

                self.dpx_list[i]    = deltaVector(self.pxr_list[i]);
                self.dpy_list[i]    = deltaVector(self.pyr_list[i]);
                self.dpz_list[i]    = deltaVector(self.pzr_list[i]);
                self.mpx_list[i]    = meanVector(self.pxr_list[i]);
                self.mpy_list[i]    = meanVector(self.pyr_list[i]);
                self.mpz_list[i]    = meanVector(self.pzr_list[i]);

                self.dE_list[i]     = deltaEVector(self.pxr_list[i],self.mass_list[i]) ./ self.dpx_list[i];

                if i == 1
                    self.momentum_species_offset[i] = 0
                else
                    self.momentum_species_offset[i] = self.momentum_species_offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
                end

            end

            if !isnothing(ElectromagneticField)
                # build electromagnetic field grids
                parameters = ElectromagneticField.parameters
                ElectromagneticFieldFunction = ElectromagneticField.ElectromagneticFieldFunction
                self.B_field, self.E_field = ElectromagneticFieldFunction(space,momentum,characteristic,self,parameters)
            else
                self.B_field = zeros(Float64,0,0,0)
                self.E_field = zeros(Float64,0,0,0)
            end

        return self

    end

end

"""
    PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,forces)

A struct for storing the phase space of the simulation.
"""
@kwdef struct PhaseSpaceStruct

    # list of particle names
    name_list::Vector{String}           = getfield(Main,Symbol("name_list"))
    Characteristic::CharacteristicStruct= CharacteristicStruct()
    Time::TimeStruct                    = TimeStruct()    
    Space::SpaceStruct                  = SpaceStruct()
    Momentum::MomentumStruct            = MomentumStruct()
    Forces::Vector{ForceType}           = getfield(Main,Symbol("Forces"))

    # interactions
    Binary_list::Vector{BinaryStruct}   = getfield(Main,Symbol("Binary_list")) 
    Emi_list::Vector{EmiStruct}         = getfield(Main,Symbol("Emi_list"))

    ElectromagneticField::Union{ElectroMagneticFieldType,Nothing} = try 
            getfield(Main,Symbol("ElectromagneticField"))
        catch
            nothing 
        end

    # grids
    Grids::GridsStruct  = GridsStruct(name_list,Time,Space,Momentum,Characteristic,ElectromagneticField)

end


"""
    BinaryMatricesStruct()

A struct for storing the big matrices associated with binary interactions in the simulation.
"""
struct BinaryMatricesStruct{T<:Union{Float32,Float64}}
    
    M_Bin::AbstractMatrix{T}  # big matrix for binary interactions
    Binary_list::Vector{BinaryStruct} # list of binary interactions
    Domain::Vector{Int64} # domain of the binary interaction in the state vector

end

"""
    EmissionMatricesStruct()

A struct for storing the big matrices associated with emission interactions and their corresponding reaction forces if applicable in the simulation.
"""
struct EmissionMatricesStruct{T<:Union{Float32,Float64}}
    
    M_Emi::AbstractMatrix{T}  # big matrix for emission interactions
    Emission_list::Vector{EmiStruct} # list of emission interactions

end

"""
    FluxMatricesStruct()

A struct for storing the flux matrices associated with the simulation.
"""
struct FluxMatricesStruct{T<:Union{Float32,Float64}}

    # time fluxes
    Ap_Flux::Vector{T}
    Am_Flux::Vector{T}
    # sum of space and momentum fluxes for speed 
    F_Flux::SparseMatrixCSC{T,Int64}
    # space time volume element vector
    Vol::Vector{T}
    # space fluxes
    B_Flux::SparseMatrixCSC{T,Int64}
    C_Flux::SparseMatrixCSC{T,Int64}
    D_Flux::SparseMatrixCSC{T,Int64} 
    # momentum fluxes
    I_Flux::SparseMatrixCSC{T,Int64}
    J_Flux::SparseMatrixCSC{T,Int64}
    K_Flux::SparseMatrixCSC{T,Int64}

end


mutable struct OutputStruct
    
    f::AbstractVector
    f_ext::AbstractVector
    t::Vector{Float64}

    function OutputStruct(f0::AbstractVector,n_save::Int64)

        self = new()
        self.f = Vector{typeof(f0)}(undef,n_save)
        self.f_ext = Vector{typeof(f0)}(undef,n_save)
        self.t = Vector{Float64}(undef,n_save)

        return self
    
    end

end
