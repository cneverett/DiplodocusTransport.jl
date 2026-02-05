"""
    Constant_ElectroMagneticField(PhaseSpace, x_index,y_index,z_index)

Returns the values of the magnetic and electric fields at a given spatial grid point (x_index,y_index,z_index) for a constant electromagnetic field with parameters defined by the `Constant_ElectroMagneticField` struct. For a constant magnetic field in the z-direction, this function returns the magnetic field strength B and zero electric field.

"""
function Constant_ElectroMagneticField(Space::SpaceStruct,Momentum::MomentumStruct,Characteristic::CharacteristicStruct,Grids::GridStruct,parameters::Vector{Float64,2}) 

    space_coords = Space.space_coords
    momentum_coords = Momentum.momentum_coords
    Characteristic = Characteristic
    Grids = Grids

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num

    L = Characteristic.CHAR_length
    B0 = parameters[1]
    E0 = parameters[2]

    B_field = zeros(Float64,x_num,y_num,z_num)
    E_field = zeros(Float64,x_num,y_num,z_num)

    # B field is ALWAYS in local momentum z-direction and E field in local momentum y-direction using local orthonormal basis (n^α,ϵ^αβγδn_βE_γB_δ,E^α,B^α)

    if space_coords isa Cartesian # B along z, # E along y
        for ix in 1:x_num, iy in 1:y_num, iz in 1:z_num
            B_field[ix,iy,iz] = B0
            E_field[ix,iy,iz] = E0
        end
    end

end