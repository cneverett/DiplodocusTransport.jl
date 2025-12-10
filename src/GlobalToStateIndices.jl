"""
    GlobalIndices_To_StateIndex(x::Int64,y::Int64,z::Int64,px::Int64,py::Int64,pz::Int64,x_num::Int64,y_num::Int64,z_num::Int64,px_num::Int64,py_num::Int64,pz_num::Int64,off_name::Int64)

Returns a the index in the state vector corresponding to the global coordinate indices `x,y,z,px,py,pz` and `species_index`. 

"""
function GlobalIndices_To_StateIndex(x::Int64,y::Int64,z::Int64,px::Int64,py::Int64,pz::Int64,species_index::Int64,PhaseSpace::PhaseSpaceStruct)

    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num
    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    px_num = px_num_list[species_index]
    py_num = py_num_list[species_index]
    pz_num = pz_num_list[species_index]
    off_name = PhaseSpace.Grids.momentum_species_offset[species_index]

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    off_space::Int64 = (x-1)*y_num*z_num+(y-1)*z_num+z-1
    off_momentum::Int64 = (pz-1)*px_num*py_num+(py-1)*px_num+px

    return off_space*n_momentum + off_name + off_momentum 

end

"""
    Location_Species_To_StateVector(StateVector,PhaseSpace;x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1,species_index::Int64=0)

Returns a `view` to the section of the `StateVector` corresponding to the momentum space of `species` at x coordinate `x_idx`, y coordinate `y_idx` and z coordinate `z_idx`. 

"""
function Location_Species_To_StateVector(StateVector::Vector{F},PhaseSpace::PhaseSpaceStruct;x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1,species_index::Int64=0) where F<:AbstractFloat

    if iszero(species_index)
        error("Species not defined")
    end
    x_num = PhaseSpace.Space.x_num
    y_num = PhaseSpace.Space.y_num
    z_num = PhaseSpace.Space.z_num

    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    offset = PhaseSpace.Grids.momentum_species_offset

    px_num = px_num_list[species_index]
    py_num = py_num_list[species_index]
    pz_num = pz_num_list[species_index]

    off_space = (x_idx-1)*y_num*z_num+(y_idx-1)*z_num+z_idx-1
    off_name = offset[species_index]

    n_space = x_num+y_num+z_num
    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    start_idx = n_momentum*off_space+off_name+1
    end_idx   = n_momentum*off_space+off_name+px_num*py_num*pz_num

    return @view(StateVector[start_idx:end_idx])

end