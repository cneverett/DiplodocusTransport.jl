# define constants with which to import data 

import BoltzmannCollisionIntegral as BCI
import BoltzmannEquationSolver as BES
using RecursiveArrayTools
using BenchmarkTools

# Set mode for simulation either "AXI" or "ISO"

    mode = "AXI"
    t_low = 0f0
    t_up = 1f5
    dt = t_up/5000
    spacetime_coords = BES.Cylindrical()
    momentum_coords = BES.Spherical()

# vector of particle names to include in the simulation i.e. [name1, name2, ...]

    name_list::Vector{String} = ["Ele",];

# vector of particle momentum and angular bin numbers i.e. [p_num1, p_num2, ...], as well as upper and lower momentum values in log10 space

    # electron synchrotron
    p_up_list::Vector{Float64} = [4.0,];
    p_low_list::Vector{Float64} = [-5.0,];
    p_grid_list::Vector{String} = ["l",];
    p_num_list::Vector{Int64} = [72,];
    u_grid_list::Vector{String} = ["u",];
    u_num_list::Vector{Int64} = [8,];

# tuple of vectors, each vector representing [1, 2, 3, 4] particles in the interactions to include in the simulation i.e. [[interaction1], [interaction2], ...]

    interaction_list_Binary::Vector{Vector{String}} = [];
    interaction_list_Emi::Vector{Vector{String}} = [];

# Collate lists into a single tuple to be passed to the solver

    Lists = BES.ListStruct(name_list,p_up_list,p_low_list,p_grid_list,p_num_list,u_grid_list,u_num_list,interaction_list_Binary,interaction_list_Emi);

    dimensions = 4
    forces::Vector{BES.ForceType} = [BES.CoordinateForce(),BES.SyncRadReact(),];
    SpaceTime = BES.SpaceTimeStruct(dimensions,spacetime_coords,momentum_coords,forces,t_up,t_low,dt);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd() * "/../BoltzmannCollisionIntegral/Data/"

# Load interaction matrices

    BigM = BES.BigMatrices(Lists);
    BES.LoadMatrices_Binary_Struct(BigM,DataDirectory,Lists;mode="AXI");

    FluxM = BES.FluxMatrices(Lists,SpaceTime)

    BES.LoadMatrices_Emi_Struct(BigM,DataDirectory,Lists;mode=mode);

# Set initial conditions

    # Power law distribution of electrons
    f1D0_Ele =  BES.Initial_PowerLaw(Lists,"Ele",40,70,1,8,2f0,1f0);

# Solver Setup and Time Stepping

    f1D0 = ArrayPartition(f1D0_Ele,);

# Run BoltzmannEquationSolver

    scheme = BES.Euler(f1D0,Lists,SpaceTime,BigM,FluxM,false)

    @time sol = BES.Solve(f1D0,scheme;save_steps=1,progress=true);

# Plot results 

# Create vector of vectors for the distribution functions of each particle and momentum space grid values

    (mass_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list) = BES.getGridValues(Lists);

# Get scaling for plots (initial values)

    (numInit_list,engInit_list,tempInit_list) = BES.getInitialScaling(sol,Lists,pr_list,ur_list,mass_list);
    
# Plotting

    num_species = length(name_list);

    BES.AllPlots_Ani(sol,num_species,name_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,"EleRadReact.mp4";istart=1,istop=50,iframe=nothing)