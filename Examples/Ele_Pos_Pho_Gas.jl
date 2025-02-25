# define constants with which to import data 

import BoltzmannCollisionIntegral as BCI
import BoltzmannEquationSolver as BES
using RecursiveArrayTools
using BenchmarkTools

# Set mode for simulation either "AXI" or "ISO"

    mode = "AXI"
    t_low = 0f0
    t_up = 1f3
    dt = t_up/500
    spacetime_coords = BES.Cylindrical()
    momentum_coords = BES.Spherical()

# vector of particle names to include in the simulation i.e. [name1, name2, ...]

    name_list::Vector{String} = ["Ele","Pos","Pho"];

# vector of particle momentum and angular bin numbers i.e. [p1_num, p2_num, ...], as well as upper and lower momentum values in log10 space

    p_up_list::Vector{Float64} = [4.0,4.0,4.0];
    p_low_list::Vector{Float64} = [-5.0,-5.0,-14.0];    
    p_grid_list::Vector{String} = ["l","l","l"];
    p_num_list::Vector{Int64} = [72,72,72];

    u_grid_list::Vector{String} = ["u","u","u"];
    u_num_list::Vector{Int64} = [8,8,8];
    

# tuple of vectors, each vector representing [1, 2, 3, 4] particles in the interactions to include in the simulation i.e. [[interaction1], [interaction2], ...]

    interaction_list_Binary::Vector{Vector{String}} = [["Ele", "Pos", "Pho", "Pho"],["Ele", "Pho", "Ele", "Pho"],["Pos", "Pho", "Pos", "Pho"],["Pho", "Pho", "Ele", "Pos"]];
    interaction_list_Emi::Vector{Vector{String}} = [];

# Collate lists into a single tuple to be passed to the solver

    Lists = BES.ListStruct(name_list,p_up_list,p_low_list,p_grid_list,p_num_list,u_grid_list,u_num_list,interaction_list_Binary,interaction_list_Emi);

    dimensions = 4
    forces::Vector{BES.ForceType} = [BES.CoordinateForce(),];
    SpaceTime = BES.SpaceTimeStruct(dimensions,spacetime_coords,momentum_coords,forces,t_up,t_low,dt);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd() * "/../BoltzmannCollisionIntegral/Data/"

# Load interaction matrices

    BigM = BES.BigMatrices(Lists);
    BES.LoadMatrices_Binary_Struct(BigM,DataDirectory,Lists;mode="AXI");
    BES.LoadMatrices_Emi_Struct(BigM,DataDirectory,Lists;mode=mode);

    FluxM = BES.FluxMatrices(Lists,SpaceTime)

# Set initial conditions

    # constant isotropic

        f_Ele = BES.Initial_Temperature(Lists,"Ele",1f10,1f0;mode=mode);
        #f_Ele = BES.Initial_Constant(Lists,"Ele",41,41,1,8,1f0 #=1cm^{-3}=#;mode=mode);
        f_Pos = BES.Initial_Temperature(Lists,"Pos",1f10,1f0;mode=mode);
        #f_Pos = BES.Initial_Constant(Lists,"Ele",41,41,1,8,1f0 #=1cm^{-3}=#;mode=mode);
        f_Pho = BES.Initial_Constant(Lists,"Ele",41,41,1,8,0f0 #=1cm^{-3}=#;mode=mode);

# Solver Setup and Time Stepping

    f_Init = ArrayPartition(f_Ele,f_Pos,f_Pho);

# Run BoltzmannEquationSolver

    scheme = BES.Euler(f_Init,Lists,SpaceTime,BigM,FluxM,false)

    @time sol = BES.Solve(f_Init,scheme;save_steps=1,progress=true);

# Plot results 

# Create vector of vectors for the distribution functions of each particle and momentum space grid values

    (mass_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list) = BES.getGridValues(Lists);

# Get scaling for plots (initial values)

    (numInit_list,engInit_list,tempInit_list) = BES.getInitialScaling(sol,Lists,pr_list,ur_list,mass_list);
    
# Plotting

    num_species = size(mass_list)[1];

    BES.AllPlots_Ani(sol,num_species,name_list,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,"test.mp4",24;istart=1,istop=50,iframe=1)

    fig_num = BES.NumPlot_AllSpecies(sol,num_species,pr_list,ur_list,p_num_list,u_num_list,mass_list,mode=mode)
    #BES.save("fig_num_3.png", fig_num)

    fig_eng = BES.EnergyPlot_AllSpecies(sol,num_species,pr_list,ur_list,p_num_list,u_num_list,mass_list,mode=mode)
    #BES.save("fig_eng_3.png", fig_eng)

    plot_dist_p = BES.PDistributionPlot_AllSpecies(sol,num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,save_dt,mode=mode,MaxwellJuttner=true,Flux=false)
    #BES.save("dist_p__sync_alpha_3.png", plot_dist_p)

    plot_dist_u = BES.uDistributionPlot_AllSpecies(sol,num_species,meanu_list,dp_list,du_list,p_num_list,u_num_list,save_dt)

    plot_dist_pu = BES.puDistributionPlot_AllSpecies(sol,num_species,meanu_list,meanp_list,dp_list,du_list,p_num_list,u_num_list,save_dt)

    fig_all = BES.AllPlots(sol,num_species,pr_list,ur_list,dp_list,du_list,meanp_list,meanu_list,numInit_list,engInit_list,tempInit_list,mass_list,p_num_list,u_num_list,save_dt,mode=mode)
    #BES.save("fig_all_two_temp.pdf", fig_all)
    #BES.save("fig_all_const.pdf", fig_all)
    #BES.save("fig_all_const_ani.pdf", fig_all)

    # animations
    BES.PDistributionPlot_AllSpecies_Ani("test_case_sync_p.mp4",sol,num_species,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,save_dt;mode="AXI")

    BES.uDistributionPlot_AllSpecies_Ani("test_case_4u.mp4",sol,num_species,meanμ_list,dp_list,p_num_list,u_num_list,save_dt;mode="AXI")
