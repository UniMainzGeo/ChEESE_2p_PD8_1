using LaMEM
using GeophysicalModelGenerator, GMT
using Plots

DirName="/local/home/iskander/projects/CHEESE2/etna_model/results/128x128x32_8CPU"
FileName="Volcano"

data, time = read_LaMEM_timestep(FileName, 40, DirName)

function get_topo_files(nx,ny,curdir,Generate_topo_files)

    if Generate_topo_files == 1
        # First time, you need to load topo:
        # Topo = import_topo([14.8,15.5,37.5,37.8], file="@earth_relief_01s")
        Topo          = import_topo([14.75,15.45,37.43,37.96], file="@earth_relief_01s")
        EGMS_velocity = load_GMG(joinpath(curdir,"data_EGMS_Etna_2019_2023"))       # INSAR data (to compare)

        proj      = ProjectionPoint(Lon=15.15, Lat=37.65)
        Topo_cart = convert2CartData(Topo, proj)
        EGMS_cart = convert2CartData(EGMS_velocity, proj)
        
        bounds   = [ minimum(Topo_cart.x.val),maximum(Topo_cart.x.val), minimum(Topo_cart.y.val),maximum(Topo_cart.y.val) ]
        x_coords = range(bounds[1],bounds[2],length=nx)
        y_coords = range(bounds[3],bounds[4],length=ny)
        Empty_cart_data = CartData(xyz_grid(x_coords,y_coords,0))

        Topo_model = project_CartData(Empty_cart_data, Topo, proj)
        EGMS_model = project_CartData(Empty_cart_data, EGMS_velocity, proj)
        
        EGMS_model = []
        save_GMG("Topo_model", Topo_model)
        save_GMG("EGMS_model", EGMS_model)

        Topo_model = load_GMG(joinpath(curdir,"Topo_model")) # Topography data from disk
        EGMS_model = load_GMG(joinpath(curdir,"EGMS_model")) # EGMS data from disk

        write_paraview(Topo_model,"Topo_cart")

    else
        EGMS_model = []
        Topo_model = load_GMG(joinpath(curdir,"Topo_model")) # Topography data from disk
        EGMS_model = load_GMG(joinpath(curdir,"EGMS_model")) # EGMS data from disk

    end

    
    # write_paraview(EGMS_model,"Topo_cart_EGMS")
    # Topo = drape_on_topo(Topo_model, EGMS_model)      # Drape the INSAR data on the topography

    return Topo_model,EGMS_model

end

function generate_setup_file(nx,ny,nz,Generate_topo_files,ParamFile_name,curdir,out_dir)

    Topo_model, EGMS_model = get_topo_files(nx,ny,curdir,Generate_topo_files)
    model = Model( 
                    Grid(               x               = [extrema(Topo_model.x.val)...],
                                        y               = [extrema(Topo_model.y.val)...],
                                        z               = [extrema(Topo_model.z.val)...],
                                        nel             = (nx,ny,nz)   ), 
                    
                    BoundaryConditions( temp_bot        = 1220.0,                     # we set temperature, but it is not used in this model
                                        temp_top        = 20.0,
                                        open_top_bound  = 1,                        # we do not want a freesurface, yet!
                                        noslip          = [0, 0, 0, 0, 0, 0]),      # [left, right, front, back, bottom, top]

                    # set timestepping parameters
                    Time(               time_end        = 10.0,                     # Time is always expressed in Myrs (input/output definition)
                                        dt              = 0.00001,                  # Target timestep, here 10k years
                                        dt_min          = 0.00000001,               # Minimum dt allowed, this is useful for more complex simulations
                                        dt_max          = 0.001,                    # max dt, here 1k years
                                        nstep_max       = 10,                       # Number of wanted timesteps
                                        nstep_out       = 5 ),                      # save output every nstep_out

                    # set solution parameters
                    SolutionParams(     eta_min         = 1e17,
                                        eta_ref         = 1e20,
                                        eta_max         = 1e23),
                    ModelSetup(
                                        msetup           =  "files",
                                        rand_noise       =  1,
                                        bg_phase         =  0,
                                        advect           =  "basic",
                                        interp           =  "stag",
                                        mark_ctrl        =  "subgrid",
                                        nmark_sub        =  2 ),

                # what will be saved in the output of the simulation

                Output(             out_file_name       = "Volcano",
                                    param_file_name     = ParamFile_name,
                                    write_VTK_setup     = false,
                                    out_dir             = out_dir,
                                    out_tot_displ       = 1,
                                    out_surf            = 1, 	
                                    out_surf_pvd        = 1,
                                    out_surf_topography = 1,
                                    out_surf_velocity   = 1,
                                    out_moment_res      = 1,
                                    out_cont_res        = 1,
                                    out_temperature     = 1,
                                    out_yield           = 1 ),) 

    #set air and water properties
    air = set_air(alpha=3e-5, G=3e10, nu=0.2, 
                    #ch=10e6, 
                    #fr=30
                    )
    #air.eta_vp = 1e21
    water=copy_phase(air, rho=1000.0, Name="water", ID=1 );

    crust      = Phase(Name="Elastic Crust", 
                    ID=2, rho=2900, 
                    alpha=3e-5, 
                    eta = 1e23,
                    #disl_prof="Mafic_Granulite-Ranalli_1995",
                    G=3e10, 
                    nu=0.2, 
                    k=3, 
                    Cp=1000, 
                    )

    crust_alt  = copy_phase(crust  , ID=3, Name="AlteredCrust")

    # to be activated with a phase transition                
    ϕ        = 0.0;
    cohesion = 3e6
    ηvp      = 1e18

    air_plastic          = copy_phase(air  ,     ID=4, Name="AirPlastic",            fr=ϕ,       ch=cohesion,    eta_vp=ηvp)
    water_plastic        = copy_phase(water,     ID=5, Name="WaterPlastic",          fr=ϕ,       ch=cohesion,    eta_vp=ηvp)
    water                = copy_phase(water_plastic, ID=1)
    air                  = copy_phase(air_plastic, ID=0)

    crust_plastic        = copy_phase(crust,     ID=6, Name="CrustPlastic",          fr=30.0,    ch=400e6,       eta_vp=ηvp)
    #crust_plastic       = copy_phase(crust,     ID=6, Name="CrustPlastic",          fr=ϕ,    ch=cohesion,       eta_vp=ηvp)
    crust_alt_plastic    = copy_phase(crust_alt, ID=7, Name="CrustAlteredPlastic",   fr=ϕ,       ch=cohesion,     eta_vp=ηvp)

    rm_phase!(model)
    add_phase!(model, air, water, crust, crust_alt, air_plastic, water_plastic, crust_plastic,  crust_alt_plastic)

    # add phase transitions to activate plasticity after an initial stress state was established (after a certain time)
    phase_transition_crust           = PhaseTransition(ID=0, 
                                                Type="Constant", 
                                                Parameter_transition="t", 
                                                ConstantValue = 4e-5,                # in Myrs
                                                PhaseBelow=[0; 1; 2; 3], 
                                                PhaseAbove=[4; 5; 6; 7], 
                                                number_phases=4, 
                                                PhaseDirection="BelowToAbove")
    model.Materials.PhaseTransitions = [phase_transition_crust]

    # write_LaMEM_inputFile(model, joinpath(curdir,model.Output.param_file_name))

    return model

end

curdir="/local/home/iskander/projects/CHEESE2/etna_model/postproc"
out_dir = ""
cd(curdir)

ParamFile_name = "Volcano.dat"



Generate_topo_files     = 1;          # Generate jld files for topography and velocity model
n_ranks                 = 8;         # Number of processors
RandomNoise             = false ;     # add random noise to particles, does not work
verbose                 = true;       # print info
num_prop                = 5;          # Number of properties, do not change
directory               = "markers" # do not change

# Resolution
nx = size(data)[1] - 1 # number of points in x direction
ny = size(data)[2] - 1 # number of points in y direction
nz = size(data)[3] - 1 # number of points in z direction

Topo_model, EGMS_model = get_topo_files(nx,ny,curdir,Generate_topo_files)

model = generate_setup_file(nx,ny,nz,Generate_topo_files,ParamFile_name,curdir,out_dir)

plot_cross_section(model, y=0, field=:phase)