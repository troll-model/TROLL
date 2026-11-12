#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"
#include "internal/lookup_tables_cache.hpp"
#include "internal/param_registry.hpp"

// ###########################################
// ###########################################
// ###########     MAIN PROGRAM    ###########
// ###########################################
// ###########################################

int main(int argc, char *argv[])
{
    Context ctx{};

    //!*********************
    //!** Initializations **
    //!*********************

    for (int argn = 1; argn < argc; argn++)
    { // Arguments of the input and output files
        if (*argv[argn] == '-')
        {
            switch (*(argv[argn] + 1))
            {
            case 'i':
                ctx.buffers.bufi = argv[argn] + 2;
                break;
            case 'd': // new v.3.0; initialisation of daytime variation in a separate file ('d' for daytime variation)
                ctx.buffers.bufi_daytimevar = argv[argn] + 2;
                break;
            case 'm': // new v.2.4; initialisation of climate parameters from separate file ('m' for meterology)
                ctx.buffers.bufi_climate = argv[argn] + 2;
                break;
            case 'p': // new v.3.0; initialisation of soil parameters from separate file ('p' for pedology)
                ctx.buffers.bufi_soil = argv[argn] + 2;
                break;
            case 's': // new v.3.0; initialisation of species parameters from separate file
                ctx.buffers.bufi_species = argv[argn] + 2;
                break;
            case 'o':
                ctx.buffers.buf = argv[argn] + 2;
                break;
            case 'f': // new v.2.3: initialisation from field, 'f' for "forest", "field data"
                ctx.buffers.bufi_data = argv[argn] + 2;
                ctx.opt._FromInventory = 1; // new v.3.1: automatic recognition of whether data sheet is provided or not
                break;
            case 'w': // new v.2.3: initialisation from field, 'f' for "forest", "field data"
                ctx.buffers.bufi_dataSWC = argv[argn] + 2;
                break;
            case 'l':
                ctx.buffers.bufi_pointcloud = argv[argn] + 2; // new v.3.1.6: output of simulated point clouds for TROLL-created stands
                ctx.opt._OUTPUT_pointcloud = 1;
                break;
            }
        }
    }

    // input files
    sprintf(ctx.fileio.inputfile, "%s", ctx.buffers.bufi);
    sprintf(ctx.fileio.inputfile_daytimevar, "%s", ctx.buffers.bufi_daytimevar);
    sprintf(ctx.fileio.inputfile_climate, "%s", ctx.buffers.bufi_climate);
    sprintf(ctx.fileio.inputfile_soil, "%s", ctx.buffers.bufi_soil);
    sprintf(ctx.fileio.inputfile_species, "%s", ctx.buffers.bufi_species);

    if (ctx.opt._OUTPUT_pointcloud)
    {
        sprintf(ctx.fileio.inputfile_pointcloud, "%s", ctx.buffers.bufi_pointcloud); // v.3.1.6
    }

    if (ctx.opt._FromInventory)
    {
        sprintf(ctx.fileio.inputfile_inventory, "%s", ctx.buffers.bufi_data);
#ifdef WATER
        sprintf(ctx.fileio.inputfile_SWC, "%s", ctx.buffers.bufi_dataSWC);
#endif
    }
    // Registers parameters (used by AssignParamFromRegistry, the AssignValueGlobal replacement)
    // see src/params/param_registry.hpp and src/params/param_registry.cpp
    RegisterParameters(ctx);
    // v.3.1: removed par output, because no single parameter sheet provided anymore (in future all separate parameter sheets could be provided as outputs as well
    ReadInputGeneral(ctx); // v.3.1 has to be done before initialisation of random number generators (ctx.opt._NONRANDOM)

    // Stuff for constant number generator
    const gsl_rng_type *Trandgsl;
    gsl_rng_env_setup();
    Trandgsl = gsl_rng_default;
    ctx.rng.gslrand = gsl_rng_alloc(Trandgsl);

    unsigned long int seed = (unsigned long int)time(NULL);

    if (ctx.opt._NONRANDOM == 1)
        seed = 1;

    gsl_rng_set(ctx.rng.gslrand, seed);

    cout << "RNG seed: " << seed << endl;
    sprintf(ctx.fileio.outputinfo, "%s_info.txt", ctx.buffers.buf);
    ctx.out.output_info.open(ctx.fileio.outputinfo, ios::out);
    if (!ctx.out.output_info)
        cerr << "ERROR with info file" << endl;

    Initialise(ctx); // Read global parameters

    InitialiseOutputStreams(ctx); // Initialise Output streams, taken outside of Initialise() function in v.3.1 to mirror AllocMem()
    AllocMem(ctx);                // Memory allocation

#ifdef Output_ABC
    InitialiseABC(ctx);
#endif

    if (ctx.opt._OUTPUT_pointcloud)
    {
        ReadInputPointcloud(ctx); // parameters for point cloud generation, v.3.1.6
    }

    if (ctx.opt._FromInventory)
    {
        ReadInputInventory(ctx); // Initial configuration of the forest, read from data
    }

    cout << "klight is: " << ctx.params.klight << endl;
    cout << "CO2 concentration is: " << ctx.params.Cair << endl;
    cout << "Number of species: " << ctx.grid.nbspp << endl
         << endl;
#ifdef WATER
    cout << "Atmospheric pressure is: " << ctx.params.PRESS << endl;
#endif

    if (ctx.opt._GPPcrown == 1)
        cout << "Activated Module: FastGPP" << endl;
    if (ctx.opt._BASICTREEFALL == 1)
        cout << "Activated Module: BASICTREEFALL" << endl;
    if (ctx.opt._NDD == 1)
        cout << "Activated Module: NDD" << endl;
    if (ctx.opt._SEEDTRADEOFF == 1)
        cout << "Activated Module: SEEDTRADEOFF" << endl;
    if (ctx.opt._FromInventory == 1)
        cout << "Activated Module: FromInventory" << endl;
    if (ctx.opt._OUTPUT_extended == 1)
        cout << "Activated Module: OUTPUT_extended" << endl;
    if (ctx.opt._OUTPUT_inventory == 1)
        cout << "Activated Module: OUTPUT_inventory" << endl;
    if (ctx.opt._OUTPUT_extended == 1 && ctx.crown.extent_visual > 0)
        cout << "Activated visualization output." << endl;
    if (ctx.opt._OUTPUT_pointcloud == 1)
        cout << "Activated Module: Point cloud output (simplified ALS simulation)" << endl; // v.3.1.6

    //!*********************
    //!** Evolution loop  **
    //!*********************

    cout << "Simulation starts with " << ctx.diag.nblivetrees << " trees." << endl;

    //** Information in file info **
    //******************************
    ctx.out.output_info << "\nTROLL simulator\n\n";
    ctx.out.output_info << "\n   2D discrete network: horizontal step = " << ctx.grid.LH
                        << " m, one tree per " << ctx.grid.LH * ctx.grid.LH << " m^2 \n\n";
    ctx.out.output_info << "\n   Tree : (t_dbh,t_height,t_CR,t_CD) \n\n";
    ctx.out.output_info << "\n            + one species label \n\n";
    ctx.out.output_info << " Number of sites      : " << ctx.grid.rows << "x" << ctx.grid.cols << "\n";
    ctx.out.output_info << " Number of iterations : " << ctx.time.nbiter << "\n";
    ctx.out.output_info << " Duration of timestep : " << ctx.time.timestep << " years\n";
    ctx.out.output_info << " Number of Species    : " << ctx.grid.nbspp << "\n\n";
    ctx.out.output_info.flush();

    // initial pattern, should be empty, unless an inventory has been provided
    if (ctx.opt._OUTPUT_extended & !ctx.opt._OUTPUT_inventory)
        OutputSnapshot(ctx, ctx.out.output_basic[1], 1, 0.01); // Initial Pattern, for trees > 0.01m DBH
    else if (ctx.opt._OUTPUT_inventory)
        OutputSnapshot(ctx, ctx.out.output_basic[1], 1, 0.001);
    else
        OutputSnapshot(ctx, ctx.out.output_basic[1], 1, 0.1); // Initial Pattern, for trees > 0.1m DBH

    double start_time, stop_time, duration = 0.0; // for simulation duration
    stop_time = clock();
    for (ctx.time.iter = 0; ctx.time.iter < ctx.time.nbiter; ctx.time.iter++)
    {
        start_time = stop_time;

        Evolution(ctx);
        stop_time = clock();
        duration += fmaxf(stop_time - start_time, 0.0);

        if (ctx.opt._OUTPUT_extended == 1 && ctx.crown.extent_visual > 0)
        {
            int timeofyear = GetTimeofyear(ctx);
            if (timeofyear == 0)
                OutputVisual(ctx);
        }

        /*if(ctx.opt._OUTPUT_pointcloud > 0 && ctx.time.iter == ctx.pc.iter_pointcloud_generation){
            ExportPointcloud(ctx.pc.mean_beam_pc, ctx.pc.sd_beam_pc, ctx.pc.klaser_pc, ctx.pc.transmittance_laser, ctx.out.output_pointcloud); // v.3.1.6
        }*/

#ifdef Output_ABC
        int timespan_abc = 10 * ctx.time.iterperyear; // every 10 years, modified in v.3.0
        int last_abc = ((ctx.time.nbiter + 5) / timespan_abc);
        last_abc *= timespan_abc;
        last_abc += -5;
        if ((ctx.time.iter + 5) % timespan_abc == 0 || ctx.time.iter == last_abc - 43)
        {
            // simulated lidar, based on data from Nouragues ALS, assuming NIR laser
            float transmittance_nir = 0.4; // transmittance of leaves in near infrared, used to calculate chances of laser hit
            float mean_beam = 12.0;        // avg number of beams that reach top of canopy
            float sd_beam = 5.0;           // standard deviation of those

            UpdateTransmittanceCHM_ABC(mean_beam, sd_beam, ctx.params.klight, transmittance_nir);
            OutputABC(ctx);
            UpdateDBHtrackingABC(ctx);
        }
#endif
    }

    cout << "Simulation ends with " << ctx.diag.nblivetrees << " trees." << endl;

    // final pattern
    if (ctx.opt._OUTPUT_extended & !ctx.opt._OUTPUT_inventory)
    {
        OutputSnapshot(ctx, ctx.out.output_basic[2], 1, 0.01); // Final Pattern, for trees > 0.01m DBH
    }
    else if (ctx.opt._OUTPUT_inventory)
    {
        OutputSnapshot(ctx, ctx.out.output_basic[2], 1, 0.001);
    }
    else
    {
        OutputSnapshot(ctx, ctx.out.output_basic[2], 1, 0.1); // Final Pattern, for trees > 0.1m DBH
    }
    if (ctx.opt._OUTPUT_extended)
    {
        OutputLAI(ctx, ctx.out.output_extended[7]);
        OutputCHM(ctx, ctx.out.output_extended[8]);
    }
    if (ctx.opt._OUTPUT_inventory)
    {
        for (int d = 0; d < ctx.grid.nbdcells; d = d + 1)
        {
            ctx.out.output_basic[3] << d;
            for (int l = 0; l < ctx.soil.nblayers_soil; l = l + 1)
            {
                ctx.out.output_basic[3] << "\t" << ctx.soil.SWC3D[l][d];
            }
            ctx.out.output_basic[3] << endl;
        }
    }

    //***********************
    //** End of simulation **
    //***********************
    float durf = duration / double(CLOCKS_PER_SEC); // output of the effective CPU time
    cout << "\n";
    ctx.out.output_info << "Average computation time: " << durf << " seconds.\n";
    ctx.out.output_info << "End of simulation.\n";
    ctx.out.output_info.flush();
    cout << "Average computation time : " << durf << " seconds.\n";
    cout << "End of simulation.\n";

    CloseOutputs(ctx); // new in v.3.1: Close and clear outputs, maybe not necessary as main function terminates shortly after, but maybe it ensures a cleaner communication with file system/within Rcpp
    FreeMem(ctx);      // Free dynamic memory  //! added in oct2013
    exit(0);
}
