#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

// ##############################################
//  Global function: calculation of the global averages every ctx.time.timestep
// ##############################################

void Average(Context &ctx)
{

#ifdef TRACK_INDIVIDUALS
    TrackingData_andOutput(ctx);
#endif

#ifdef Output_ABC
    UpdateMovingAveragesABC(ctx);
#endif

    int site, spp;
    float sum1 = 0.0, sum10 = 0.0, sum30 = 0.0, ba = 0.0, npp = 0.0, gpp = 0.0, ba10 = 0.0, agb = 0.0, rday = 0.0, rnight = 0.0, rstem = 0.0, litterfall = 0.0;

    float inbcells = 1.0 / float(ctx.grid.sites * mpi_size);
    float inbhectares = inbcells * ctx.grid.NH * ctx.grid.NH * 10000.0;

#ifdef WATER
    ctx.soil.abund_phi_root = 0.0;
    ctx.soil.abund10_phi_root = 0.0;
    ctx.soil.agb_phi_root = 0.0;
#endif

    if (!mpi_rank)
    {
        for (spp = 1; spp <= ctx.grid.nbspp; spp++)
            ctx.S[spp].s_sum10 = ctx.S[spp].s_sum30 = ctx.S[spp].s_ba = ctx.S[spp].s_ba10 = ctx.S[spp].s_agb = ctx.S[spp].s_gpp = ctx.S[spp].s_npp = ctx.S[spp].s_rday = ctx.S[spp].s_rnight = ctx.S[spp].s_rstem = ctx.S[spp].s_litterfall = 0;

        for (site = 0; site < ctx.grid.sites; site++)
            ctx.T[site].Average(ctx);

        for (spp = 1; spp <= ctx.grid.nbspp; spp++)
        {
            float s_sum1 = float(ctx.S[spp].s_nbind) * inbhectares;
            ctx.S[spp].s_sum10 *= inbhectares;
            ctx.S[spp].s_sum30 *= inbhectares;
            ctx.S[spp].s_ba *= inbhectares;
            ctx.S[spp].s_ba10 *= inbhectares;
            ctx.S[spp].s_agb *= inbhectares;
            ctx.S[spp].s_gpp *= inbhectares;
            ctx.S[spp].s_npp *= inbhectares;
            ctx.S[spp].s_rday *= inbhectares;
            ctx.S[spp].s_rnight *= inbhectares;
            ctx.S[spp].s_rstem *= inbhectares;
            ctx.S[spp].s_litterfall *= inbhectares;

            sum1 += float(ctx.S[spp].s_nbind) * inbhectares;
            sum10 += ctx.S[spp].s_sum10;
            sum30 += ctx.S[spp].s_sum30;
            ba += ctx.S[spp].s_ba;
            ba10 += ctx.S[spp].s_ba10;
            agb += ctx.S[spp].s_agb;
            gpp += ctx.S[spp].s_gpp;
            npp += ctx.S[spp].s_npp;
            rday += ctx.S[spp].s_rday;
            rnight += ctx.S[spp].s_rnight;
            rstem += ctx.S[spp].s_rstem;
            litterfall += ctx.S[spp].s_litterfall;

            // if(ctx.opt._OUTPUT_extended){
            ctx.out.output_extended[0] << ctx.time.iter << "\t" << ctx.S[spp].s_name << "\t" << s_sum1 << "\t" << ctx.S[spp].s_sum10 << "\t" << ctx.S[spp].s_sum30 << "\t" << ctx.S[spp].s_ba << "\t" << ctx.S[spp].s_ba10 << "\t" << ctx.S[spp].s_agb << "\t" << ctx.S[spp].s_gpp << "\t" << ctx.S[spp].s_npp << "\t" << ctx.S[spp].s_rday << "\t" << ctx.S[spp].s_rnight << "\t" << ctx.S[spp].s_rstem << "\t" << ctx.S[spp].s_litterfall << endl;
            //}
        }

        ctx.out.output_basic[0] << ctx.time.iter << "\t" << sum1 << "\t" << sum10 << "\t" << sum30 << "\t" << ba << "\t" << ba10 << "\t" << agb << "\t" << gpp << "\t" << npp << "\t" << rday << "\t" << rnight << "\t" << rstem << "\t" << litterfall << endl;

#ifdef MIP_Lichstein
        if ((!ctx.opt._FromInventory && ctx.time.iter >= (ctx.time.nbiter - 100 * ctx.time.iterperyear)) || ctx.opt._FromInventory)
        {
            ctx.out.output_MIP_eco << ctx.time.iter << "\t" << ctx.time.iter << "\t" << ctx.time.iter << "\t" << gpp * 100 << "\t" << npp * 100 << "\t";
        }
#endif

        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(2);

#ifdef WATER
        cout << ctx.time.iter << "\tTrees (1/ha): " << sum1 << " | " << sum10 << " | " << sum30 << " *** nbdead (%): " << 100.0 * ctx.diag.nbdead_n1 * inbhectares / sum1 << " | " << 100.0 * ctx.diag.nbdead_n10 * inbhectares / sum10 << " | " << 100.0 * ctx.diag.nbdead_n30 * inbhectares / sum30 << " *** AGB (t/ha): " << round(agb / 1000.0) << " GPP (MgC/ha/yr) " << gpp * ctx.time.iterperyear << " NPP " << npp * ctx.time.iterperyear << " litterfall (Mg/ha/yr) " << litterfall * ctx.time.iterperyear << " *** ctx.soil.Transpiration (mm): ";
#else

        cout << ctx.time.iter << "\tTrees (1/ha): " << sum1 << " | " << sum10 << " | " << sum30 << " *** nbdead (%): " << 100.0 * ctx.diag.nbdead_n1 * inbhectares / sum1 << " | " << 100.0 * ctx.diag.nbdead_n10 * inbhectares / sum10 << " | " << 100.0 * ctx.diag.nbdead_n30 * inbhectares / sum30 << " *** AGB (t/ha): " << round(agb / 1000.0) << " GPP (MgC/ha/yr) " << gpp * ctx.time.iterperyear << " NPP " << npp * ctx.time.iterperyear << " litterfall (Mg/ha/yr) " << litterfall * ctx.time.iterperyear << endl;
#endif

        if (ctx.opt._OUTPUT_extended)
        {
            float tototest = 0.0, tototest2 = 0.0, flux;
            for (int site = 0; site < ctx.grid.sites; site++)
            {
                flux = ctx.climate.WDailyMean * exp(-fmaxf(ctx.field.LAI3D[0][site + ctx.grid.SBORD], 0.0) * ctx.params.kpar);
                tototest += flux;
                tototest2 += flux * flux;
            }
            tototest /= float(ctx.grid.sites * ctx.grid.LH * ctx.grid.LH); // Average light flux (PPFD) on the ground
            tototest2 /= float(ctx.grid.sites * ctx.grid.LH * ctx.grid.LH);
            if (ctx.time.iter)
                ctx.out.output_extended[1] << ctx.time.iter << "\tMean PPFDground\t" << tototest << "\t" << sqrt(tototest2 - tototest * tototest) << "\n";

            if (ctx.opt._BASICTREEFALL)
                ctx.out.output_extended[2] << ctx.time.iter << "\t" << ctx.diag.nbdead_n1 * inbhectares << "\t" << ctx.diag.nbdead_n10 * inbhectares << "\t" << ctx.diag.nbTreefall1 * inbhectares << "\t" << ctx.diag.nbTreefall10 * inbhectares << endl;
            else
                ctx.out.output_extended[2] << ctx.time.iter << "\t" << ctx.diag.nbdead_n1 * inbhectares << "\t" << ctx.diag.nbdead_n10 * inbhectares << endl;
        }
    }

    if (ctx.opt._NDD)
        ctx.crown.BAtot = ba;

    for (int site = 0; site < ctx.grid.sites; site++)
    {
        if (ctx.T[site].t_age > 0)
        {
            if (ctx.T[site].t_NPP <= 0.0)
            {
                ctx.diag.nbtrees_carbstarv_n1++;
                if (ctx.T[site].t_dbh >= 0.1)
                    ctx.diag.nbtrees_carbstarv_n10++;
                if (ctx.T[site].t_dbh >= 0.3)
                    ctx.diag.nbtrees_carbstarv_n30++;
            }
        }
    }

#ifdef CHECK_CARBON
    float carbon_stored_leaves_previous, carbon_stored_trunk_previous = 0.0, carbon_stored_free_previous = 0.0, carbon_assimilated_total_previous = 0.0, carbon_net_total_previous = 0.0;

    if (ctx.time.iter == 0)
    {
        carbon_stored_leaves_previous = 0.0;
        carbon_stored_trunk_previous = 0.0;
        carbon_stored_free_previous = 0.0;
    }
    else
    {
        carbon_stored_leaves_previous = carbon_stored_leaves;
        carbon_stored_trunk_previous = carbon_stored_trunk;
        carbon_stored_free_previous = carbon_stored_free;
    }
    carbon_stored_leaves = 0.0;
    carbon_stored_trunk = 0.0;
    carbon_stored_free = 0.0;
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        if (ctx.T[s].t_age > 0)
        {
            float agb = 1000.0 * ctx.T[s].CalcAGB(); // convert to g
            float carbon_trunk = agb * 0.5;
            carbon_stored_trunk += carbon_trunk;

            float carbon_free = ctx.T[s].t_carbon_storage;
            carbon_stored_free += carbon_free;

            float carbon_leaves = ctx.T[s].t_LA * ctx.T[s].t_LMA * 0.5;
            carbon_stored_leaves += carbon_leaves;
        }
    }

    float factor_weight = 0.000001; // factor to convert carbon from g to tons

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(5);

    cout << ctx.time.iter << "\tTrunkC: " << carbon_stored_trunk * factor_weight << " LeavesC: " << carbon_stored_leaves * factor_weight << " FreeC: " << carbon_stored_free * factor_weight << " Total AssimC: " << carbon_assimilated_total * factor_weight << " Total NetC: " << carbon_net_total * factor_weight << endl;
    cout << ctx.time.iter << "\tTrunkC change: " << (carbon_stored_trunk - carbon_stored_trunk_previous) * factor_weight << " LeavesC: " << (carbon_stored_leaves - carbon_stored_leaves_previous) * factor_weight << " FreeC: " << (carbon_stored_free - carbon_stored_free_previous) * factor_weight << endl;

#endif

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(5);

#ifdef WATER

    float evapo = 0.0, runoff = 0.0, leak = 0.0, interception = 0.0, throughfall = 0.0, lai = 0.0;
    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        evapo += ctx.soil.Evaporation[d];         // in m3
        interception += ctx.soil.Interception[d]; // in mm, as rainfall
        throughfall += ctx.soil.Throughfall[d];   // in m3
        runoff += ctx.soil.Runoff[d];             // in m3
        leak += ctx.soil.Leakage[d];              // in m3
        lai += ctx.soil.LAI_DCELL[0][d];
    }

    float isites = 1.0 / float(ctx.grid.sites * ctx.grid.LH * ctx.grid.LH);
    float icells = 1.0 / float(ctx.grid.nbdcells);
    evapo *= isites;                // in ctx.params.m
    runoff *= isites;               // in ctx.params.m
    leak *= isites;                 // in ctx.params.m
    throughfall *= isites;          // in ctx.params.m
    interception *= isites * 0.001; // in ctx.params.m
    lai *= icells;
    ctx.soil.transpiration_1016 *= isites;

    ctx.out.output[11] << ctx.time.iter << "\t" << ctx.climate.precip << "\t" << interception << "\t" << throughfall << "\t" << runoff << "\t" << leak << "\t" << evapo << "\t";
    ctx.out.output[21] << ctx.time.iter << "\t" << lai << endl;

#ifdef MIP_Lichstein
    float transpitot = 0.0;
#endif

    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        float transpi = 0.0;
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            transpi += ctx.soil.Transpiration[l][d]; // in m3
        }
        transpi *= isites; // in ctx.params.m
        ctx.out.output[11] << transpi << "\t";
        cout << transpi * 1000 << " | "; // in mm

#ifdef MIP_Lichstein
        transpitot += transpi;
#endif
    }

    ctx.out.output[11] << ctx.soil.transpiration_1016 << "\t";

    cout << ctx.soil.transpiration_1016 * 1000 << " | " << endl;

#ifdef MIP_Lichstein
    if ((!ctx.opt._FromInventory && ctx.time.iter >= (ctx.time.nbiter - 100 * ctx.time.iterperyear)) || ctx.opt._FromInventory)
    {
        ctx.out.output_MIP_eco << (transpitot + evapo) * 1000 << "\t" << lai << "\t" << litterfall * 0.5 * 100 << "\t";
    }

    float SW1 = 0.0, SW2 = 0.0, SW3 = 0.0, SW4 = -9999; // Guyaflux
    // float SW1=0.0, SW2=0.0, SW3=0.0, SW4=0.0; // Tapajos
    float LT1 = 0.0, LT2 = 0.0, LT3 = 0.0, LT4 = 0.0;
#endif

    float layer_depth_previous = 0.0;
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        float soilWC = 0.0;
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            soilWC += ctx.soil.SWC3D[l][d]; // in m3
        }
        float layer_depth_current = ctx.soil.layer_depth[l];
        float layer_thickness = layer_depth_current - layer_depth_previous;
        soilWC *= isites / layer_thickness; // in m3/m3
        ctx.out.output[11] << soilWC << "\t";

#ifdef MIP_Lichstein
        if (l == 0 || l == 1)
        { // Guyaflux
            SW1 += soilWC * layer_thickness;
            LT1 += layer_thickness;
        }
        if (l == 2 || l == 3)
        {
            SW2 += soilWC * layer_thickness;
            LT2 += layer_thickness;
        }
        if (l == 4)
        {
            SW3 += soilWC * layer_thickness;
            LT3 += layer_thickness;
        }

        /*if (l==0 || l==1) {             //Tapajos
            SW1+=soilWC*layer_thickness;
            LT1+=layer_thickness;
        }
        if (l==2) {
            SW2+=soilWC*layer_thickness;
            LT2+= layer_thickness;
        }
        if (l==3) {
            SW3+=soilWC*layer_thickness;
            LT3+= layer_thickness;
        }
        if (l==4) {
            SW4+=soilWC*layer_thickness;
            LT4+= layer_thickness;
        }*/
#endif

        layer_depth_previous = layer_depth_current;
    }
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        float soilPhi = 0.0;
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            soilPhi += ctx.soil.soil_phi3D[l][d]; // in MPa
        }
        soilPhi *= icells; // in MPa
        ctx.out.output[11] << soilPhi << "\t";
    }

    ctx.out.output[11] << "\n";

#ifdef MIP_Lichstein
    SW1 /= LT1;
    SW2 /= LT2;
    SW3 /= LT3;
    // SW4/=LT4;     //Tapajos

    if ((!ctx.opt._FromInventory && ctx.time.iter >= (ctx.time.nbiter - 100 * ctx.time.iterperyear)) || ctx.opt._FromInventory)
    {
        ctx.out.output_MIP_eco << SW1 << "\t" << SW2 << "\t" << SW3 << "\t" << SW4 << endl;
    }
#endif

    if (ctx.time.iter == (ctx.time.nbiter - 90) || ctx.time.iter == (ctx.time.nbiter - 45) || ctx.time.iter == (ctx.time.nbiter - 1))
    {

        int o_swc, o_swp, o_wfluxes;
        if (ctx.time.iter == (ctx.time.nbiter - 90))
        {
            o_swc = 1;
            o_swp = 4;
            o_wfluxes = 18;
        }
        if (ctx.time.iter == (ctx.time.nbiter - 45))
        {
            o_swc = 2;
            o_swp = 5;
            o_wfluxes = 19;
        }
        if (ctx.time.iter == (ctx.time.nbiter - 1))
        {
            o_swc = 3;
            o_swp = 6;
            o_wfluxes = 20;
        }

        ctx.out.output[o_wfluxes] << "LAI" << "\t";
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            ctx.out.output[o_wfluxes] << ctx.soil.LAI_DCELL[0][d] << "\t";
        }
        ctx.out.output[o_wfluxes] << endl;

        ctx.out.output[o_wfluxes] << "Evaporation" << "\t";
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            ctx.out.output[o_wfluxes] << ctx.soil.Evaporation[d] * ctx.grid.i_sites_per_dcell << "\t"; // in ctx.params.m
        }
        ctx.out.output[o_wfluxes] << endl;

        layer_depth_previous = 0.0;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            float layer_depth_current = ctx.soil.layer_depth[l];
            float layer_thickness = layer_depth_current - layer_depth_previous;
            float norm = ctx.grid.i_sites_per_dcell / layer_thickness;
            ctx.out.output[o_swc] << l << "\t";
            ctx.out.output[o_swp] << l << "\t";
            ctx.out.output[o_wfluxes] << "Transpiration_" << l << "\t";
            for (int d = 0; d < ctx.grid.nbdcells; d++)
            {
                ctx.out.output[o_swc] << ctx.soil.SWC3D[l][d] * norm << "\t"; // in m3/m3
                ctx.out.output[o_swp] << ctx.soil.soil_phi3D[l][d] << "\t";
                ctx.out.output[o_wfluxes] << ctx.soil.Transpiration[l][d] * ctx.grid.i_sites_per_dcell << "\t";
            }
            layer_depth_previous = layer_depth_current;
            ctx.out.output[o_swc] << endl;
            ctx.out.output[o_swp] << endl;
            ctx.out.output[o_wfluxes] << endl;
        }
    }

    ctx.out.output[22] << ctx.time.iter << "\t";
    ctx.out.output[23] << ctx.time.iter << "\t";
    ctx.out.output[24] << ctx.time.iter << "\t";
    for (int l = 0; l < ctx.grid.HEIGHT + 1; l++)
    {
        ctx.soil.LAI_young[l] *= isites;
        ctx.soil.LAI_mature[l] *= isites;
        ctx.soil.LAI_old[l] *= isites;
        ctx.out.output[22] << ctx.soil.LAI_young[l] << "\t";
        ctx.out.output[23] << ctx.soil.LAI_mature[l] << "\t";
        ctx.out.output[24] << ctx.soil.LAI_old[l] << "\t";
        ctx.soil.LAI_young[l] = 0.0;
        ctx.soil.LAI_mature[l] = 0.0;
        ctx.soil.LAI_old[l] = 0.0;
    }
    ctx.out.output[22] << endl;
    ctx.out.output[23] << endl;
    ctx.out.output[24] << endl;

    ctx.soil.abund_phi_root *= inbhectares / sum1;
    ctx.soil.abund10_phi_root *= inbhectares / sum10;
    ctx.soil.agb_phi_root *= inbhectares / agb;

    ctx.out.output[31] << ctx.time.iter << "\t" << ctx.soil.abund_phi_root << "\t" << ctx.soil.abund10_phi_root << "\t" << ctx.soil.agb_phi_root << endl;

#endif

#ifdef MPI
    // This section corresponds to the parallel version of the reporting of the global diagnostic variables. Since much work has been done on routine Average over the past years, this would need a full rewrite, !!!!Action 20/01/2016: rework the parallel version of function Average!!!!

    //     MPI_Reduce(&(ctx.S[spp].s_nbind),&sind,1,
    //     MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    //     MPI_Reduce(ctx.S[spp].s_output_field,ctx.S[spp].s_output_field,5,
    //     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    //     MPI_Reduce(Mortality,Mortality,4,
    //     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    //     MPI_Reduce(&ctx.S[spp].s_output_field[6],&ctx.S[spp].s_output_field[6],5,
    //     MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
#endif
    cout.flush();
}

// ##############################################
//  Global function: output of the field variables every ctx.time.timestep
// ##############################################
void OutputField(Context &ctx)
{
    int site, h;
    if ((ctx.time.nbout) && ((ctx.time.iter % ctx.time.freqout) == ctx.time.freqout - 1))
    {
        // output fields, ctx.time.nbout times during simulation (every ctx.time.freqout iterations)
        int d;
        for (d = 0; d < ctx.grid.dbhmaxincm; d++)
            ctx.diag.nbdbh[d] = 0;
        for (site = 0; site < ctx.grid.sites; site++)
            ctx.T[site].histdbh(ctx);

        for (h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        {
            ctx.diag.layer[h] = 0;
            for (site = 0; site < ctx.grid.sites; site++)
                ctx.diag.layer[h] += ctx.field.LAI3D[h][site + ctx.grid.SBORD];
        }

#ifdef MPI
        MPI_Status status;
        MPI_Reduce(ctx.diag.nbdbh, ctx.diag.nbdbh, ctx.grid.dbhmaxincm, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(ctx.diag.layer, ctx.diag.layer, ctx.grid.HEIGHT, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        if (!mpi_rank)
        {
            // output of the dbh histograms (ctx.out.output[31])
            for (d = 1; d < ctx.grid.dbhmaxincm; d++)
                ctx.out.output[31] << d << "\t" << ctx.diag.nbdbh[d] << "\n";
            ctx.out.output[31] << "\n";
            // output of the mean LAI per height class (ctx.out.output[32])
            float norm = 1.0 / float(ctx.grid.sites * ctx.grid.LH * ctx.grid.LH * mpi_size);
            for (h = 0; h < (ctx.grid.HEIGHT + 1); h++)
                ctx.out.output[32] << ctx.time.iter << "\t" << h * ctx.grid.LV << "\t" << ctx.diag.layer[h] * norm << "\n";
            ctx.out.output[32] << "\n";
        }
    }
}

// ##############################################
//  Global function: output snapshots of the scene at one point in time
// ##############################################
//! - This can be used to take snapshots of the forest in more detail and track its development over time.
void OutputSnapshot(Context &ctx, fstream &output, bool header, float dbh_limit)
{
    cout << "Writing snapshot of forest to file." << endl;
    if (header == 1)
    {

#ifdef WATER
        output << "iter\tcol\trow\tfrom_Data\tsp_lab\tsite\tCrownDisplacement\tPmass\tNmass\tLMA\tleafarea\ttlp\twsg\tRdark\tVcmax\tJmax\tleaflifespan\tlambda_young\tlambda_mature\tlambda_old\tdbhmature\tdbhmax\thmax\tah\tCt\tLAImax\tfraction_filled\tmult_height\tmult_CR\tmult_CD\tmult_P\tmult_N\tmult_LMA\tmult_dbhmax\tmult_leafarea\tmult_tlp\tdev_wsg\tage\tdbh\tsapwood_area\theight\tCD\tCR\tGPP\tNPP\tRday\tRnight\tRstem\tLAmax\tLA\tyoungLA\tmatureLA\toldLA\tLAI\tlitter\tcarbon_storage\tcarbon_biometry\tmultiplier_seed\thurt\tNPPneg";

        output << "\troot_depth\tphi_root\tWSF\tWSF_A\ttranspiration\tg1_0\tg1";
#ifdef PHENO_DROUGHT
        output << "\tNdays_dry\tNdays_wet\tpheno_factor";
#endif
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            output << "\troot_biomass" << l;
        }
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            output << "\tsoil_layer_weight" << l;
        }

#else // WATER
        output << "iter\tcol\trow\tfrom_Data\tsp_lab\tsite\tCrownDisplacement\tPmass\tNmass\tLMA\twsg\tRdark\tVcmax\tJmax\tleaflifespan\tlambda_young\tlambda_mature\tlambda_old\tdbhmature\tdbhmax\thmax\tah\tCt\tLAImax\tfraction_filled\tmult_height\tmult_CR\tmult_CD\tmult_P\tmult_N\tmult_LMA\tmult_dbhmax\tdev_wsg\tage\tdbh\tsapwood_area\theight\tCD\tCR\tGPP\tNPP\tRday\tRnight\tRstem\tLAmax\tLA\tyoungLA\tmatureLA\toldLA\tLAI\tlitter\tcarbon_storage\tcarbon_biometry\tmultiplier_seed\thurt\tNPPneg";

#endif

#ifdef Output_ABC
        output << "\tdbh_previous";
#endif
        output << "\tAGB\ts_name" << endl;
    }

    // reset the canopy to make sure that it is well-constructed
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        for (int sbsite = 0; sbsite < ctx.grid.sites + 2 * ctx.grid.SBORD; sbsite++)
            ctx.field.LAI3D[h][sbsite] = 0.0;

    for (int site = 0; site < ctx.grid.sites; site++)
    { // Each tree contribues to ctx.field.LAI3D
        CalcLAI(ctx, ctx.T[site]);
    }

    for (int h = ctx.grid.HEIGHT; h > 0; h--)
    { // LAI is computed by summing LAI from the canopy top to the ground
        for (int site = 0; site < ctx.grid.sites; site++)
        {
            int sbsite = site + ctx.grid.SBORD;
            ctx.field.LAI3D[h - 1][sbsite] += ctx.field.LAI3D[h][sbsite];
        }
    }

    // output.setf(ios::fixed,ios::floatfield);
    if (ctx.opt._OUTPUT_inventory)
    {
        output.precision(10); // needed for correct restart
    }
    else
    {
        output.precision(5);
    }
    for (int row = 0; row < ctx.grid.rows; row++)
    {
        for (int col = 0; col < ctx.grid.cols; col++)
        {
            int site = col + ctx.grid.cols * row;
            if (ctx.T[site].t_age > 0 && ctx.T[site].t_dbh >= dbh_limit)
            {
                // recalculate photosynthesis and respiration
                // ctx.T[site].CalcRespGPP();
                // ctx.T[site].CalcNPP();
                // output all tree variables, this is potentially a very large file
                // we currently do not output the t_NDDfield vector, as it is too large

#ifdef WATER
                output << ctx.time.iter << "\t" << col << "\t" << row << "\t" << ctx.T[site].t_from_Data << "\t" << ctx.T[site].t_sp_lab << "\t" << site << "\t" << ctx.T[site].t_CrownDisplacement << "\t" << ctx.T[site].t_Pmass << "\t" << ctx.T[site].t_Nmass << "\t" << ctx.T[site].t_LMA << "\t" << ctx.T[site].t_leafarea << "\t" << ctx.T[site].t_tlp << "\t" << ctx.T[site].t_wsg << "\t" << ctx.T[site].t_Rdark << "\t" << ctx.T[site].t_Vcmax << "\t" << ctx.T[site].t_Jmax << "\t" << ctx.T[site].t_leaflifespan << "\t" << ctx.T[site].t_lambda_young << "\t" << ctx.T[site].t_lambda_mature << "\t" << ctx.T[site].t_lambda_old << "\t" << ctx.T[site].t_dbhmature << "\t" << ctx.T[site].t_dbhmax << "\t" << ctx.T[site].t_hmax << "\t" << ctx.T[site].t_ah << "\t" << ctx.T[site].t_Ct << "\t" << ctx.T[site].t_LAImax << "\t" << ctx.T[site].t_fraction_filled << "\t" << ctx.T[site].t_mult_height << "\t" << ctx.T[site].t_mult_CR << "\t" << ctx.T[site].t_mult_CD << "\t" << ctx.T[site].t_mult_P << "\t" << ctx.T[site].t_mult_N << "\t" << ctx.T[site].t_mult_LMA << "\t" << ctx.T[site].t_mult_dbhmax << "\t" << ctx.T[site].t_mult_leafarea << "\t" << ctx.T[site].t_mult_tlp << "\t" << ctx.T[site].t_dev_wsg << "\t" << ctx.T[site].t_age << "\t" << ctx.T[site].t_dbh << "\t" << ctx.T[site].t_sapwood_area << "\t" << ctx.T[site].t_height << "\t" << ctx.T[site].t_CD << "\t" << ctx.T[site].t_CR << "\t" << ctx.T[site].t_GPP << "\t" << ctx.T[site].t_NPP << "\t" << ctx.T[site].t_Rday << "\t" << ctx.T[site].t_Rnight << "\t" << ctx.T[site].t_Rstem << "\t" << ctx.T[site].t_LAmax << "\t" << ctx.T[site].t_LA << "\t" << ctx.T[site].t_youngLA << "\t" << ctx.T[site].t_matureLA << "\t" << ctx.T[site].t_oldLA << "\t" << ctx.T[site].t_LAI << "\t" << ctx.T[site].t_litter << "\t" << ctx.T[site].t_carbon_storage << "\t" << ctx.T[site].t_carbon_biometry << "\t" << ctx.T[site].t_multiplier_seed << "\t" << ctx.T[site].t_hurt << "\t" << ctx.T[site].t_NPPneg;

                output << "\t" << ctx.T[site].t_root_depth << "\t" << ctx.T[site].t_phi_root << "\t" << ctx.T[site].t_WSF << "\t" << ctx.T[site].t_WSF_A << "\t" << ctx.T[site].t_transpiration << "\t" << ctx.T[site].t_g1 << "\t" << ctx.T[site].t_g1_0;
#ifdef PHENO_DROUGHT
                output << "\t" << ctx.T[site].t_Ndays_dry << "\t" << ctx.T[site].t_Ndays_wet << "\t" << ctx.T[site].t_pheno_factor;

#endif
                for (int l = 0; l < ctx.soil.nblayers_soil; l++)
                {
                    output << "\t" << ctx.T[site].t_root_biomass[l];
                }
                for (int l = 0; l < ctx.soil.nblayers_soil; l++)
                {
                    output << "\t" << ctx.T[site].t_soil_layer_weight[l];
                }

#else
                output << ctx.time.iter << "\t" << col << "\t" << row << "\t" << ctx.T[site].t_from_Data << "\t" << ctx.T[site].t_sp_lab << "\t" << site << "\t" << ctx.T[site].t_CrownDisplacement << "\t" << ctx.T[site].t_Pmass << "\t" << ctx.T[site].t_Nmass << "\t" << ctx.T[site].t_LMA << "\t" << ctx.T[site].t_wsg << "\t" << ctx.T[site].t_Rdark << "\t" << ctx.T[site].t_Vcmax << "\t" << ctx.T[site].t_Jmax << "\t" << ctx.T[site].t_leaflifespan << "\t" << ctx.T[site].t_lambda_young << "\t" << ctx.T[site].t_lambda_mature << "\t" << ctx.T[site].t_lambda_old << "\t" << ctx.T[site].t_dbhmature << "\t" << ctx.T[site].t_dbhmax << "\t" << ctx.T[site].t_hmax << "\t" << ctx.T[site].t_ah << "\t" << ctx.T[site].t_Ct << "\t" << ctx.T[site].t_LAImax << "\t" << ctx.T[site].t_fraction_filled << "\t" << ctx.T[site].t_mult_height << "\t" << ctx.T[site].t_mult_CR << "\t" << ctx.T[site].t_mult_CD << "\t" << ctx.T[site].t_mult_P << "\t" << ctx.T[site].t_mult_N << "\t" << ctx.T[site].t_mult_LMA << "\t" << ctx.T[site].t_mult_dbhmax << "\t" << ctx.T[site].t_dev_wsg << "\t" << ctx.T[site].t_age << "\t" << ctx.T[site].t_dbh << "\t" << ctx.T[site].t_sapwood_area << "\t" << ctx.T[site].t_height << "\t" << ctx.T[site].t_CD << "\t" << ctx.T[site].t_CR << "\t" << ctx.T[site].t_GPP << "\t" << ctx.T[site].t_NPP << "\t" << ctx.T[site].t_Rday << "\t" << ctx.T[site].t_Rnight << "\t" << ctx.T[site].t_Rstem << "\t" << ctx.T[site].t_LAmax << "\t" << ctx.T[site].t_LA << "\t" << ctx.T[site].t_youngLA << "\t" << ctx.T[site].t_matureLA << "\t" << ctx.T[site].t_oldLA << "\t" << ctx.T[site].t_LAI << "\t" << ctx.T[site].t_litter << "\t" << ctx.T[site].t_carbon_storage << "\t" << ctx.T[site].t_carbon_biometry << "\t" << ctx.T[site].t_multiplier_seed << "\t" << ctx.T[site].t_hurt << "\t" << ctx.T[site].t_NPPneg;
#endif

#ifdef Output_ABC
                output << "\t" << ctx.T[site].t_dbh_previous;
#endif

                // we add a few tree-based variables that are derived or environment-related, but not directly kept track of
                float AGB = ctx.T[site].CalcAGB(ctx);

                output << "\t" << AGB << "\t" << ctx.S[ctx.T[site].t_sp_lab].s_name << endl;
            }
        }
    }
}

#ifdef CHM_SPIKEFREE
// ##########################
// ## Make a spikefree CHM ##
// ##########################
void MakeCHMspikefree(Context &ctx, vector<int> &chm_spikefree)
{
    chm_spikefree.clear();
    chm_spikefree.reserve(ctx.grid.sites);
    for (int s = 0; s < ctx.grid.sites; s++)
        chm_spikefree.push_back(0);

#ifdef Output_ABC // IM2023 added, otherwise row_start, row_end, col_start, col_end are not necessarily decalred, and this blocks code building.

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
            {
#ifdef CROWN_UMBRELLA
                float height = ctx.T[s].t_height;
                float CR = ctx.T[s].t_CR;
                float CD = ctx.T[s].t_CD;
                float fraction_filled_target = 1.0; // to remove spikes, we assume a fully filled crown for all trees
                int shell_fromtop = 0;              // toplayer
                float noinput = 0.0;

                LoopLayerUpdateCrownStatistic_template(ctx, r, c, height, CR, CD, fraction_filled_target, shell_fromtop, [&ctx](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); }, noinput, chm_spikefree, KeepFloatAsIs, UpdateCHMvector);
#else
                int crown_top = int(ctx.T[s].t_height);
                int crown_intarea = GetCrownIntarea(ctx.T[s].t_CR);

                int row_crowncenter = s / ctx.grid.cols;
                int col_crowncenter = s % ctx.grid.cols;

                for (int i = 0; i < crown_intarea; i++)
                {
                    int site_relative = ctx.lookup.LookUp_Crown_site[i];
                    int row = row_crowncenter + site_relative / 51 - 25;
                    int col = col_crowncenter + site_relative % 51 - 25;
                    if (row >= 0 && row < ctx.grid.rows && col >= 0 && col < ctx.grid.cols)
                    {
                        int site = col + row * ctx.grid.cols;
                        if (chm_spikefree[site] < crown_top)
                            chm_spikefree[site] = crown_top;
                    }
                }
#endif
            }
        }
    }

#endif
}
#endif

// ################################
// ### Output for visualization ###
// ################################
void OutputVisual(Context &ctx)
{
    // first simple chm output
#ifdef CHM_SPIKEFREE
    vector<int> chm_spikefree;
    MakeCHMspikefree(ctx, chm_spikefree);

    for (int col = ctx.crown.mincol_visual; col < ctx.crown.maxcol_visual; col++)
    {
        for (int row = ctx.crown.minrow_visual; row < ctx.crown.maxrow_visual; row++)
        {
            int site = col + row * ctx.grid.cols;
            int height_canopy = 0;
            for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            {
                if (ctx.field.LAI3D[h][site + ctx.grid.SBORD] > 0.0)
                    height_canopy = max(h, height_canopy);
            }
            ctx.out.output_visual[0] << ctx.time.iter << "\t" << row << "\t" << col << "\t" << height_canopy + 1 << "\t" << chm_spikefree[site] << "\t" << ctx.field.LAI3D[0][site + ctx.grid.SBORD] << endl;
        }
    }
#else
    for (int col = ctx.crown.mincol_visual; col < ctx.crown.maxcol_visual; col++)
    {
        for (int row = ctx.crown.minrow_visual; row < ctx.crown.maxrow_visual; row++)
        {
            int site = col + row * ctx.grid.cols;
            int height_canopy = 0;
            for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            {
                if (ctx.field.LAI3D[h][site + ctx.grid.SBORD] > 0.0)
                    height_canopy = max(h, height_canopy);
            }
            ctx.out.output_visual[0] << ctx.time.iter << "\t" << row << "\t" << col << "\t" << height_canopy + 1 << "\t" << ctx.field.LAI3D[0][site + ctx.grid.SBORD] << endl;
        }
    }
#endif

    // now the sliced output
    for (int row = ctx.crown.minrow_visual_slice; row < ctx.crown.maxrow_visual_slice; row++)
    {
        for (int col = 0; col < ctx.grid.cols; col++)
        {
            int s = col + row * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
            {
                int row_slice = row;

                float height = ctx.T[s].t_height;
                float CD = ctx.T[s].t_CD;
                int crown_top = int(height);
                float CR = ctx.T[s].t_CR;
                int crown_base = int(height - CD);

                vector<float> output_statistics;
                output_statistics.reserve(4);
                output_statistics.push_back(ctx.T[s].t_sp_lab);

                float ratio_height_Ct;
                if (ctx.T[s].t_Ct > 0.0)
                    ratio_height_Ct = ctx.T[s].t_height / ctx.T[s].t_Ct;
                else
                    ratio_height_Ct = 0.0;
                output_statistics.push_back(ratio_height_Ct);

                float ratio_NPP_GPP;
                if (ctx.T[s].t_GPP > 0.0)
                    ratio_NPP_GPP = ctx.T[s].t_NPP / ctx.T[s].t_GPP;
                else
                    ratio_NPP_GPP = 0.0;
                output_statistics.push_back(ratio_NPP_GPP);

#ifdef CROWN_UMBRELLA
                float fraction_filled_target = 1.0;                  // we assume a fully filled crown for all trees
                int max_shells = min(crown_top - crown_base + 1, 4); // since the new crown shapes

                for (int h = 0; h <= crown_top - max_shells; h++)
                    OutputCrownSliced(ctx, h, s, row_slice, output_statistics);

                for (int shell_fromtop = 0; shell_fromtop < max_shells; shell_fromtop++)
                {
                    LoopLayerUpdateCrownStatistic_template(ctx, row, col, height, CR, CD, fraction_filled_target, shell_fromtop, [&ctx](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); }, row_slice, output_statistics, KeepIntAsIs, [&ctx](int h, int s, int rs, vector<float> &os){ OutputCrownSliced(ctx, h, s, rs, os); });
                }
#else
                int crown_intarea = GetCrownIntarea(CR);

                for (int h = 0; h <= crown_base; h++)
                    OutputCrownSliced(ctx, h, s, row_slice, output_statistics);

                for (int h = crown_base; h <= crown_top; h++)
                {
                    for (int i = 0; i < crown_intarea; i++)
                    {
                        int site_relative = ctx.lookup.LookUp_Crown_site[i];
                        int row_crown = row + site_relative / 51 - 25;
                        int col_crown = col + site_relative % 51 - 25;
                        int site_crown = col_crown + row_crown * ctx.grid.cols;
                        OutputCrownSliced(ctx, h, s, row_slice, output_statistics);
                    }
                }
#endif
            }
        }
    }
}

// ##################
// ### Output CHM ###
// ##################
void OutputCHM(Context &ctx, fstream &output_CHM)
{
#ifdef CHM_SPIKEFREE
    vector<int> chm_spikefree;
    MakeCHMspikefree(ctx, chm_spikefree);

    output_CHM << "site" << "\t" << "row" << "\t" << "col" << "\t" << "height" << "\t" << "height_spikefree" << "\t" << "LAI" << endl;
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        int height_canopy = 0;
        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        {
            if (ctx.field.LAI3D[h][s + ctx.grid.SBORD] > 0.0)
                height_canopy = max(h, height_canopy);
        }
        output_CHM << s << "\t" << int(s / ctx.grid.cols) << "\t" << int(s % ctx.grid.cols) << "\t" << height_canopy + 1 << "\t" << chm_spikefree[s] << "\t" << ctx.field.LAI3D[0][s + ctx.grid.SBORD] << endl;
    }

#else
    output_CHM << "site" << "\t" << "row" << "\t" << "col" << "\t" << "height" << "\t" << "LAI" << endl;
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        int height_canopy = 0;
        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            if (ctx.field.LAI3D[h][s + ctx.grid.SBORD] > 0.0)
                height_canopy = max(h, height_canopy);
        output_CHM << s << "\t" << int(s / ctx.grid.cols) << "\t" << int(s % ctx.grid.cols) << "\t" << height_canopy + 1 << "\t" << ctx.field.LAI3D[0][s + ctx.grid.SBORD] << endl;
    }
#endif
}

// ##############################################
//  Global function: writes the whole 3D LAI voxel field to file
// ##############################################
void OutputLAI(Context &ctx, fstream &output_transmLAI3D)
{
    output_transmLAI3D << "s\trow\tcol\th\tLAI3D" << endl;
    for (int s = 0; s < ctx.grid.sites; s++)
        for (int h = 0; h < ctx.grid.HEIGHT; h++)
            output_transmLAI3D << s << "\t" << int(s / ctx.grid.cols) << "\t" << int(s % ctx.grid.cols) << "\t" << h << "\t" << ctx.field.LAI3D[h][s + ctx.grid.SBORD] << endl;
}

// ##############################################
//  Global function: transforming TROLL ctx.field.LAI3D field into point cloud and export as .las file
// ##############################################
//! - Basic ALS simulation, same as for calculating the TROLL transmittance field, added in v.3.1.6; !!!: TODO: merge with UpdateTransmittanceCHM_ABC
//! - First draws from a distribution to calculate the sampling density, i.e. the number of beams per voxel column
//! - Then loops over beams per voxel column from top to bottom and calculates the number of hits given the leaf area density (LAD) of the respective voxel
//! - If beams hit ground, all produce guaranteed returns
//! - The parameters used is the k assumed for the laser, which is only based on leaf geometry and should thus be taken equal to ctx.params.klight
//! - The difference to ctx.params.kpar is accounted for by parameterizing a probability of a further return after a hit, i.e. a transmittance of the leaves that are hit that is larger than for visible light and corresponds roughly to empirical fractions of second returns (~0.3-0.4). If taken together, a ctx.params.klight = 0.5 and a transmittance probability after a hit of 0.4 would yield an effective knir = 0.2, which is reasonable for most laser scanners.
//! - This approach makes lots of simplifying assumptions, among which are: no flightline, no angles, no beam diameter/divergence, Lambert-Beer extinction, including assumptions about the conversion between transmittance probability and energy extinction. Importantly, a TROLL forest contains no woody material and does not represent topography at the moment, all of which may influence sampling densities as well.
//! - One further question concerns the probability of obtaining a ground return: there is a long discussion in the literature about backscatter ratios between vegetation and ground. However, this seems to be neither a constant ratio nor seems there to be agreement on how it varies, likely due to dependence on footprint size and composition of reflecting surfaces. Cf., for example, Ni-Meister et al. 2001, IEEE Transactions: they use a ratio of mostly 1.0 (the implicit assumption in this simulation as well), but state that it should vary from site to site. Knapp et al. 2021, Remote Sensing, on the other hand: "The reflectance of the forest ground voxels was down-weighted by dividing by 2.5 in order to account for the lower reflectivity of the ground vs. vegetation." This specifically concerns GEDI waveform simulations, so large-footprint lasers, but the source publication does not seem to give a direct reference for the value. Furthermore, it would imply that the ratio of ground to vegetation reflectance is 0.4 for large-footprint full-waveform lidar, which seems to be in stark contrast to Chen et al. 2014, Remote Sensing of Environment, who looked at footprints of up to 8m and found a constant ratio ground/vegetation of 1.7 (they measured vegetation/ground ratio as ~0.57, so 1.7 = 1.0/0.57). However, their Figure 11 seems to indicate an inverse relationship for small footprints (~0.4m) with a ground to vegetation backscatter ratio of ca. 0.8. In the future, it would be worth running TROLL voxel output through an explicit raytracing simulator, with a much more detailed parameterization of surfaces (e.g. DART or Helios++, https://github.com/3dgeo-heidelberg/helios) and compare our simplified simulations to its outcomes/adjust accordingly.

void GenerateVoxelreturnsALS(Context &ctx, vector<int> &beams, vector<float> &beams_returns, float mean_beam, float sd_beam, float klaser, float transmittance_laser)
{

    int beams_expected = int(ctx.grid.sites * mean_beam * 1.01); // reserve a bit more
    int returns_maximum = beams_expected * 5;           // if every laser had 5 returns
    beams.reserve(beams_expected);
    beams_returns.reserve(returns_maximum);

    // loop through the whole array by site first, then by height, then by return number
    for (int r = 0; r < ctx.grid.rows; r++)
    {
        for (int c = 0; c < ctx.grid.cols; c++)
        {
            int site = c + r * ctx.grid.cols;
            int nbbeams = int(mean_beam + gsl_ran_gaussian(ctx.rng.gslrand, sd_beam)); // always rounding up
            nbbeams = max(nbbeams, 1);

            for (int beam = 0; beam < nbbeams; beam++)
            {
                // loop over the field from maximum height to 0 and iteratively update voxels from top to bottom, following the beam. An alternative version, also allowing for ground returns, can be activated to extending the loop to h >= -1. In this case, when a beam is not extinguished before it reaches the ground (h >= 0), then it is counted as a ground return

                beams.push_back(site); // store the location of the current beam

                int beam_continues = 1;
                int beam_return = 0;
                int h = ctx.grid.HEIGHT - 1; // starting point to sample from

                while (h >= -1 && beam_continues && beam_return < 5)
                {
                    // there will only be a return, if the beam still exists

                    if (h >= 0)
                    {
                        // returns due to vegetation
                        float LAI_above = ctx.field.LAI3D[h + 1][site + ctx.grid.SBORD];
                        float LAI_current = ctx.field.LAI3D[h][site + ctx.grid.SBORD];

                        float LAD = LAI_current - LAI_above;

                        float prob_hit;
                        if (LAD > 0.0)
                            prob_hit = 1.0 - exp(-klaser * LAD);
                        else
                            prob_hit = 0.0;

                        int is_hit = gsl_ran_bernoulli(ctx.rng.gslrand, prob_hit);

                        if (is_hit)
                        {
                            beam_continues = gsl_ran_bernoulli(ctx.rng.gslrand, transmittance_laser);
                            beam_return++;
                            float z_hit = float(h) + gsl_rng_uniform(ctx.rng.gslrand);
                            beams_returns.push_back(z_hit);
                        }
                    }
                    else
                    {
                        // ground returns
                        beams_returns.push_back(0.0);
                        beam_return++;
                        beam_continues = 0;
                    }

                    h--;
                }

                beams.push_back(beam_return); // add the return number
            }
        }
    }
}
// !!!: At the moment, we simply assume a little-endian system for output (most personal computers, but not servers)
// BELOW some hints on how we would check for endianness and implement conversion of big to little endian
// the other way round (little to big) would be easy, as this is implemented in c++ via htol function

// check endianness,cf. https://stackoverflow.com/questions/4181951/how-to-check-whether-a-system-is-big-endian-or-little-endian
// int n = 1;
//// little endian if true
// if(*(char *)&n == 1) {...}

// or with <bit> library, cf. https://en.cppreference.com/w/cpp/types/endian, but: C++20
// #include <bit>
// #include <iostream>
//
// int main() {
//
//    if constexpr (std::endian::native == std::endian::big)
//        std::cout << "big-endian\n";
//    else if constexpr (std::endian::native == std::endian::little)
//        std::cout << "little-endian\n";
//    else std::cout << "mixed-endian\n";
//}

// swapping template (then we just need to check the endianness of the system)
// https://mklimenko.github.io/english/2018/08/22/robust-endian-swap/
// template <typename T>
// void SwapEndian(T &val) {
//    union U {
//        T val;
//        std::array<std::uint8_t, sizeof(T)> raw;
//    } src, dst;
//
//    src.val = val;
//    std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
//    val = dst.val;
//}

/* void ExportPointcloudHeader(vector<int> &beams, fstream& ctx.out.output_pointcloud){
    // las files are defined as little endian
    // for the moment, we assume a little endian system and that chars actually have 8 bits (1 byte)
    // all names are just LAS definition names with underscores
    // future versions should upgrade format to 1.4, and include ways to parameterize coordinate reference system, etc.

    char file_signature[5] = "LASF";
    ctx.out.output_pointcloud.write(file_signature, sizeof(file_signature) - 1); // remove terminating NULL in char

    uint16_t file_source_id = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&file_source_id), sizeof(file_source_id));

    uint16_t global_encoding = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&global_encoding), sizeof(global_encoding));

    uint32_t project_id_guid_data_1 = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&project_id_guid_data_1), sizeof(project_id_guid_data_1));

    uint16_t project_id_guid_data_2 = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&project_id_guid_data_2), sizeof(project_id_guid_data_2));

    uint16_t project_id_guid_data_3 = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&project_id_guid_data_3), sizeof(project_id_guid_data_3));

    unsigned char project_id_guid_data_4[9] = "";
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&project_id_guid_data_4), sizeof(project_id_guid_data_4) - 1);

    unsigned char version_major = 1;    // this is slightly weird: why were these not simply defined as uint8_t in LAS specification? Or is "unsigned char" short for uint8_t? Because there is also an ascii char reserved for numbers (e.g. "1" corresponds to char = 49)
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&version_major), sizeof(version_major));

    unsigned char version_minor = 2;    // this is slightly weird: why were these not simply defined as uint8_t in LAS specification? Or is "unsigned char" short for uint8_t? Because there is also an ascii char reserved for numbers (e.g. "1" corresponds to char = 49)
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&version_minor), sizeof(version_minor));

    char system_identifier[33];
    sprintf(system_identifier,"ALS simulator");
    ctx.out.output_pointcloud.write(system_identifier, sizeof(system_identifier) - 1); // remove terminating NULL in char

    char generating_software[33] = "TROLL v.3.1.6+ forest simulator";
    ctx.out.output_pointcloud.write(generating_software, sizeof(generating_software) - 1); // remove terminating NULL in char

    // get current day and year, cf. https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
    time_t now = time(0);
    tm *ltm = localtime(&now);

    uint16_t file_creation_day_of_year = ltm->tm_yday;
    uint16_t file_creation_year = 1900 + ltm->tm_year;

    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&file_creation_day_of_year), sizeof(file_creation_day_of_year));
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&file_creation_year), sizeof(file_creation_year));

    uint16_t header_size = 227;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&header_size), sizeof(header_size));

    uint32_t offset_to_point_data = 227;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&offset_to_point_data), sizeof(offset_to_point_data));

    uint32_t number_of_variable_length_records = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&number_of_variable_length_records), sizeof(number_of_variable_length_records));

    unsigned char point_data_format_ID = 0; // no GPS time needed
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&point_data_format_ID), sizeof(point_data_format_ID));

    // each record in format 0 has 20 bytes (12 for coordinates, 2 for intensity, 6 for other information)
    uint16_t point_data_record_length = 20;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&point_data_record_length), sizeof(point_data_record_length));

    // calculate the number of returns per return number
    int nb_perreturn[5] = {0};
    int nb_beams = int(beams.size()/2);
    int nb_returns = 0;

    for(int beam = 0; beam < nb_beams; beam++){
        int nb_returns_beam = beams[1 + beam * 2];
        for(int rtrn = 0; rtrn < nb_returns_beam; rtrn++){
            nb_perreturn[rtrn]++;
            nb_returns++;
        }
    }

    cout << "Sampled " << nb_beams << " pulses, creating " << nb_returns << " returns." << endl;

    uint32_t number_of_point_records = nb_returns;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&number_of_point_records), sizeof(number_of_point_records));

    for(int i = 0; i < 5; i++){
        uint32_t number_of_points_by_return = nb_perreturn[i];
        ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&number_of_points_by_return), sizeof(number_of_points_by_return));
    }

    // there is no fixed-width type for floating-point numbers, so we assume that the 8 byte required by the las specification are fulfilled
    for(int i = 0; i < 3; i++){
        double xyz_scale_factor = 0.01;
        ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&xyz_scale_factor), 8); // hardcoded 8 bytes
    }

    for(int i = 0; i < 3; i++){
        double xyz_offset = 0.0;
        ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&xyz_offset), 8); // hardcoded 8 bytes
    }

    double max_x = ctx.grid.cols;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&max_x), 8); // hardcoded 8 bytes

    double min_x = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&min_x), 8); // hardcoded 8 bytes

    double max_y = ctx.grid.rows;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&max_y), 8); // hardcoded 8 bytes

    double min_y = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&min_y), 8); // hardcoded 8 bytes

    double max_z = ctx.grid.HEIGHT;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&max_z), 8); // hardcoded 8 bytes

    double min_z = 0;
    ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&min_z), 8); // hardcoded 8 bytes
}

 void ExportPointcloud(float mean_beam, float sd_beam, float klaser, float ctx.pc.transmittance_laser, fstream& ctx.out.output_pointcloud){
    cout << "Point cloud generation." << endl;

    vector<int> beams;
    vector<float> beams_returns;

    // three options
    GenerateVoxelreturnsALS(beams, beams_returns, mean_beam, sd_beam, klaser, ctx.pc.transmittance_laser);

    ExportPointcloudHeader(beams, ctx.out.output_pointcloud);

    // now write point cloud records to file
    int nb_beams = int(beams.size()/2);
    int index_return = 0;

    for(int beam = 0; beam < nb_beams; beam++){
        int site_beam = beams[0 + beam * 2];
        int nb_returns_beam = beams[1 + beam * 2];

        int row = site_beam/ctx.grid.cols;
        int col = site_beam%ctx.grid.cols;

        // long format + 0.01 scaling (1cm precision)
        int32_t x_hit = round((float(col) + gsl_rng_uniform(ctx.rng.gslrand)) * 100.0);
        int32_t y_hit = round((float(row) + gsl_rng_uniform(ctx.rng.gslrand)) * 100.0);

        for(int rtrn = 0; rtrn < nb_returns_beam; rtrn++){
            int32_t z_hit = round(beams_returns[index_return] * 100.0);
            index_return++;

            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&x_hit), sizeof(x_hit));
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&y_hit), sizeof(y_hit));
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&z_hit), sizeof(z_hit));

            uint16_t intensity = 0;
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&intensity), sizeof(intensity));

            // the LAS-format makes use of sub-byte level information for return number / number of returns / scan direction / edge of flight line
            // all are together in one single byte, stored in, respectively, 3 bits, 3 bits, 1 bit, 1 bit
            // we need to reconstruct this here
            // cf. comment in LAS defintion: "[A note on Bit Fields – The LAS storage format is “Little Endian.” This means that multi-byte data fields are stored in memory from least significant byte at the low address to most significant byte at the high address. Bit fields are always interpreted as bit 0 set to 1 equals 1, bit 1 set to 1 equals 2, bit 2 set to 1 equals 4 and so forth.]"; here bit 0 seems to mean the bit with the exponent 0, NOT the location/position 0 from left to right
            // essentially (tested through trial and error), the system is:
            // bit location 0 (exponent 7, "bit 7") = edge of flight line
            // bit location 1 (exponent 6, "bit 6") = scan direction flag
            // bit location 2-4 (exponents 5-3) = number of returns
            // bit location 5-7 (exponents 2-0) = return number

            // we only have 5 cases for both return number and number of returns (both max 5), which are:
            // 1: 001
            // 2: 010
            // 3: 011
            // 4: 100
            // 5: 101

            vector<int> information_bitlevel;
            information_bitlevel.reserve(8);

            // we need to start with bits 7 and 6
            information_bitlevel.push_back(0);  // edge of flight line (no such thing here)
            information_bitlevel.push_back(0);  // scan direction (unimportant here)

            // the number of returns
            if(nb_returns_beam <= 3){
                information_bitlevel.push_back(0);
                if(nb_returns_beam == 1){
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(1);
                } else {
                    information_bitlevel.push_back(1);
                    if(nb_returns_beam == 2){
                        information_bitlevel.push_back(0);
                    } else {
                        information_bitlevel.push_back(1);
                    }
                }
            } else {
                information_bitlevel.push_back(1);
                if(nb_returns_beam == 4){ // i.e. 4th return
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(0);
                } else {
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(1);
                }
            }

            // and return number
            int rtrn_actual = rtrn + 1;
            if(rtrn_actual <= 3){
                information_bitlevel.push_back(0);
                if(rtrn_actual == 1){
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(1);
                } else {
                    information_bitlevel.push_back(1);
                    if(rtrn_actual == 2){
                        information_bitlevel.push_back(0);
                    } else {
                        information_bitlevel.push_back(1);
                    }
                }
            } else {
                information_bitlevel.push_back(1);
                if(rtrn_actual == 4){
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(0);
                } else {
                    information_bitlevel.push_back(0);
                    information_bitlevel.push_back(1);
                }
            }

            // now put combined binary information into a single decimal number
            // e.g. first return out of one return total should be: 00 001 001 (which is 9 in decimal notation)
            // e.g. fourth return out of five returns total should be: 00 101 100 (which is 44 in decimal notation)
            int8_t return_info = 0;
            for(int i = 0; i < 8; i++){
                return_info += information_bitlevel[i] * pow(2,7-i);
            }
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&return_info), sizeof(return_info));

            unsigned char classification = 0;
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&classification), sizeof(classification));

            char scan_angle_rank = 0;
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&scan_angle_rank), sizeof(scan_angle_rank));

            unsigned char user_data = 0;
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&user_data), sizeof(user_data));

            uint16_t point_source_id = 1;
            ctx.out.output_pointcloud.write(reinterpret_cast<const char *>(&point_source_id), sizeof(point_source_id));

        }
    }
} */

#ifdef TRACK_INDIVIDUALS
// ##############################################
//  Global function: tree level tracking of key variables
// ##############################################
void TrackingData_andOutput(Context &ctx)
{
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        // we start the accounting the year after the trees have been born
        if (ctx.T[site].t_age > 0)
        {
            int timeofyear_born = ctx.T[site].t_timeofyear_born;

            if (timeofyear_born >= 0)
            {
                // these are the yearly figures, reset to zero every year
                ctx.T[site].t_GPP_sumyear += ctx.T[site].t_GPP;
                ctx.T[site].t_GPPsquared_sumyear += ctx.T[site].t_GPP * ctx.T[site].t_GPP;
                ctx.T[site].t_NPP_sumyear += ctx.T[site].t_NPP;
                ctx.T[site].t_NPPsquared_sumyear += ctx.T[site].t_NPP * ctx.T[site].t_NPP;
                ctx.T[site].t_Rday_sumyear += ctx.T[site].t_Rday;
                ctx.T[site].t_Rnight_sumyear += ctx.T[site].t_Rnight;
                ctx.T[site].t_Rstem_sumyear += ctx.T[site].t_Rstem;
                ctx.T[site].t_carbon_storage_avgyear += ctx.T[site].t_carbon_storage * ctx.time.timestep;
                // these are the whole lifetime cumulated figures, never reset to zero and put out at tree death
                ctx.T[site].t_GPPcum += ctx.T[site].t_GPP;
                ctx.T[site].t_NPPcum += ctx.T[site].t_NPP;
                ctx.T[site].t_GPPsquared_cum += ctx.T[site].t_GPP * ctx.T[site].t_GPP;
                ctx.T[site].t_NPPsquared_cum += ctx.T[site].t_NPP * ctx.T[site].t_NPP;

                float agb = 1000.0 * ctx.T[site].CalcAGB(ctx);

                // write to output every year and then reset to zero
                int timeofyear = GetTimeofyear(ctx);
                if (timeofyear == timeofyear_born)
                {
                    // write to output
                    if (ctx.T[site].t_dbh >= 0.1)
                    {
                        ctx.out.output_track[1] << ctx.T[site].t_site << "\t" << timeofyear_born << "\t" << ctx.time.iter << "\t" << ctx.T[site].t_age << "\t" << ctx.T[site].t_seedsproduced_sumyear << "\t" << ctx.T[site].t_seedsproduced << "\t" << ctx.T[site].t_time_carbonstarvation_year << "\t" << ctx.T[site].t_time_carbonstarvation << "\t" << ctx.T[site].t_dbh << "\t" << ctx.T[site].t_dbh - ctx.T[site].t_dbh_tracked << "\t" << ctx.T[site].t_height << "\t" << ctx.T[site].t_height - ctx.T[site].t_height_tracked << "\t" << ctx.T[site].t_CR << "\t" << ctx.T[site].t_CR - ctx.T[site].t_CR_tracked << "\t" << agb << "\t" << agb - ctx.T[site].t_agb_tracked << "\t" << ctx.T[site].t_GPP_sumyear << "\t" << ctx.T[site].t_GPPsquared_sumyear << "\t" << ctx.T[site].t_NPP_sumyear << "\t" << ctx.T[site].t_NPPsquared_sumyear << "\t" << ctx.T[site].t_Rday_sumyear << "\t" << ctx.T[site].t_Rnight_sumyear << "\t" << ctx.T[site].t_Rstem_sumyear << "\t" << ctx.T[site].t_LAIabove_effavgyear << "\t" << ctx.T[site].t_carbon_storage_avgyear << endl;
                    }
                    // reset
                    ctx.T[site].t_time_carbonstarvation_year = 0;
                    ctx.T[site].t_seedsproduced_sumyear = 0;
                    ctx.T[site].t_GPP_sumyear = 0.0;
                    ctx.T[site].t_GPPsquared_sumyear = 0.0;
                    ctx.T[site].t_NPP_sumyear = 0.0;
                    ctx.T[site].t_NPPsquared_sumyear = 0.0;
                    ctx.T[site].t_Rday_sumyear = 0.0;
                    ctx.T[site].t_Rnight_sumyear = 0.0;
                    ctx.T[site].t_Rstem_sumyear = 0.0;
                    ctx.T[site].t_LAIabove_effavgyear = 0.0;
                    ctx.T[site].t_carbon_storage_avgyear = 0.0;
                    ctx.T[site].t_dbh_tracked = ctx.T[site].t_dbh;
                    ctx.T[site].t_height_tracked = ctx.T[site].t_height;
                    ctx.T[site].t_CR_tracked = ctx.T[site].t_CR;
                    ctx.T[site].t_agb_tracked = agb;
                }
            }
        }
    }
}
#endif

#ifdef Output_ABC
// ##############################################
//  Global ABC function: yearly statistics
// ##############################################
//! - Statistics averaged over one year for ten years and then used for the point estimates in the OutputABC(...) routine
//! - !!!: TODO needs to be checked for different timesteps than 1 month
void UpdateMovingAveragesABC(Context &ctx)
{

    int nbtrees_abc = 0;
    float inbhectares_abc = 10000.0 / float(sites_abc);
    float GPP_abc = 0.0, litter_abc = 0.0;

    for (int row = row_start; row < row_end; row++)
    {
        for (int col = col_start; col < col_end; col++)
        {
            int site = col + row * ctx.grid.cols;
            GPP_abc += ctx.T[site].t_GPP * 1.0e-6;
            litter_abc += ctx.T[site].t_litter * 1.0e-6;
            if (ctx.T[site].t_dbh >= 0.1)
                nbtrees_abc++;
        }
    }

    GPP_abc *= inbhectares_abc;
    litter_abc *= inbhectares_abc;

    float mortality_abc;
    float treefall_abc;
    if (nbtrees_abc > 0)
    {
        mortality_abc = float(nbdead_n10_abc) / float(nbtrees_abc);
        treefall_abc = float(nbTreefall10_abc) / float(nbtrees_abc);
    }
    else
    {
        mortality_abc = 0.0;
        treefall_abc = 0.0;
    }
    if (ctx.time.iter < 120)
    {
        GPP_MA[ctx.time.iter] = GPP_abc;
        Litterfall_MA[ctx.time.iter] = litter_abc;
        Mortality_MA[ctx.time.iter] = mortality_abc;
        Treefall_MA[ctx.time.iter] = treefall_abc;
    }
    else
    {
        // move values down
        for (int i = 0; i < 119; i++)
        {
            GPP_MA[i] = GPP_MA[i + 1];
            Litterfall_MA[i] = Litterfall_MA[i + 1];
            Mortality_MA[i] = Mortality_MA[i + 1];
            Treefall_MA[i] = Treefall_MA[i + 1];
        }
        // new value at the top
        GPP_MA[119] = GPP_abc;
        Litterfall_MA[119] = litter_abc;
        Mortality_MA[119] = mortality_abc;
        Treefall_MA[119] = treefall_abc;
    }
}

// ##############################################
//  Global ABC function: update DBH function for ABC routines
// ##############################################
void UpdateDBHtrackingABC(Context &ctx)
{
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
                ctx.T[s].t_dbh_previous = ctx.T[s].t_dbh;
        }
    }
}

// ##############################################
//  Global function: calculating the TROLL transmittance field from simulated LiDAR
// ##############################################
//! - First draws from a distribution to calculate the sampling density, i.e. the number of beams per voxel column
//! - Then loop over voxel column from top to bottom and calculate the number of hits given the density of the respective voxel
//! - If beams are all extinct, NAs are returned (-1), if beams hit ground, all produce guaranteed returns
//! - The parameters used is the k assumed for the laser, which is only based on leaf geometry and should thus be taken equal to ctx.params.klight
//! - The difference to ctx.params.kpar is accounted for by parameterizing the transmittance of leaves, which, in the NIR spectrum, is much larger than for visible light (0.4 vs. 0.1)
//! - Results are saved in transmittance_simulatedALS_sampling for the number of beams, transmittance_simulatedALS for the transmittance
void UpdateTransmittanceCHM_ABC(float mean_beam, float sd_beam, float klaser, float ctx.pc.transmittance_laser)
{
    // loop over the ctx.field.LAI3D field
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int site = c + r * ctx.grid.cols;
            int nbbeams = int(mean_beam + gsl_ran_gaussian(ctx.rng.gslrand, sd_beam)); // always rounding up
            nbbeams = max(nbbeams, 1);

            // loop over the field from maximum height to 0 and iteratively update voxels from top to bottom, following the beam. An alternative version, also allowing for ground returns, can be activated to extending the loop to h >= -1. In this case, when a beam is not extinguished before it reaches the ground (h >= 0), then it is counted as a ground return
            for (int h = ctx.grid.HEIGHT - 1; h >= 0; h--)
            {
                // for(int h = ctx.grid.HEIGHT - 1; h >= -1; h--){
                int hits;
                float transmittance;
                transmittance_simulatedALS_sampling[h][site] = nbbeams;
                if (nbbeams == 0)
                {
                    // If there is no beam reaching the voxel, transmittance and hits are set to NA (i.e. -1.0)
                    hits = 0;
                    transmittance = -1.0;
                }
                else
                {
                    if (h >= 0)
                    {
                        // returns due to vegetation
                        float LAI_above = ctx.field.LAI3D[h + 1][site + ctx.grid.SBORD];
                        float LAI_current = ctx.field.LAI3D[h][site + ctx.grid.SBORD];

                        float prob_hit;
                        if (LAI_above == 100.0 & LAI_current == 100.0)
                        {
                            // stem returns
                            hits = nbbeams;
                            nbbeams = 0;
                            transmittance = 0.0;
                        }
                        else
                        {
                            // leaf/twig returns
                            float LAD = LAI_current - LAI_above;
                            if (LAD > 0.0)
                                prob_hit = 1.0 - exp(-klaser * LAD);
                            else
                                prob_hit = 0.0;
                            hits = gsl_ran_binomial(ctx.rng.gslrand, prob_hit, nbbeams);
                            // transmittance = exp(-klaser * LAD);
                            if (hits == 0)
                            {
                                transmittance = 1.0;
                            }
                            else
                            {
                                transmittance = float(nbbeams - hits) / float(nbbeams);
                                nbbeams -= hits;
                                // nbbeams += int(0.1*float(hits));                            // 10% of intercepted beams are not getting extinct
                                //  now simulate transmittance of beam through the leaves
                                int hits_notextinct = gsl_ran_binomial(ctx.rng.gslrand, ctx.pc.transmittance_laser, hits);
                                nbbeams += hits_notextinct;
                            }
                        }
                    }
                    else
                    {
                        // ground returns
                        hits = nbbeams;
                        // nbbeams = 0;
                        transmittance = -1.0;
                    }
                }
                transmittance_simulatedALS[h][site] = transmittance;
            }
        }
    }
    // Also output the "direct"/"actual" transmittance of each voxel, i.e. simply based on inversing the Beer Lambert law and estimating transmittance from the leaf area density
    // This can be used to compare actual vs. lidar-derived transmittance estimates
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int s = c + r * ctx.grid.cols;
                float LAD = 0.0;
                if (h < ctx.grid.HEIGHT)
                    LAD = ctx.field.LAI3D[h][s + ctx.grid.SBORD] - ctx.field.LAI3D[h + 1][s + ctx.grid.SBORD];
                transmittance_direct[h][s] = exp(-ctx.params.klight * LAD);
            }
        }
    }

    // Now calculate CHM fields, both with and without ALS simulation
    // Set previous CHM fields
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            chm_field_previous[s] = chm_field_current[s];
            chm_field_previous_ALS[s] = chm_field_current_ALS[s];
            chm_field_current[s] = 0;
            chm_field_current_ALS[s] = 0;
        }
    }

#ifdef CHM_SPIKEFREE
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
            {
#ifdef CROWN_UMBRELLA
                float height = ctx.T[s].t_height;
                float CR = ctx.T[s].t_CR;
                float CD = ctx.T[s].t_CD;
                float fraction_filled_target = 1.0; // to remove spikes, we assume a fully filled crown for all trees
                int shell_fromtop = 0;              // toplayer
                float noinput = 0.0;

                LoopLayerUpdateCrownStatistic_template(ctx, r, c, height, CR, CD, fraction_filled_target, shell_fromtop, [](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); }, noinput, chm_field_current, KeepFloatAsIs, UpdateCHM);
#else
                int crown_top = int(ctx.T[s].t_height);

                float crown_area = PI * ctx.T[s].t_CR * ctx.T[s].t_CR; // floor of crown_area to bound area accumulation
                int crown_intarea = int(crown_area);           // floor of crown_area to bound area accumulation
                crown_intarea = max(crown_intarea, 1);         // minimum area of crown (1)
                crown_intarea = min(crown_intarea, 1963);      // maximum area of crown (radius 25), int(3.14*25*25)

                int row_crowncenter = s / ctx.grid.cols;
                int col_crowncenter = s % ctx.grid.cols;

                for (int i = 0; i < crown_intarea; i++)
                {
                    int site_relative = ctx.lookup.LookUp_Crown_site[i];
                    int row = row_crowncenter + site_relative / 51 - 25;
                    int col = col_crowncenter + site_relative % 51 - 25;
                    if (row >= 0 && row < ctx.grid.rows && col >= 0 && col < ctx.grid.cols)
                    {
                        int site = col + row * ctx.grid.cols;
                        if (chm_field_current[site] < crown_top)
                            chm_field_current[site] = crown_top;
                    }
                }
#endif
            }
        }
    }
#else

    // calculate current CHM fields
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_canopy = 0;
            for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            {
                if (ctx.field.LAI3D[h][s + ctx.grid.SBORD] > 0.0)
                    height_canopy = max(h, height_canopy);
            }
            chm_field_current[s] = height_canopy;
        }
    }
#endif

    // calculate CHM based on simulated lidar
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_canopy_ALS = 0;
            for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            {
                if (transmittance_simulatedALS[h][s] >= 0.0 && transmittance_simulatedALS[h][s] < 1.0)
                    height_canopy_ALS = max(h, height_canopy_ALS);
            }
            chm_field_current_ALS[s] = height_canopy_ALS;
        }
    }

    // compute CHM changes, both for full simulation and simulated ALS
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        chm_field_changes[s] = chm_field_previous[s] - chm_field_current[s];
        chm_field_changes_ALS[s] = chm_field_previous_ALS[s] - chm_field_current_ALS[s];
    }
}

// ##############################################
//  Global ABC function: output general ABC statistics
// ##############################################
void OutputABC(Context &ctx)
{
    cout << " ABC: Conservation of Traits " << endl;
    OutputABCConservationTraits(ctx, ctx.out.output[11]);
    cout << " ABC: Ground data " << endl;
    OutputABC_ground(ctx, ctx.out.output[12]);
    cout << " ABC: CHM simulation " << endl;
    OutputABC_CHM(ctx, ctx.out.output[13], ctx.out.output[14], ctx.out.output[19]);
    cout << " ABC: Transmittance simulation " << endl;
    OutputABC_transmittance(ctx, ctx.out.output[15], ctx.out.output[16]);
    cout << " ABC: Species outputs " << endl;
    OutputABC_species(ctx, ctx.out.output[23], ctx.out.output[24], ctx.out.output[25], ctx.out.output[26], ctx.out.output[27]);
}

// ##############################################
//  Global ABC function: write headers for ABC outputs
// ##############################################
void OutputABCWriteHeaders(fstream &output_traitconservation, fstream &output_field, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_transmittance, fstream &output_transmittance_ALS, fstream &output_LAIfield, fstream &output_LAIfield_ALS, fstream &output_chmpotential, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass)
{
    // Write headers for trait conservation metrics
    output_traitconservation << "Iter\t" << "mean_ran\t" << "sd_ran\t" << "Height_output\t" << "Heightsd_output\t" << "CR_output\t" << "CRsd_output\t" << "CD_output\t" << "CDsd_output\t" << "P_output\t" << "Psd_output\t" << "N_output\t" << "Nsd_output\t" << "LMA_output\t" << "LMAsd_output\t" << "wsg_output\t" << "wsgsd_output\t" << "dmax_output\t" << "dmaxsd_output\t" << "Height_input\t" << "Heightsd_input\t" << "CR_input\t" << "CRsd_input\t" << "CD_input\t" << "CDsd_input\t" << "P_input\t" << "Psd_input\t" << "N_input\t" << "Nsd_input\t" << "LMA_input\t" << "LMAsd_input\t" << "wsg_input\t" << "wsgsd_input\t" << "dmax_input\t" << "dmaxsd_input" << endl;
    // Write headers for ground metrics
    output_field << "Iter\t" << "Nbsites\t" << "NBspecies_realized10\t" << "NBspecies\t" << "NBspecies10\t" << "Shannon\t" << "Shannon10\t" << "Simpson\t" << "Simpson10\t" << "Abu\t" << "Abu10\t" << "Abu30\t" << "Abu10_retained\t" << "Abu30_retained\t" << "AGB\t" << "AGB10\t" << "BA\t" << "BA10\t" << "LoreyH\t" << "LoreyH10\t" << "YearlyGPP\t" << "YearlyLitterfall\t" << "YearlyMortality\t" << "YearlyTreefall\t" << "mean_LMA\t" << "mean_Nmass\t" << "mean_Pmass\t" << "mean_wsg\t" << "mean_CR\t" << "mean_LMA10\t" << "mean_Nmass10\t" << "mean_Pmass10\t" << "mean_wsg10\t" << "mean_CR10\t" << "DBH_mean\t" << "DBH_sd";

    for (int d = 0; d < 50; d++)
        output_field << "\tDBH_hist";
    output_field << "\tDBHgrowth_mean" << "\tDBHgrowth_sd";
    for (int d = 0; d < 50; d++)
        output_field << "\tDBHgrowth_hist";
    output_field << "\tDBHgrowth_yearly_mean" << "\tDBHgrowth_yearly_sd";
    for (int d = 0; d < 50; d++)
        output_field << "\tDBHgrowth_yearly_hist";
    for (int d = -50; d < 0; d++)
        output_field << "\tDBHshrinkage_hist_werr";
    for (int d = 0; d < 50; d++)
        output_field << "\tDBHgrowth_hist_werr";
    for (int d = -50; d < 0; d++)
        output_field << "\tDBHshrinkage_yearly_hist_werr";
    for (int d = 0; d < 50; d++)
        output_field << "\tDBHgrowth_yearly_hist_werr";
    output_field << endl;

    // Write second row, giving additional information on binning etc.
    output_field << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t";
    for (int d = 0; d < 50; d++)
        output_field << "\t" << float(d) / 20.0;
    output_field << "\t" << "\t";
    for (int d = 0; d < 50; d++)
        output_field << "\t" << float(d) / 250.0;
    output_field << "\t" << "\t";
    for (int d = 0; d < 50; d++)
        output_field << "\t" << float(d) / 1000.0;
    for (int d = -50; d < 0; d++)
        output_field << "\t" << float(d) / 250.0;
    for (int d = 0; d < 50; d++)
        output_field << "\t" << float(d) / 250.0;
    for (int d = -50; d < 0; d++)
        output_field << "\t" << float(d) / 1000.0;
    for (int d = 0; d < 50; d++)
        output_field << "\t" << float(d) / 1000.0;
    output_field << endl;

    // Write headers for chm metrics
    output_CHM << "Iter\t" << "Sites_abc\t" << "CHM_mean\t" << "CHM_sd\t";
    for (int h = 0; h < 70; h++)
        output_CHM << "CHM_hist\t";
    output_CHM << "CHMchange_mean\t" << "CHMchange_sd\t";
    for (int h = 0; h < 140; h++)
        output_CHM << "CHMchange_hist\t";
    output_CHM << endl;

    // Second line, information on bin size etc.
    output_CHM << "\t" << "\t" << "\t" << "\t";
    for (int h = 0; h < 70; h++)
        output_CHM << h << "\t";
    output_CHM << "\t" << "\t";
    for (int h = 0; h < 140; h++)
        output_CHM << h << "\t";
    output_CHM << endl;

    // Same for ALS derived metrics
    output_CHM_ALS << "Iter\t" << "Sites_abc\t" << "CHM_mean\t" << "CHM_sd\t";
    for (int h = 0; h < 70; h++)
        output_CHM_ALS << "CHM_hist\t";
    output_CHM_ALS << "CHMchange_mean\t" << "CHMchange_sd\t";
    for (int h = 0; h < 140; h++)
        output_CHM_ALS << "CHMchange_hist\t";
    output_CHM_ALS << endl;

    // Second line, information on bin size etc.
    output_CHM_ALS << "\t" << "\t" << "\t" << "\t";
    for (int h = 0; h < 70; h++)
        output_CHM_ALS << h << "\t";
    output_CHM_ALS << "\t" << "\t";
    for (int h = 0; h < 140; h++)
        output_CHM_ALS << h << "\t";
    output_CHM_ALS << endl;

    // Write headers for transmittance metrics
    output_transmittance << "Iter\t" << "Sites_abc\t";
    output_transmittance << "maxheight\t" << "voxcanopy_total\t" << "voxmaxheight_total\t" << "voxcanopy_filled_total\t" << "voxcanopy_greater2\t" << "voxmaxheight_greater2\t" << "voxcanopy_filled_greater2\t" << "voxcrown_total\t" << "voxcrown_greater2\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "voxcrown\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "transm_nogaps_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "transm_incanopy_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "transm_full_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "voxcanopy_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "voxcanopy_na_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "voxcanopy_empty_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << "voxcanopy_filled_z\t";

    for (int d = 0; d < 70; d++)
        output_transmittance << "transm_nogaps_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << "transm_full_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << "voxna_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << "voxempty_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << "voxfilled_d\t";

    output_transmittance << "nb_aggregates\t";
    for (int bin = 0; bin < 80; bin++)
        output_transmittance << "LAIhist\t";
    output_transmittance << endl;
    output_transmittance << "\t" << "\t";
    output_transmittance << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t";
    for (int i = 0; i < 8; i++)
        for (int h = 0; h < 70; h++)
            output_transmittance << h << "\t";
    for (int i = 0; i < 5; i++)
        for (int d = 0; d < 70; d++)
            output_transmittance << d << "\t";
    output_transmittance << "\t";
    for (int bin = 0; bin < 80; bin++)
        output_transmittance << bin << "\t";
    output_transmittance << endl;
    output_transmittance_ALS << "Iter\t" << "Sites_abc\t";

    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "transm_nogaps_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "transm_incanopy_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "transm_full_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "voxcanopy_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "voxcanopy_na_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "voxcanopy_empty_z\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << "voxcanopy_filled_z\t";

    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << "transm_nogaps_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << "transm_full_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << "voxna_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << "voxempty_d\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << "voxfilled_d\t";

    output_transmittance_ALS << "nb_aggregates\t";
    for (int bin = 0; bin < 80; bin++)
        output_transmittance_ALS << "LAIhist\t";
    output_transmittance_ALS << endl;
    output_transmittance_ALS << "\t" << "\t";
    for (int i = 0; i < 7; i++)
        for (int h = 0; h < 70; h++)
            output_transmittance_ALS << h << "\t";
    for (int i = 0; i < 5; i++)
        for (int d = 0; d < 70; d++)
            output_transmittance_ALS << d << "\t";
    output_transmittance_ALS << "\t";
    for (int bin = 0; bin < 80; bin++)
        output_transmittance_ALS << bin << "\t";
    output_transmittance_ALS << endl;

    // Now the species headers
    output_species << "Iter";
    output_species10 << "Iter";
    //    output_traits << "Iter\tTrait";
    //    output_traits10 << "Iter\tTrait";
    output_biomass << "Iter";

    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
    {
        output_species << "\t" << ctx.S[spp].s_name;
        output_species10 << "\t" << ctx.S[spp].s_name;
        //        output_traits << "\t" << ctx.S[spp].s_name;
        //        output_traits10 << "\t" << ctx.S[spp].s_name;
        output_biomass << "\t" << ctx.S[spp].s_name;
    }

    output_species << endl;
    output_species10 << endl;
    //    output_traits << endl;
    //    output_traits10 << endl;
    output_biomass << endl;
}

// ##############################################
//  Global ABC function: assesses the conservation of traits between input and output in TROLL
// ##############################################
//! - for example, it may be that a particular species is more successful than expected from its relative frequency at the parameterized site, and thus its abundance will exceed the real abundance. This could then lead to important shifts in the overall community trait patterns between the trait pattern that is provided to TROLL and the one that is found in the simulated community
//! - There is a potentially important ecological question behind it, namely: to what extent do the trait distributions that we observe for seeds/seedlings/saplings (i.e. input from seedrain) differ from those of the mature community
void OutputABCConservationTraits(Context &ctx, fstream &output_traitconservation)
{

    // Trait input/output variation
    int nb_trees_counted = 0;

    // Summary statistics for initial, random configuration (i.e. randomly drawn variance around trait means)
    float mu_random = 0.0, moment2_random = 0.0;

    float mu_height_varinput = 0.0, mu_CR_varinput = 0.0, mu_CD_varinput = 0.0, mu_P_varinput = 0.0, mu_N_varinput = 0.0, mu_LMA_varinput = 0.0, mu_wsg_varinput = 0.0, mu_dbhmax_varinput = 0.0;
    float mu_height_varoutput = 0.0, mu_CR_varoutput = 0.0, mu_CD_varoutput = 0.0, mu_P_varoutput = 0.0, mu_N_varoutput = 0.0, mu_LMA_varoutput = 0.0, mu_wsg_varoutput = 0.0, mu_dbhmax_varoutput = 0.0;

    float moment2_height_varinput = 0.0, moment2_CR_varinput = 0.0, moment2_CD_varinput = 0.0, moment2_P_varinput = 0.0, moment2_N_varinput = 0.0, moment2_LMA_varinput = 0.0, moment2_wsg_varinput = 0.0, moment2_dbhmax_varinput = 0.0;
    float moment2_height_varoutput = 0.0, moment2_CR_varoutput = 0.0, moment2_CD_varoutput = 0.0, moment2_P_varoutput = 0.0, moment2_N_varoutput = 0.0, moment2_LMA_varoutput = 0.0, moment2_wsg_varoutput = 0.0, moment2_dbhmax_varoutput = 0.0;

    // Means and 2nd moment for calculation of sd later on
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        if (ctx.T[s].t_age > 0 && ctx.T[s].t_dbh >= 0.1)
        {
            int dev_rand = int(gsl_rng_uniform_int(ctx.rng.gslrand, 10000)); // modified FF, v.3.1.5

            nb_trees_counted++;
            mu_random += dev_rand;

            mu_height_varinput += log(ctx.intra.d_intraspecific_height[dev_rand]);
            mu_CR_varinput += log(ctx.intra.d_intraspecific_CR[dev_rand]);
            mu_CD_varinput += log(ctx.intra.d_intraspecific_CD[dev_rand]);
            mu_P_varinput += log(ctx.intra.d_intraspecific_P[dev_rand]);
            mu_N_varinput += log(ctx.intra.d_intraspecific_N[dev_rand]);
            mu_LMA_varinput += log(ctx.intra.d_intraspecific_LMA[dev_rand]);
            mu_dbhmax_varinput += log(ctx.intra.d_intraspecific_dbhmax[dev_rand]);
            mu_wsg_varinput += log(ctx.intra.d_intraspecific_wsg[dev_rand]);

            mu_height_varoutput += log(ctx.T[s].t_mult_height);
            mu_CR_varoutput += log(ctx.T[s].t_mult_CR);
            mu_CD_varoutput += log(ctx.T[s].t_mult_CD);
            mu_P_varoutput += log(ctx.T[s].t_mult_P);
            mu_N_varoutput += log(ctx.T[s].t_mult_N);
            mu_LMA_varoutput += log(ctx.T[s].t_mult_LMA);
            mu_dbhmax_varoutput += log(ctx.T[s].t_mult_dbhmax);
            mu_wsg_varoutput += log(ctx.T[s].t_dev_wsg);

            moment2_random += dev_rand;
            moment2_height_varinput += log(ctx.intra.d_intraspecific_height[dev_rand]) * log(ctx.intra.d_intraspecific_height[dev_rand]);
            moment2_CR_varinput += log(ctx.intra.d_intraspecific_CR[dev_rand]) * log(ctx.intra.d_intraspecific_CR[dev_rand]);
            moment2_CD_varinput += log(ctx.intra.d_intraspecific_CD[dev_rand]) * log(ctx.intra.d_intraspecific_CD[dev_rand]);
            moment2_P_varinput += log(ctx.intra.d_intraspecific_P[dev_rand]) * log(ctx.intra.d_intraspecific_P[dev_rand]);
            moment2_N_varinput += log(ctx.intra.d_intraspecific_N[dev_rand]) * log(ctx.intra.d_intraspecific_N[dev_rand]);
            moment2_LMA_varinput += log(ctx.intra.d_intraspecific_LMA[dev_rand]) * log(ctx.intra.d_intraspecific_LMA[dev_rand]);
            moment2_dbhmax_varinput += log(ctx.intra.d_intraspecific_dbhmax[dev_rand]) * log(ctx.intra.d_intraspecific_dbhmax[dev_rand]);
            moment2_wsg_varinput += log(ctx.intra.d_intraspecific_wsg[dev_rand]) * log(ctx.intra.d_intraspecific_wsg[dev_rand]);

            moment2_height_varoutput += log(ctx.T[s].t_mult_height) * log(ctx.T[s].t_mult_height);
            moment2_CR_varoutput += log(ctx.T[s].t_mult_CR) * log(ctx.T[s].t_mult_CR);
            moment2_CD_varoutput += log(ctx.T[s].t_mult_CD) * log(ctx.T[s].t_mult_CD);
            moment2_P_varoutput += log(ctx.T[s].t_mult_P) * log(ctx.T[s].t_mult_P);
            moment2_N_varoutput += log(ctx.T[s].t_mult_N) * log(ctx.T[s].t_mult_N);
            moment2_LMA_varoutput += log(ctx.T[s].t_mult_LMA) * log(ctx.T[s].t_mult_LMA);
            moment2_dbhmax_varoutput += log(ctx.T[s].t_mult_dbhmax) * log(ctx.T[s].t_mult_dbhmax);
            moment2_wsg_varoutput += log(ctx.T[s].t_dev_wsg) * log(ctx.T[s].t_dev_wsg);
        }
    }

    // Calculation of standard deviation
    float sd_random;

    float sd_height_varinput = 0.0, sd_CR_varinput = 0.0, sd_CD_varinput = 0.0, sd_P_varinput = 0.0, sd_N_varinput = 0.0, sd_LMA_varinput = 0.0, sd_wsg_varinput = 0.0, sd_dbhmax_varinput = 0.0;
    float sd_height_varoutput = 0.0, sd_CR_varoutput = 0.0, sd_CD_varoutput = 0.0, sd_P_varoutput = 0.0, sd_N_varoutput = 0.0, sd_LMA_varoutput = 0.0, sd_wsg_varoutput = 0.0, sd_dbhmax_varoutput = 0.0;

    if (nb_trees_counted > 0)
    {
        mu_random *= 1.0 / float(nb_trees_counted);

        mu_height_varinput *= 1.0 / float(nb_trees_counted);
        mu_CR_varinput *= 1.0 / float(nb_trees_counted);
        mu_CD_varinput *= 1.0 / float(nb_trees_counted);
        mu_P_varinput *= 1.0 / float(nb_trees_counted);
        mu_N_varinput *= 1.0 / float(nb_trees_counted);
        mu_LMA_varinput *= 1.0 / float(nb_trees_counted);
        mu_wsg_varinput *= 1.0 / float(nb_trees_counted);
        mu_dbhmax_varinput *= 1.0 / float(nb_trees_counted);

        mu_height_varoutput *= 1.0 / float(nb_trees_counted);
        mu_CR_varoutput *= 1.0 / float(nb_trees_counted);
        mu_CD_varoutput *= 1.0 / float(nb_trees_counted);
        mu_P_varoutput *= 1.0 / float(nb_trees_counted);
        mu_N_varoutput *= 1.0 / float(nb_trees_counted);
        mu_LMA_varoutput *= 1.0 / float(nb_trees_counted);
        mu_dbhmax_varoutput *= 1.0 / float(nb_trees_counted);
        mu_wsg_varoutput *= 1.0 / float(nb_trees_counted);

        sd_random = sqrt(moment2_random / float(nb_trees_counted) - mu_random * mu_random);

        sd_height_varinput = sqrt(moment2_height_varinput / float(nb_trees_counted) - mu_height_varinput * mu_height_varinput);
        sd_CR_varinput = sqrt(moment2_CR_varinput / float(nb_trees_counted) - mu_CR_varinput * mu_CR_varinput);
        sd_CD_varinput = sqrt(moment2_CD_varinput / float(nb_trees_counted) - mu_CD_varinput * mu_CD_varinput);
        sd_P_varinput = sqrt(moment2_P_varinput / float(nb_trees_counted) - mu_P_varinput * mu_P_varinput);
        sd_N_varinput = sqrt(moment2_N_varinput / float(nb_trees_counted) - mu_N_varinput * mu_N_varinput);
        sd_LMA_varinput = sqrt(moment2_LMA_varinput / float(nb_trees_counted) - mu_LMA_varinput * mu_LMA_varinput);
        sd_wsg_varinput = sqrt(moment2_wsg_varinput / float(nb_trees_counted) - mu_wsg_varinput * mu_wsg_varinput);
        sd_dbhmax_varinput = sqrt(moment2_dbhmax_varinput / float(nb_trees_counted) - mu_dbhmax_varinput * mu_dbhmax_varinput);

        sd_height_varoutput = sqrt(moment2_height_varoutput / float(nb_trees_counted) - mu_height_varoutput * mu_height_varoutput);
        sd_CR_varoutput = sqrt(moment2_CR_varoutput / float(nb_trees_counted) - mu_CR_varoutput * mu_CR_varoutput);
        sd_CD_varoutput = sqrt(moment2_CD_varoutput / float(nb_trees_counted) - mu_CD_varoutput * mu_CD_varoutput);
        sd_P_varoutput = sqrt(moment2_P_varoutput / float(nb_trees_counted) - mu_P_varoutput * mu_P_varoutput);
        sd_N_varoutput = sqrt(moment2_N_varoutput / float(nb_trees_counted) - mu_N_varoutput * mu_N_varoutput);
        sd_LMA_varoutput = sqrt(moment2_LMA_varoutput / float(nb_trees_counted) - mu_LMA_varoutput * mu_LMA_varoutput);
        sd_wsg_varoutput = sqrt(moment2_wsg_varoutput / float(nb_trees_counted) - mu_wsg_varoutput * mu_wsg_varoutput);
        sd_dbhmax_varoutput = sqrt(moment2_dbhmax_varoutput / float(nb_trees_counted) - mu_dbhmax_varoutput * mu_dbhmax_varoutput);
    }
    else
    {
        sd_random = 0.0;
    }
    // Write to output file
    output_traitconservation << ctx.time.iter << "\t" << mu_random << "\t" << sd_random << "\t" << mu_height_varoutput << "\t" << sd_height_varoutput << "\t" << mu_CR_varoutput << "\t" << sd_CR_varoutput << "\t" << mu_CD_varoutput << "\t" << sd_CD_varoutput << "\t" << mu_P_varoutput << "\t" << sd_P_varoutput << "\t" << mu_N_varoutput << "\t" << sd_N_varoutput << "\t" << mu_LMA_varoutput << "\t" << sd_LMA_varoutput << "\t" << mu_wsg_varoutput << "\t" << sd_wsg_varoutput << "\t" << mu_dbhmax_varoutput << "\t" << sd_dbhmax_varoutput << "\t" << mu_height_varinput << "\t" << sd_height_varinput << "\t" << mu_CR_varinput << "\t" << sd_CR_varinput << "\t" << mu_CD_varinput << "\t" << sd_CD_varinput << "\t" << mu_P_varinput << "\t" << sd_P_varinput << "\t" << mu_N_varinput << "\t" << sd_N_varinput << "\t" << mu_LMA_varinput << "\t" << sd_LMA_varinput << "\t" << mu_wsg_varinput << "\t" << sd_wsg_varinput << "\t" << mu_dbhmax_varinput << "\t" << sd_dbhmax_varinput << "\t" << endl;
}

// ##############################################
//  Global ABC function: returns ABC ground outputs
// ##############################################
//! - field measured summary statistics: total number of species, abundances > 10cm, >30cm (provided for convenience, information already included in distribution below), AGB
void OutputABC_ground(Context &ctx, fstream &output_field)
{

    // Self-explanatory
    int NBspecies_realized10 = 0;
    int NBspecies = 0, NBspecies10 = 0;
    float Shannon = 0.0, Shannon10 = 0.0, Simpson = 0.0, Simpson10 = 0.0;
    float Abu = 0, Abu10 = 0, Abu30 = 0, Abu10_prev = 0, Abu30_prev = 0, Abu10_retained = 0, Abu30_retained = 0;
    float AGB = 0.0, AGB10 = 0.0, ba = 0.0, ba10 = 0.0, LoreyH = 0.0, LoreyH10 = 0.0;

    float mean_LMA = 0.0, mean_Nmass = 0.0, mean_Pmass = 0.0, mean_wsg = 0.0, mean_CR = 0.0;
    float mean_LMA10 = 0.0, mean_Nmass10 = 0.0, mean_Pmass10 = 0.0, mean_wsg10 = 0.0, mean_CR10 = 0.0;

    // This is just a check whether the initialized species are preserved in the model
    for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
        if (ctx.S[spp].s_dbhmax_realized > 0.1)
            NBspecies_realized10++;
    // and calculate the number of species
    for (int spp = 0; spp < ctx.grid.nbspp + 1; spp++)
    {
        abundances_species[spp] = 0;
        abundances_species10[spp] = 0;
    }
    // Compute summary statistics, but only over the respective area
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
            {
                Abu++;
                float dbh_tree = ctx.T[s].t_dbh;
                float dbh_track = ctx.T[s].t_dbh_previous;
                float height_tree = ctx.T[s].t_height;
                float ba_tree = dbh_tree * dbh_tree * 0.25 * PI;
                float AGBtree;
                AGBtree = 0.001 * ctx.T[s].CalcAGB(); // convert from kg to tons
                AGB += AGBtree;
                LoreyH += ba_tree * height_tree;
                ba += ba_tree;
                mean_LMA += ctx.T[s].t_LMA;
                mean_Nmass += ctx.T[s].t_Nmass;
                mean_Pmass += ctx.T[s].t_Pmass;
                mean_wsg += ctx.T[s].t_wsg;
                mean_CR += ctx.T[s].t_CR;

                int sp_lab = ctx.T[s].t_sp_lab;
                abundances_species[sp_lab]++;

                if (dbh_tree >= 0.1)
                {
                    abundances_species10[sp_lab]++;
                    Abu10++;
                    LoreyH10 += ba_tree * height_tree;
                    ba10 += ba_tree;
                    AGB10 += AGBtree;

                    mean_LMA10 += ctx.T[s].t_LMA;
                    mean_Nmass10 += ctx.T[s].t_Nmass;
                    mean_Pmass10 += ctx.T[s].t_Pmass;
                    mean_wsg10 += ctx.T[s].t_wsg;
                    mean_CR10 += ctx.T[s].t_CR;

                    if (dbh_tree >= 0.3)
                        Abu30++;
                }
                if (dbh_track >= 0.1)
                {
                    Abu10_prev++;
                    if (dbh_track >= 0.3)
                        Abu30_prev++;
                }
            }
        }
    }

    // Divide means by abundances
    float inv_abu, inv_abu10;
    if (Abu > 0)
        inv_abu = 1.0 / Abu;
    else
        inv_abu = 0.0;
    if (Abu10 > 0)
        inv_abu10 = 1.0 / Abu10;
    else
        inv_abu10 = 0.0;

    mean_LMA *= inv_abu;
    mean_Nmass *= inv_abu;
    mean_Pmass *= inv_abu;
    mean_wsg *= inv_abu;
    mean_CR *= inv_abu;

    mean_LMA10 *= inv_abu10;
    mean_Nmass10 *= inv_abu10;
    mean_Pmass10 *= inv_abu10;
    mean_wsg10 *= inv_abu10;
    mean_CR10 *= inv_abu10;

    // Denominator of Lorey's height
    if (ba > 0.0)
        LoreyH *= 1.0 / ba;
    if (ba10 > 0.0)
        LoreyH10 *= 1.0 / ba10;

    // Abundance growth
    Abu10_retained = Abu10 - Abu10_prev;
    Abu30_retained = Abu30 - Abu30_prev;

    // Now calculate the first three orders of diversity
    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
    {
        int abu_spp = abundances_species[spp];
        int abu_spp10 = abundances_species10[spp];

        if (abu_spp > 0)
        {
            float prop_spp = float(abu_spp) / float(Abu);

            NBspecies++;
            Shannon -= prop_spp * log(prop_spp);
            Simpson += prop_spp * prop_spp;
        }

        if (abu_spp10 > 0)
        {
            float prop_spp10 = float(abu_spp10) / float(Abu10);

            NBspecies10++;
            Shannon10 -= prop_spp10 * log(prop_spp10);
            Simpson10 += prop_spp10 * prop_spp10;
        }
    }

    // Standardize output metrics
    float isites_abc = 10000.0 / float(sites_abc);
    AGB *= isites_abc;
    AGB10 *= isites_abc;
    Abu10 *= isites_abc;
    Abu30 *= isites_abc;
    Abu10_retained *= isites_abc;
    Abu30_retained *= isites_abc;
    ba *= isites_abc;
    ba10 *= isites_abc;

    // Now add the process metrics from moving averages
    float YearlyGPP = 0.0, YearlyLitterfall = 0.0, YearlyMortality = 0.0, YearlyTreefall = 0.0;
    float nb_years = 120.0 / float(ctx.time.iterperyear);
    float inb_years = 1.0 / nb_years;

    for (int i = 0; i < 120; i++)
    {
        YearlyGPP += GPP_MA[i];
        YearlyLitterfall += Litterfall_MA[i];
        YearlyMortality += Mortality_MA[i];
        YearlyTreefall += Treefall_MA[i];
    }

    YearlyGPP *= inb_years;
    YearlyLitterfall *= inb_years;
    YearlyMortality *= inb_years;
    YearlyTreefall *= inb_years;

    // Create DBH histogram, mean and sd
    int dbh_abc[50] = {0};
    float mean_dbh = 0.0;
    float sd_dbh = 0.0;
    float icount_dbh = 0.0;

    if (Abu10 > 1)
    {
        icount_dbh = 1.0 / float(Abu10);

        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int s = c + r * ctx.grid.cols;
                float dbh_tree = ctx.T[s].t_dbh;
                if (ctx.T[s].t_age > 0 && dbh_tree >= 0.1)
                {
                    mean_dbh += dbh_tree;
                    int dbh_bin = int(dbh_tree * 20.0);
                    if (dbh_bin > 49)
                        dbh_bin = 49;
                    dbh_abc[dbh_bin]++;
                }
            }
        }
        mean_dbh *= icount_dbh;
        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int s = c + r * ctx.grid.cols;
                float dbh_tree = ctx.T[s].t_dbh;
                if (ctx.T[s].t_age > 0 && dbh_tree >= 0.1)
                {
                    sd_dbh += (dbh_tree - mean_dbh) * (dbh_tree - mean_dbh);
                }
            }
        }
        float icount_dbh_sd = 0.0;
        icount_dbh_sd = 1.0 / float(Abu10 - 1);
        sd_dbh = sqrt(sd_dbh * icount_dbh_sd);
    }

    // Create DBH growth histogram, mean and sd / CHM analogous period and yearly values
    int dbhgrowth_abc[50] = {0};
    int dbhgrowth_abc_yearly[50] = {0};
    float mean_dbhgrowth = 0.0;
    float mean_dbhgrowth_yearly = 0.0;
    float sd_dbhgrowth = 0.0;
    float sd_dbhgrowth_yearly = 0.0;
    float icount_dbhgrowth = 0;

    // The dbhgrowth witherror variables also simulate errors in measurement, taken from Chave et al. 2004 (all given in cm). See also Réjou-Méchain 2017
    //  this can be used to compare distributions with negative diameter growth
    //  this is, however, empirical, and does not account for actual negative diameter growth, i.e. trees shrinking slightly under different environmental conditions etc.
    int dbhgrowth_witherror_abc[50] = {0};
    int dbhgrowth_witherror_abc_yearly[50] = {0};
    int dbhshrinkage_witherror_abc[50] = {0};
    int dbhshrinkage_witherror_abc_yearly[50] = {0};
    int count_dbhgrowth_witherror = 0;
    float icount_dbhgrowth_witherror = 0;

    if (Abu10_prev > 1)
    {
        icount_dbhgrowth = 1.0 / float(Abu10_prev);

        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int s = c + r * ctx.grid.cols;
                float dbh_previous = ctx.T[s].t_dbh_previous;
                // Errors taken from Chave et al. 2004 (all given in cm). See also Réjou-Méchain 2017
                //  slight modification by assuming that measurement errors cannot go beyond a third of the diameter, e.g. for a stem of 10cm not be above or below 3.3 cm, and for a stem of 100cm not be below or above 33cm), also prevents negative diameters
                float sd1 = (0.0062 * dbh_previous * 100.0 + 0.0904) * 0.01;
                float sd2 = 0.0464;

                float prob = gsl_rng_uniform(ctx.rng.gslrand);
                float error;
                if (prob < 0.95)
                    error = gsl_ran_gaussian(ctx.rng.gslrand, sd1);
                else
                    error = gsl_ran_gaussian(ctx.rng.gslrand, sd2);

                float dbh_previous_witherror = fmaxf(dbh_previous + error, 0.66 * dbh_previous);
                dbh_previous_witherror = fminf(dbh_previous_witherror, 1.33 * dbh_previous);

                if (dbh_previous >= 0.1)
                {
                    float dbh = ctx.T[s].t_dbh;
                    float dbh_growth = dbh - dbh_previous;
                    float dbh_growth_yearly = dbh_growth * 12.0 / float(chmchange_iter);

                    mean_dbhgrowth += dbh_growth;
                    mean_dbhgrowth_yearly += dbh_growth_yearly;
                    int dbhgrowth_bin = int(dbh_growth * 250.0);                // with a limit at 50 bins, this implies maximum growth of 20cm (in census intervals of 10 years)
                    int dbhgrowth_bin_yearly = int(dbh_growth_yearly * 1000.0); // with a limit at 50 bins, this implies that there is no further increase above 0.05m, i.e. maximum yearly growth is 5cm

                    if (dbhgrowth_bin >= 50)
                        dbhgrowth_bin = 50 - 1;
                    if (dbhgrowth_bin_yearly >= 50)
                        dbhgrowth_bin_yearly = 50 - 1;
                    dbhgrowth_abc[dbhgrowth_bin]++;
                    dbhgrowth_abc_yearly[dbhgrowth_bin_yearly]++;
                }

                if (dbh_previous_witherror >= 0.1)
                {
                    count_dbhgrowth_witherror++;
                    float dbh = ctx.T[s].t_dbh;
                    // Errors taken from Chave et al. 2004 (all given in cm). See also Réjou-Méchain 2017
                    //  slight modification by assuming that measurement errors cannot go beyond a third of the diameter, e.g. for a stem of 10cm not be above or below 3.3 cm, and for a stem of 100cm not be below or above 33cm), also prevents negative diameters
                    float sd1 = (0.0062 * dbh_previous * 100.0 + 0.0904) * 0.01;
                    float sd2 = 0.0464;

                    float prob = gsl_rng_uniform(ctx.rng.gslrand);
                    float error;
                    if (prob < 0.95)
                        error = gsl_ran_gaussian(ctx.rng.gslrand, sd1);
                    else
                        error = gsl_ran_gaussian(ctx.rng.gslrand, sd2);

                    float dbh_witherror = fmaxf(dbh + error, 0.66 * dbh);
                    dbh_witherror = fminf(dbh_witherror, 1.33 * dbh);

                    float dbh_growth_witherror = dbh_witherror - dbh_previous_witherror;
                    float dbh_growth_yearly_witherror = dbh_growth_witherror * 12.0 / float(chmchange_iter);

                    if (dbh_growth_witherror >= 0.0)
                    {
                        int dbhgrowth_witherror_bin = int(dbh_growth_witherror * 250.0);                // with a limit at 50 bins, this implies maximum growth of 50cm (in census intervals of 10 years)
                        int dbhgrowth_witherror_bin_yearly = int(dbh_growth_yearly_witherror * 1000.0); // with a limit at 50 bins, this implies that there is no further increase above 0.05m, i.e. maximum yearly growth is 5cm

                        if (dbhgrowth_witherror_bin >= 50)
                            dbhgrowth_witherror_bin = 50 - 1;
                        if (dbhgrowth_witherror_bin_yearly >= 50)
                            dbhgrowth_witherror_bin_yearly = 50 - 1;
                        dbhgrowth_witherror_abc[dbhgrowth_witherror_bin]++;
                        dbhgrowth_witherror_abc_yearly[dbhgrowth_witherror_bin_yearly]++;
                    }
                    else
                    {
                        dbh_growth_witherror *= -1.0;
                        dbh_growth_yearly_witherror *= -1.0;

                        int dbhgrowth_witherror_bin = int(dbh_growth_witherror * 250.0);                // with a limit at 50 bins, this implies maximum growth of 50cm (in census intervals of 10 years)
                        int dbhgrowth_witherror_bin_yearly = int(dbh_growth_yearly_witherror * 1000.0); // with a limit at 50 bins, this implies that there is no further increase above 0.05m, i.e. maximum yearly growth is 5cm

                        if (dbhgrowth_witherror_bin >= 50)
                            dbhgrowth_witherror_bin = 50 - 1;
                        if (dbhgrowth_witherror_bin_yearly >= 50)
                            dbhgrowth_witherror_bin_yearly = 50 - 1;
                        dbhshrinkage_witherror_abc[dbhgrowth_witherror_bin]++;
                        dbhshrinkage_witherror_abc_yearly[dbhgrowth_witherror_bin_yearly]++;
                    }
                }
            }
        }

        if (count_dbhgrowth_witherror > 0)
            icount_dbhgrowth_witherror = 1.0 / float(count_dbhgrowth_witherror);
        else
            icount_dbhgrowth_witherror = 0.0;

        mean_dbhgrowth *= icount_dbhgrowth;
        mean_dbhgrowth_yearly *= icount_dbhgrowth;

        // Now standard deviation of dbh growth
        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int s = c + r * ctx.grid.cols;

                float dbh_previous = ctx.T[s].t_dbh_previous;

                if (dbh_previous >= 0.1)
                {
                    float dbh_growth = ctx.T[s].t_dbh - dbh_previous;
                    float dbh_growth_yearly = dbh_growth * 12.0 / float(chmchange_iter);

                    sd_dbhgrowth += (dbh_growth - mean_dbhgrowth) * (dbh_growth - mean_dbhgrowth);
                    sd_dbhgrowth_yearly += (dbh_growth_yearly - mean_dbhgrowth_yearly) * (dbh_growth - mean_dbhgrowth_yearly);
                }
            }
        }

        float icount_dbhgrowth_sd = 0.0;
        icount_dbhgrowth_sd = 1.0 / float(Abu10_prev - 1);

        sd_dbhgrowth = sqrt(sd_dbhgrowth * icount_dbhgrowth_sd);
        sd_dbhgrowth_yearly = sqrt(sd_dbhgrowth_yearly * icount_dbhgrowth_sd);
    }

    // Ripley's K(r), transformed to L estimator through L(r) = sqrt(K/Pi) - r
    //    float itrees_ripley = 0.0;
    //    int area_ripley = ctx.grid.cols * ctx.grid.rows;
    //    int sum_ripley[25] = {0};
    //
    //    if(Abu10 > 0){
    //        int trees_ripley = Abu10;
    //        itrees_ripley = 1.0/trees_ripley;
    //
    //        for(int r=row_start;r<row_end;r++){
    //            for(int c=col_start;c<col_end;c++){
    //                int site = c + r*ctx.grid.cols;
    //                if(ctx.T[site].t_age > 0 && ctx.T[site].t_dbh > 0.1){
    //
    //                    for(int r = 0; r < 25; r++){
    //                        int t = r + 1;
    //                        int t_squared = t * t;
    //                        for(int site_compare = 0; site_compare < ctx.grid.sites; site_compare++){
    //                            if(ctx.T[site_compare].t_age > 0 && ctx.T[site_compare].t_dbh > 0.1 && site != site_compare){
    //                                int row = site/ctx.grid.cols;
    //                                int col = site%ctx.grid.cols;
    //                                int row_compare = site_compare/ctx.grid.cols;
    //                                int col_compare = site_compare%ctx.grid.cols;
    //                                int dist_squared = (row - row_compare) * (row - row_compare) + (col - col_compare) * (col - col_compare);
    //                                if(dist_squared < t_squared){
    //                                    sum_ripley[r]++;
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }

    // Write to file
    output_field << ctx.time.iter << "\t" << sites_abc << "\t" << NBspecies_realized10 << "\t" << NBspecies << "\t" << NBspecies10 << "\t" << Shannon << "\t" << Shannon10 << "\t" << Simpson << "\t" << Simpson10 << "\t" << Abu << "\t" << Abu10 << "\t" << Abu30 << "\t" << Abu10_retained << "\t" << Abu30_retained << "\t" << AGB << "\t" << AGB10 << "\t" << ba << "\t" << ba10 << "\t" << LoreyH << "\t" << LoreyH10 << "\t" << YearlyGPP << "\t" << YearlyLitterfall << "\t" << YearlyMortality << "\t" << YearlyTreefall << "\t" << mean_LMA << "\t" << mean_Nmass << "\t" << mean_Pmass << "\t" << mean_wsg << "\t" << mean_CR << "\t" << mean_LMA10 << "\t" << mean_Nmass10 << "\t" << mean_Pmass10 << "\t" << mean_wsg10 << "\t" << mean_CR10 << "\t" << mean_dbh << "\t" << sd_dbh;

    for (int d = 0; d < 50; d++)
    {
        float density_DBH = float(dbh_abc[d]) * icount_dbh;
        output_field << "\t" << density_DBH;
    }

    // Write to file
    output_field << "\t" << mean_dbhgrowth << "\t" << sd_dbhgrowth;

    for (int d = 0; d < 50; d++)
    {
        float density_dbhgrowth = float(dbhgrowth_abc[d]) * icount_dbhgrowth;
        output_field << "\t" << density_dbhgrowth;
    }

    output_field << "\t" << mean_dbhgrowth_yearly << "\t" << sd_dbhgrowth_yearly;

    for (int d = 0; d < 50; d++)
    {
        float density_dbhgrowth_yearly = float(dbhgrowth_abc_yearly[d]) * icount_dbhgrowth;
        output_field << "\t" << density_dbhgrowth_yearly;
    }

    for (int d = 50 - 1; d >= 0; d--)
    {
        float density_dbhshrinkage_witherror = float(dbhshrinkage_witherror_abc[d]) * icount_dbhgrowth_witherror;
        output_field << "\t" << density_dbhshrinkage_witherror;
    }

    for (int d = 0; d < 50; d++)
    {
        float density_dbhgrowth_witherror = float(dbhgrowth_witherror_abc[d]) * icount_dbhgrowth_witherror;
        output_field << "\t" << density_dbhgrowth_witherror;
    }

    for (int d = 50 - 1; d >= 0; d--)
    {
        float density_dbhshrinkage_yearly_witherror = float(dbhshrinkage_witherror_abc_yearly[d]) * icount_dbhgrowth_witherror;
        output_field << "\t" << density_dbhshrinkage_yearly_witherror;
    }

    for (int d = 0; d < 50; d++)
    {
        float density_dbhgrowth_yearly_witherror = float(dbhgrowth_witherror_abc_yearly[d]) * icount_dbhgrowth_witherror;
        output_field << "\t" << density_dbhgrowth_yearly_witherror;
    }

    //    for(int r = 0; r < 25; r++){
    //        int t = r + 1;
    //        float K_ripley = float(area_ripley) * float(sum_ripley[r]) * itrees_ripley * itrees_ripley;
    //        output_field << sqrt(K_ripley / PI) - t << "\t";
    //    }

    output_field << endl;
}

// ##############################################
//  Global ABC function: returns ABC outputs for species
// ##############################################
void OutputABC_species(Context &ctx, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass)
{

    // Empty the vectors
    for (int spp = 0; spp < ctx.grid.nbspp + 1; spp++)
    {
        abundances_species[spp] = 0;
        abundances_species10[spp] = 0;
        biomass_species[spp] = 0.0;
        for (int trait = 0; trait < 10; trait++)
        {
            traits_species[spp][trait] = 0.0;
            traits_species10[spp][trait] = 0.0;
        }
    }

    // Now calculate
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            if (ctx.T[s].t_age > 0)
            {
                // First abundance, then mean and standard deviation of the variation in the 5 traits (leaf level traits, wood density, crown radius), mean and standard deviation are calculated in one pass from sum of traits and sum of squares of traits
                int sp_lab = ctx.T[s].t_sp_lab;
                abundances_species[sp_lab]++;

                float agb_tree = 0.001 * ctx.T[s].CalcAGB(); // convert from kg to tons
                biomass_species[sp_lab] += agb_tree;

                float LMA = log(ctx.T[s].t_mult_LMA);
                float LMA2 = LMA * LMA;
                float Nmass = log(ctx.T[s].t_mult_N);
                float Nmass2 = Nmass * Nmass;
                float Pmass = log(ctx.T[s].t_mult_P);
                float Pmass2 = Pmass * Pmass;
                float wsg = ctx.T[s].t_dev_wsg;
                float wsg2 = wsg * wsg;
                float crown = log(ctx.T[s].t_mult_CR);
                float crown2 = crown * crown;

                traits_species[sp_lab][0] += LMA;
                traits_species[sp_lab][1] += LMA2;
                traits_species[sp_lab][2] += Nmass;
                traits_species[sp_lab][3] += Nmass2;
                traits_species[sp_lab][4] += Pmass;
                traits_species[sp_lab][5] += Pmass2;
                traits_species[sp_lab][6] += wsg;
                traits_species[sp_lab][7] += wsg2;
                traits_species[sp_lab][8] += crown;
                traits_species[sp_lab][9] += crown2;

                if (ctx.T[s].t_dbh >= 0.1)
                {
                    abundances_species10[sp_lab]++;

                    traits_species10[sp_lab][0] += LMA;
                    traits_species10[sp_lab][1] += LMA2;
                    traits_species10[sp_lab][2] += Nmass;
                    traits_species10[sp_lab][3] += Nmass2;
                    traits_species10[sp_lab][4] += Pmass;
                    traits_species10[sp_lab][5] += Pmass2;
                    traits_species10[sp_lab][6] += wsg;
                    traits_species10[sp_lab][7] += wsg2;
                    traits_species10[sp_lab][8] += crown;
                    traits_species10[sp_lab][9] += crown2;
                }
            }
        }
    }

    // Create a vector of traits
    vector<string> trait_names;

    trait_names.push_back("LMA");
    trait_names.push_back("LMA_sd");
    trait_names.push_back("Nmass");
    trait_names.push_back("Nmass_sd");
    trait_names.push_back("Pmass");
    trait_names.push_back("Pmass_sd");
    trait_names.push_back("wsg");
    trait_names.push_back("wsg_sd");
    trait_names.push_back("crown");
    trait_names.push_back("crown_sd");

    // Now calculate the means and sds properly
    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
    {
        int abu = abundances_species[spp];
        int abu10 = abundances_species10[spp];
        float inv_abu, inv_abu10;
        if (abu > 0)
            inv_abu = 1.0 / float(abu);
        else
            inv_abu = 0.0;
        if (abu10 > 0)
            inv_abu10 = 1.0 / float(abu10);
        else
            inv_abu10 = 0.0;

        // First divide both quantities by the number of individuals
        for (int trait = 0; trait < 10; trait++)
        {
            traits_species[spp][trait] *= inv_abu;
            traits_species10[spp][trait] *= inv_abu10;
        }

        // Then subtract the squared means from the standard deviations
        traits_species[spp][1] -= traits_species[spp][0] * traits_species[spp][0];
        traits_species[spp][3] -= traits_species[spp][2] * traits_species[spp][2];
        traits_species[spp][5] -= traits_species[spp][4] * traits_species[spp][4];
        traits_species[spp][7] -= traits_species[spp][6] * traits_species[spp][6];
        traits_species[spp][9] -= traits_species[spp][8] * traits_species[spp][8];

        traits_species10[spp][1] -= traits_species10[spp][0] * traits_species10[spp][0];
        traits_species10[spp][3] -= traits_species10[spp][2] * traits_species10[spp][2];
        traits_species10[spp][5] -= traits_species10[spp][4] * traits_species10[spp][4];
        traits_species10[spp][7] -= traits_species10[spp][6] * traits_species10[spp][6];
        traits_species10[spp][9] -= traits_species10[spp][8] * traits_species10[spp][8];
    }

    // Now create outputs
    output_species << ctx.time.iter;
    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
        output_species << "\t" << abundances_species[spp];
    output_species << endl;

    output_species10 << ctx.time.iter;
    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
        output_species10 << "\t" << abundances_species10[spp];
    output_species10 << endl;

    output_biomass << ctx.time.iter;
    for (int spp = 1; spp < ctx.grid.nbspp + 1; spp++)
        output_biomass << "\t" << biomass_species[spp];
    output_biomass << endl;

    //    for(int trait = 0; trait < 10; trait++){
    //
    //        output_traits << ctx.time.iter << "\t" << trait_names[trait];
    //        for(int spp = 1; spp < ctx.grid.nbspp+1; spp++) output_traits  << "\t" << traits_species[spp][trait];
    //        output_traits << endl;
    //
    //        output_traits10 << ctx.time.iter << "\t" << trait_names[trait];
    //        for(int spp = 1; spp < ctx.grid.nbspp+1; spp++) output_traits10  << "\t" << traits_species10[spp][trait];
    //        output_traits10 << endl;
    //    }
}

// ##############################################
//  Global function: returns ABC outputs for canopy height model (CHM)
// ##############################################
void OutputABC_CHM(Context &ctx, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_chmpotential)
{
    // Compute CHM changes
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        chm_field_changes[s] = chm_field_previous[s] - chm_field_current[s];
        chm_field_changes_ALS[s] = chm_field_previous_ALS[s] - chm_field_current_ALS[s];
    }

    // Compute CHM distribution
    int chm_abc[70] = {0}, chm_abc_ALS[70] = {0};
    float mean_chm = 0.0, mean_chm_ALS = 0.0;
    float sd_chm = 0.0, sd_chm_ALS = 0.0;

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_chm = chm_field_current[s];
            mean_chm += float(height_chm);
            chm_abc[height_chm]++;
            int height_chm_ALS = chm_field_current_ALS[s];
            mean_chm_ALS += float(height_chm_ALS);
            chm_abc_ALS[height_chm_ALS]++;
        }
    }

    mean_chm *= isites_abc;
    mean_chm_ALS *= isites_abc;

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_chm = chm_field_current[s];
            sd_chm += (height_chm - mean_chm) * (height_chm - mean_chm);
            int height_chm_ALS = chm_field_current_ALS[s];
            sd_chm_ALS += (height_chm_ALS - mean_chm_ALS) * (height_chm_ALS - mean_chm_ALS);
        }
    }

    sd_chm = sqrt(sd_chm * isites_abc);
    sd_chm_ALS = sqrt(sd_chm_ALS * isites_abc);

    // Compute CHM change distribution
    int chmchange_abc[140] = {0}, chmchange_abc_ALS[140] = {0};
    float mean_chmchange = 0.0, mean_chmchange_ALS = 0.0;
    float sd_chmchange = 0.0, sd_chmchange_ALS = 0.0;

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_chmchange = chm_field_changes[s];
            mean_chmchange += float(height_chmchange);
            chmchange_abc[height_chmchange + 70]++;
            int height_chmchange_ALS = chm_field_changes_ALS[s];
            mean_chmchange_ALS += float(height_chmchange_ALS);
            chmchange_abc_ALS[height_chmchange + 70]++;
        }
    }

    mean_chmchange *= isites_abc;
    mean_chmchange_ALS *= isites_abc;

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int s = c + r * ctx.grid.cols;
            int height_chmchange = chm_field_changes[s];
            sd_chmchange += (height_chmchange - mean_chmchange) * (height_chmchange - mean_chmchange);
            int height_chmchange_ALS = chm_field_changes_ALS[s];
            sd_chmchange_ALS += (height_chmchange_ALS - mean_chmchange_ALS) * (height_chmchange_ALS - mean_chmchange_ALS);
        }
    }

    sd_chmchange = sqrt(sd_chmchange * isites_abc);
    sd_chmchange_ALS = sqrt(sd_chmchange * isites_abc);

    //    //! create variogram for both chmchange fields (final and change)
    //
    //    float chm_variogram[50] = {0.0}, chm_variogram_ALS[50] = {0.0};
    //    int chm_variogram_nb[50] = {0};
    //
    //    float chmchange_variogram[50] = {0.0}, chmchange_variogram_ALS[50] = {0.0};
    //    int chmchange_variogram_nb[50] = {0};
    //
    //    for(int r=row_start;r<row_end;r++){
    //        for(int c=col_start;c<col_end;c++){
    //            int height_canopy = chm_field_current[c + r * ctx.grid.cols];
    //            int height_canopychange = chm_field_changes[c + r * ctx.grid.cols];
    //            int height_canopy_ALS = chm_field_current_ALS[c + r * ctx.grid.cols];
    //            int height_canopychange_ALS = chm_field_changes_ALS[c + r * ctx.grid.cols];
    //
    //            for(int r_compare=row_start;r_compare<row_end;r_compare++){
    //                for(int c_compare=col_start;c_compare<col_end;c_compare++){
    //                    //! only compare if you have not already compared
    //                    if(r_compare > r && c_compare > c){
    //                        int height_canopy_compare = chm_field_current[c_compare + r_compare * ctx.grid.cols];
    //                        int height_canopychange_compare = chm_field_changes[c_compare + r_compare * ctx.grid.cols];
    //                        int height_canopy_compare_ALS = chm_field_current_ALS[c_compare + r_compare * ctx.grid.cols];
    //                        int height_canopychange_compare_ALS = chm_field_changes_ALS[c_compare + r_compare * ctx.grid.cols];
    //
    //                        int dist = int(sqrt(float((r_compare - r) * (r_compare - r) + (c_compare - c) * (c_compare - c))));
    //
    //                        //! only compute semivariance for distances < 50
    //                        if(dist < 50){
    //                            chm_variogram[dist] += float((height_canopy_compare - height_canopy) * (height_canopy_compare - height_canopy));
    //                            chm_variogram_ALS[dist] += float((height_canopy_compare_ALS - height_canopy_ALS) * (height_canopy_compare_ALS - height_canopy_ALS));
    //                            chm_variogram_nb[dist]++;
    //
    //                            chmchange_variogram[dist] += float((height_canopychange_compare - height_canopychange) * (height_canopychange_compare - height_canopychange));
    //                            chmchange_variogram_ALS[dist] += float((height_canopychange_compare_ALS - height_canopychange_ALS) * (height_canopychange_compare_ALS - height_canopychange_ALS));
    //                            chmchange_variogram_nb[dist]++;
    //                        }
    //
    //                    }
    //                }
    //            }
    //        }
    //    }

    // Write to file
    output_CHM << ctx.time.iter << "\t" << sites_abc << "\t" << mean_chm << "\t" << sd_chm << "\t";

    for (int h = 0; h < 70; h++)
    {
        float density_CHM = float(chm_abc[h]) * isites_abc;
        output_CHM << density_CHM << "\t";
    }

    output_CHM << mean_chmchange << "\t" << sd_chmchange << "\t";

    for (int h = 0; h < 140; h++)
    {
        float density_chmchange = float(chmchange_abc[h]) * isites_abc;
        output_CHM << density_chmchange << "\t";
    }

    output_CHM << endl;
    output_CHM_ALS << ctx.time.iter << "\t" << sites_abc << "\t" << mean_chm_ALS << "\t" << sd_chm_ALS << "\t";

    for (int h = 0; h < 70; h++)
    {
        float density_CHM = chm_abc_ALS[h] * isites_abc;
        output_CHM_ALS << density_CHM << "\t";
    }
    output_CHM_ALS << mean_chmchange_ALS << "\t" << sd_chmchange_ALS << "\t";

    for (int h = 0; h < 140; h++)
    {
        float density_chmchange = chmchange_abc_ALS[h] * isites_abc;
        output_CHM_ALS << density_chmchange << "\t";
    }
    output_CHM_ALS << endl;
    // In the future: compute canopy extension (i.e. lateral growth, vertical growth etc.) to test model of Kellner & Asner 2014
}

// ##############################################
//  Global ABC function: ABC outputs
// ##############################################
//! - this function creates various transmittance and lidar metrics, from the LAI field and the simulated lidar
//! - the output variables are: avg transmittance per ctx.diag.layer within all voxels of that ctx.diag.layer (i.e. transmittance_full), only within voxels of that ctx.diag.layer that are also inside the canopy (i.e. transmittance_incanopy, voxels below the CHM, if there is a gap with an average canopy height of 3m for example, then voxelsat 10m won't be counted for the transmittance avg at 10m), and transmittance only in filled voxes (i.e. only voxels with leaves, transmittance_nogaps)
//! - additional metrics are the total number of voxels that have been considered for the averages
//! - all metrics will be computed normalized to the ground (suffix z) and normalized to the canopy (suffix d), and both for the actual 3D canopy (no additional suffix) and a simulated lidar (additional suffix ALS)
//! - !!!: TODO, detailed documentation
void OutputABC_transmittance(Context &ctx, fstream &output_transmittance, fstream &output_transmittance_ALS)
{

    // Compute height normalized metrics (z stands for height dimension) from it, both directly and with ALS simulation
    //  filled and empty voxels are only taken within the canopy
    float transmittance_nogaps_z[70] = {0.0}, transmittance_nogaps_zALS[70] = {0.0};
    float transmittance_incanopy_z[70] = {0.0}, transmittance_incanopy_zALS[70] = {0.0};
    float transmittance_full_z[70] = {0.0}, transmittance_full_zALS[70] = {0.0};
    int voxcrown[70] = {0};
    int voxcanopy_z[70] = {0}, voxcanopy_zALS[70] = {0};
    int voxcanopy_empty_z[70] = {0}, voxcanopy_empty_zALS[70] = {0};
    int voxcanopy_filled_z[70] = {0}, voxcanopy_filled_zALS[70] = {0};
    int voxcanopy_na_z[70] = {0}, voxcanopy_na_zALS[70] = {0};

    // Calculate the volume filled by crowns (for packing densities)
    //  to be consistent with other estimates, we do not only use leaf-filled area, but the whole area used up by the crown (even when overlapping)
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        if (ctx.T[site].t_age > 0)
        {
            int row_center = site / ctx.grid.cols;
            int col_center = site % ctx.grid.cols;

            AddCrownVolumeLayer(ctx, row_center, col_center, ctx.T[site].t_height, ctx.T[site].t_CR, ctx.T[site].t_CD, voxcrown);
        }
    }

    // Compute the voxels for heights greater 2 separately, as the lowest ctx.diag.layer in TROLL can be a bit difficult to interpret (lots of seedlings by default)
    int voxcrown_total = voxcrown[0] + voxcrown[1];
    int voxcrown_greater2 = 0;

    for (int h = 2; h < 70; h++)
    {
        voxcrown_total += voxcrown[h];
        voxcrown_greater2 += voxcrown[h];
    }

    // Update statistics based on 3D fields, including the simulated transmittance field
    for (int h = min(70, ctx.grid.HEIGHT) - 1; h >= 0; h--)
    {
        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int site = c + r * ctx.grid.cols;
                int height_canopy = chm_field_current[site];
                int height_canopy_ALS = chm_field_current_ALS[site];

                float transmittance = transmittance_direct[h][site];
                float transmittanceALS = transmittance_simulatedALS[h][site];

                if (h <= height_canopy)
                {
                    // Without ALS simulation, "real"
                    if (transmittance < 0.0)
                    {
                        voxcanopy_na_z[h]++;
                    }
                    else
                    {
                        if (transmittance == 1.0)
                        {
                            voxcanopy_empty_z[h]++;
                        }
                        else
                        {
                            voxcanopy_filled_z[h]++;
                            transmittance_nogaps_z[h] += transmittance;
                        }
                        transmittance_incanopy_z[h] += transmittance;
                        voxcanopy_z[h]++;
                    }
                    // From ALS simulation
                }

                if (h <= height_canopy_ALS)
                {
                    if (transmittanceALS < 0.0)
                    {
                        voxcanopy_na_zALS[h]++;
                    }
                    else
                    {
                        if (transmittanceALS == 1.0)
                        {
                            voxcanopy_empty_zALS[h]++;
                        }
                        else
                        {
                            voxcanopy_filled_zALS[h]++;
                            transmittance_nogaps_zALS[h] += transmittanceALS;
                        }

                        transmittance_incanopy_zALS[h] += transmittanceALS;
                        voxcanopy_zALS[h]++;
                    }
                }
            }
        }

        // Now calculate the mean
        if (voxcanopy_filled_z[h] > 0)
            transmittance_nogaps_z[h] *= 1.0 / float(voxcanopy_filled_z[h]);
        else
            transmittance_nogaps_z[h] = -1.0;
        if ((voxcanopy_filled_z[h] + voxcanopy_empty_z[h]) > 0)
            transmittance_incanopy_z[h] *= 1.0 / float(voxcanopy_filled_z[h] + voxcanopy_empty_z[h]);
        else
            transmittance_incanopy_z[h] = -1.0;

        if (voxcanopy_filled_zALS[h] > 0)
            transmittance_nogaps_zALS[h] *= 1.0 / float(voxcanopy_filled_zALS[h]);
        else
            transmittance_nogaps_zALS[h] = -1.0;
        if ((voxcanopy_filled_zALS[h] + voxcanopy_empty_zALS[h]) > 0)
            transmittance_incanopy_zALS[h] *= 1.0 / float(voxcanopy_filled_zALS[h] + voxcanopy_empty_zALS[h]);
        else
            transmittance_incanopy_zALS[h] = -1.0;

        // Now the full LAI profile, including non-canopy areas
        int voxeltotal = 0;

        for (int r = row_start; r < row_end; r++)
        {
            for (int c = col_start; c < col_end; c++)
            {
                int site = c + r * ctx.grid.cols;

                float transmittance = transmittance_direct[h][site];
                float transmittanceALS = transmittance_simulatedALS[h][site];

                if (transmittance < 0.0)
                    transmittance_full_z[h] += transmittance_incanopy_z[h];
                else
                    transmittance_full_z[h] += transmittance;

                if (transmittanceALS < 0.0)
                    transmittance_full_zALS[h] += transmittance_incanopy_zALS[h];
                else
                    transmittance_full_zALS[h] += transmittanceALS;

                voxeltotal++;
            }
        }

        if (voxeltotal > 0)
        {
            transmittance_full_z[h] *= 1.0 / float(voxeltotal);
            transmittance_full_zALS[h] *= 1.0 / float(voxeltotal);
        }
        else
        {
            transmittance_full_z[h] = -1.0;
            transmittance_full_zALS[h] = -1.0;
        }
    }

    // Now compute metrics normalized with regard to top of canopy (d for depth)
    float transmittance_nogaps_d[70] = {0.0}, transmittance_nogaps_dALS[70] = {0.0};
    float transmittance_full_d[70] = {0.0}, transmittance_full_dALS[70] = {0.0};
    int voxempty_d[70] = {0}, voxempty_dALS[70] = {0};
    int voxfilled_d[70] = {0}, voxfilled_dALS[70] = {0};
    int voxna_d[70] = {0}, voxna_dALS[70] = {0};

    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int site = c + r * ctx.grid.cols;
            int height_canopy = chm_field_current[site];
            int height_canopy_ALS = chm_field_current_ALS[site];

            for (int d = 0; d <= height_canopy; d++)
            {
                float transmittance = transmittance_direct[height_canopy - d][site];

                if (transmittance < 0.0)
                {
                    voxna_d[d]++;
                }
                else
                {
                    if (transmittance == 1.0)
                    {
                        voxempty_d[d]++;
                    }
                    else
                    {
                        voxfilled_d[d]++;
                        transmittance_nogaps_d[d] += transmittance;
                    }
                    transmittance_full_d[d] += transmittance;
                }
            }

            for (int d = 0; d <= height_canopy_ALS; d++)
            {
                float transmittanceALS = transmittance_simulatedALS[height_canopy_ALS - d][site];

                if (transmittanceALS < 0.0)
                {
                    voxna_dALS[d]++;
                }
                else
                {
                    if (transmittanceALS == 1.0)
                    {
                        voxempty_dALS[d]++;
                    }
                    else
                    {
                        voxfilled_dALS[d]++;
                        transmittance_nogaps_dALS[d] += transmittanceALS;
                    }
                    transmittance_full_dALS[d] += transmittanceALS;
                }
            }
        }
    }

    for (int d = 0; d < 70; d++)
    {
        if (voxfilled_d[d] > 0)
            transmittance_nogaps_d[d] *= 1.0 / float(voxfilled_d[d]);
        else
            transmittance_nogaps_d[d] = -1.0;

        int voxnona = voxfilled_d[d] + voxempty_d[d];
        if (voxnona > 0)
            transmittance_full_d[d] *= 1.0 / float(voxnona);
        else
            transmittance_full_d[d] = -1.0;

        if (voxfilled_dALS[d] > 0)
            transmittance_nogaps_dALS[d] *= 1.0 / float(voxfilled_dALS[d]);
        else
            transmittance_nogaps_dALS[d] = -1.0;

        int voxnonaALS = voxfilled_dALS[d] + voxempty_dALS[d];
        if (voxnonaALS > 0)
            transmittance_full_dALS[d] *= 1.0 / float(voxnonaALS);
        else
            transmittance_full_dALS[d] = -1.0;
    }

    // Finally we create PAI distributions, similar to efforts to estimate PAI from lidar simulations (cf. Greg's work)
    //  interpretation outside of an ABC context is a bit problematic, since we have many NA values that need to be imputated (otherwise PAI estimate is biased downwards), and since mean(log(x)) != log(mean(x))
    //  we fill up voxels within canopy with mean canopy ctx.diag.layer value, this decreases variance, but it does so for both empirical and simulated distributions, and the effect can be investigated given that we also have the direct measurements form ctx.field.LAI3D field
    //  furthermore, the canopy height will probably have a stronger effect on the distribution of PAI values than local variation within one ctx.diag.layer
    //  as above, we will calculate the actual distribution and the ALS-inferred one
    int stepsize = 10;
    int PAIhist[80] = {0}, PAIhist_ALS[80] = {0}; // for values going from 0 to 20, anything beyond is put into the 20 bin
    int nb_aggregates = 0;

    for (int r = row_start; r < row_end; r += stepsize)
    {
        for (int c = col_start; c < col_end; c += stepsize)
        {
            float PAIcolumn = 0.0, PAIcolumn_ALS = 0.0;
            nb_aggregates++;

            // Loop over the whole canopy except for 0 height where there is no transmittance value
            for (int h = ctx.grid.HEIGHT; h > 0; h--)
            {
                int nbvoxels = 0;
                float transmittance_avg = 0.0, transmittance_avgALS = 0.0;

                for (int r_avg = r; r_avg < r + stepsize; r_avg++)
                {
                    for (int c_avg = c; c_avg < c + stepsize; c_avg++)
                    {

                        if (r_avg >= row_start && r_avg < row_end && c_avg >= col_start && c_avg < col_end)
                        {
                            int site_avg = c_avg + r_avg * ctx.grid.cols;
                            float transmittance = transmittance_direct[h][site_avg];
                            float transmittanceALS = transmittance_simulatedALS[h][site_avg];

                            if (transmittance < 0.0)
                                transmittance_avg += transmittance_incanopy_z[h];
                            else
                                transmittance_avg += transmittance;
                            if (transmittanceALS < 0.0)
                                transmittance_avgALS += transmittance_incanopy_zALS[h];
                            else
                                transmittance_avgALS += transmittanceALS;

                            nbvoxels++;
                        }
                    }
                }
                if (nbvoxels > 0)
                {
                    if (transmittance_avg >= 0.0)
                    {
                        transmittance_avg *= 1.0 / float(nbvoxels);
                        float PAI = -log(transmittance_avg) / ctx.params.klight;
                        PAIcolumn += PAI;
                    }
                    if (transmittance_avgALS >= 0.0)
                    {
                        transmittance_avgALS *= 1.0 / float(nbvoxels);
                        float PAI_ALS = -log(transmittance_avgALS) / ctx.params.klight;
                        PAIcolumn_ALS += PAI_ALS;
                    }
                }
            }
            // We bin in bins of 0.25 LAI difference
            int PAIbin = min(int(4.0 * PAIcolumn), 79);
            int PAIbin_ALS = min(int(4.0 * PAIcolumn_ALS), 79);
            PAIhist[PAIbin]++;
            PAIhist_ALS[PAIbin_ALS]++;
        }
    }

    // Calculate summary statistics for crown packing analysis, i.e. the crown volume per voxel, average per ctx.diag.layer and canopy height
    //  first compute maxheight
    int height_max = 0;
    for (int r = row_start; r < row_end; r++)
    {
        for (int c = col_start; c < col_end; c++)
        {
            int site = c + r * ctx.grid.cols;
            int height_canopy = chm_field_current[site];
            if (height_max < height_canopy)
                height_max = height_canopy;
        }
    }

    // number of voxels to consider
    int voxmaxheight_total = height_max * sites_abc, voxmaxheight_greater2 = (height_max - 2) * sites_abc;
    int voxcanopy_total = 0, voxcanopy_greater2 = 0;
    int voxcanopy_filled_total = 0, voxcanopy_filled_greater2 = 0;

    voxcanopy_total += voxcanopy_z[0] + voxcanopy_z[1];
    voxcanopy_filled_total += voxcanopy_filled_z[0] + voxcanopy_filled_z[1];

    for (int h = 2; h < 70; h++)
    {
        voxcanopy_total += voxcanopy_z[h];
        voxcanopy_greater2 += voxcanopy_z[h];
        voxcanopy_filled_total += voxcanopy_filled_z[h];
        voxcanopy_filled_greater2 += voxcanopy_filled_z[h];
    }

    // Write to file
    output_transmittance << ctx.time.iter << "\t" << sites_abc << "\t";

    // Summary statistics
    output_transmittance << height_max << "\t" << voxcanopy_total << "\t" << voxmaxheight_total << "\t" << voxcanopy_filled_total << "\t" << voxcanopy_greater2 << "\t" << voxmaxheight_greater2 << "\t" << voxcanopy_filled_greater2 << "\t" << voxcrown_total << "\t" << voxcrown_greater2 << "\t";

    // Now the more detailed metrics
    //  from the ground upwards
    for (int h = 0; h < 70; h++)
        output_transmittance << voxcrown[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << transmittance_nogaps_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << transmittance_incanopy_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << transmittance_full_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << voxcanopy_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << voxcanopy_na_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << voxcanopy_empty_z[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance << voxcanopy_filled_z[h] << "\t";

    // From the canopy top downwards
    for (int d = 0; d < 70; d++)
        output_transmittance << transmittance_nogaps_d[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << transmittance_full_d[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << voxna_d[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << voxempty_d[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance << voxfilled_d[d] << "\t";

    output_transmittance << nb_aggregates << "\t";

    for (int bin = 0; bin < 80; bin++)
    {
        if (nb_aggregates > 0)
            output_transmittance << PAIhist[bin] / float(nb_aggregates) << "\t";
        else
            output_transmittance << "NA\t";
    }

    output_transmittance << endl;

    // Same as before, but for the simulated ALS scan
    output_transmittance_ALS << ctx.time.iter << "\t" << sites_abc << "\t";

    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << transmittance_nogaps_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << transmittance_incanopy_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << transmittance_full_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << voxcanopy_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << voxcanopy_na_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << voxcanopy_empty_zALS[h] << "\t";
    for (int h = 0; h < 70; h++)
        output_transmittance_ALS << voxcanopy_filled_zALS[h] << "\t";

    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << transmittance_nogaps_dALS[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << transmittance_full_dALS[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << voxna_dALS[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << voxempty_dALS[d] << "\t";
    for (int d = 0; d < 70; d++)
        output_transmittance_ALS << voxfilled_dALS[d] << "\t";

    output_transmittance_ALS << nb_aggregates << "\t";

    for (int bin = 0; bin < 80; bin++)
    {
        if (nb_aggregates > 0)
            output_transmittance_ALS << PAIhist_ALS[bin] / float(nb_aggregates) << "\t";
        else
            output_transmittance_ALS << "NA\t";
    }

    output_transmittance_ALS << endl;
}
#endif

//! Close outputs
void CloseOutputs(Context &ctx)
{
    ctx.out.output_info.close();
    ctx.out.output_info.clear();

    for (int i = 0; i < 4; i++)
    {
        ctx.out.output_basic[i].close();
        ctx.out.output_basic[i].clear();
    }

    if (ctx.opt._OUTPUT_extended == 1)
    {
        for (int i = 0; i < 9; i++)
        {
            ctx.out.output_extended[i].close();
            ctx.out.output_extended[i].clear();
        }
        if (ctx.crown.extent_visual > 0)
        {
            for (int i = 0; i < 2; i++)
            {
                ctx.out.output_visual[i].close();
                ctx.out.output_visual[i].clear();
            }
        }
    }

#ifdef Output_ABC
    for (int i = 0; i < 11; i++)
    {
        ctx.out.output_abc[i].close();
        ctx.out.output_abc[i].clear();
    }
#endif
#ifdef WATER
    for (int i = 0; i < 40; i++)
    {
        ctx.out.output[i].close();
        ctx.out.output[i].clear();
    }
#endif
#ifdef TRACK_INDIVIDUALS
    for (int i = 0; i < 3; i++)
    {
        ctx.out.output_track[i].close();
        ctx.out.output_track[i].clear();
    }
#endif
}

void OutputCrownSliced(Context &ctx, int height, int site, int row_slice, vector<float> &output_statistics)
{
    int row_current = site / ctx.grid.cols;
    int col_current = site % ctx.grid.cols;
    if (row_current == row_slice && col_current >= ctx.crown.mincol_visual && col_current < ctx.crown.maxcol_visual)
    {
        ctx.out.output_visual[1] << ctx.time.iter << "\t" << row_current << "\t" << col_current << "\t" << height;
        for (int i = 0; i < output_statistics.size(); i++)
        {
            ctx.out.output_visual[1] << "\t" << output_statistics[i];
        } // we only output tree parts that fall in the current slice extent
        ctx.out.output_visual[1] << endl;
    }
};


