#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

//! Percentage num/den*100, or 0 when den is 0 (avoids printing "-nan" when a size class is empty)
//! @param num numerator
//! @param den denominator; 0 gives 0 instead of NaN
static double PercentOrZero(double num, double den)
{
    return den > 0.0 ? 100.0 * num / den : 0.0;
}

// ##############################################
//  Global function: calculation of the global averages every ctx.time.timestep
// ##############################################

//! @param ctx simulation context (all model state)
void Average(Context &ctx)
{

#ifdef TRACK_INDIVIDUALS
    TrackingData_andOutput(ctx);
#endif


    int site, spp;
    float sum1 = 0.0, sum10 = 0.0, sum30 = 0.0, ba = 0.0, npp = 0.0, gpp = 0.0, ba10 = 0.0, agb = 0.0, rday = 0.0, rnight = 0.0, rstem = 0.0, litterfall = 0.0;

    float inbcells = 1.0 / float(ctx.grid.sites);
    float inbhectares = inbcells * ctx.grid.NH * ctx.grid.NH * 10000.0;

#ifdef WATER
    ctx.soil.abund_phi_root = 0.0;
    ctx.soil.abund10_phi_root = 0.0;
    ctx.soil.agb_phi_root = 0.0;
#endif

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
        cout << ctx.time.iter << "\tTrees (1/ha): " << sum1 << " | " << sum10 << " | " << sum30 << " *** nbdead (%): " << PercentOrZero(ctx.diag.nbdead_n1 * inbhectares, sum1) << " | " << PercentOrZero(ctx.diag.nbdead_n10 * inbhectares, sum10) << " | " << PercentOrZero(ctx.diag.nbdead_n30 * inbhectares, sum30) << " *** AGB (t/ha): " << round(agb / 1000.0) << " GPP (MgC/ha/yr) " << gpp * ctx.time.iterperyear << " NPP " << npp * ctx.time.iterperyear << " litterfall (Mg/ha/yr) " << litterfall * ctx.time.iterperyear << " *** ctx.soil.Transpiration (mm): ";
#else

        cout << ctx.time.iter << "\tTrees (1/ha): " << sum1 << " | " << sum10 << " | " << sum30 << " *** nbdead (%): " << PercentOrZero(ctx.diag.nbdead_n1 * inbhectares, sum1) << " | " << PercentOrZero(ctx.diag.nbdead_n10 * inbhectares, sum10) << " | " << PercentOrZero(ctx.diag.nbdead_n30 * inbhectares, sum30) << " *** AGB (t/ha): " << round(agb / 1000.0) << " GPP (MgC/ha/yr) " << gpp * ctx.time.iterperyear << " NPP " << npp * ctx.time.iterperyear << " litterfall (Mg/ha/yr) " << litterfall * ctx.time.iterperyear << endl;
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
            float agb = 1000.0 * ctx.T[s].CalcAGB(ctx); // convert to g
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

    cout.flush();
}

// ##############################################
//  Global function: output of the field variables every ctx.time.timestep
// ##############################################
//! @param ctx simulation context (all model state)
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

        {
            // output of the dbh histograms (ctx.out.output[31])
            for (d = 1; d < ctx.grid.dbhmaxincm; d++)
                ctx.out.output[31] << d << "\t" << ctx.diag.nbdbh[d] << "\n";
            ctx.out.output[31] << "\n";
            // output of the mean LAI per height class (ctx.out.output[32])
            float norm = 1.0 / float(ctx.grid.sites * ctx.grid.LH * ctx.grid.LH);
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
//! @param ctx simulation context (all model state)
//! @param output output stream (file of the snapshot)
//! @param header write the header line first if true
//! @param dbh_limit minimum DBH of the trees written (m)
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

        output << "\tAGB\ts_name" << endl;
    }

    // reset the canopy to make sure that it is well-constructed
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        for (int sbsite = 0; sbsite < ctx.grid.sites + 2 * ctx.grid.SBORD; sbsite++)
            ctx.field.LAI3D[h][sbsite] = 0.0;

    for (int site = 0; site < ctx.grid.sites; site++)
    { // Each tree contribues to ctx.field.LAI3D
        ctx.T[site].CalcLAI(ctx);
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
//! Compute the spike-free canopy height of each site into chm_spikefree (used for CHM.txt; CHM_SPIKEFREE option)
//! @param ctx simulation context (all model state)
//! @param chm_spikefree [out] spike-free canopy height per site (m)
void MakeCHMspikefree(Context &ctx, vector<int> &chm_spikefree)
{
    chm_spikefree.clear();
    chm_spikefree.reserve(ctx.grid.sites);
    for (int s = 0; s < ctx.grid.sites; s++)
        chm_spikefree.push_back(0);

}
#endif

// ################################
// ### Output for visualization ###
// ################################
//! @param ctx simulation context (all model state)
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
//! @param ctx simulation context (all model state)
//! @param output_CHM output stream of CHM.txt
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
//! @param ctx simulation context (all model state)
//! @param output_transmLAI3D output stream of LAI.txt
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

//! @param ctx simulation context (all model state)
//! @param beams [in,out] number of laser beams per site
//! @param beams_returns [in,out] laser returns per voxel
//! @param mean_beam mean number of beams per site
//! @param sd_beam standard deviation of the number of beams per site
//! @param klaser extinction coefficient of the laser in foliage
//! @param transmittance_laser transmittance of the laser through empty voxels
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
    snprintf(system_identifier, sizeof(system_identifier),"ALS simulator");
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
//! @param ctx simulation context (all model state)
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



//! Close outputs
namespace
{
    void CloseStream(std::fstream &stream)
    {
        if (stream.is_open())
            stream.close();
        stream.clear();
    }

    template <std::size_t N>
    void CloseStreams(std::fstream (&streams)[N])
    {
        for (auto &stream : streams)
            CloseStream(stream);
    }
} // namespace

//! Close (and thereby flush) every output stream, whether or not it was opened.
//! Previously only a subset was closed and the rest relied on exit() running global
//! destructors; with a local Context those buffers could be lost.
//! @param ctx simulation context (all model state)
void CloseOutputs(Context &ctx)
{
    CloseStream(ctx.out.output_info);
    CloseStreams(ctx.out.output_basic);
    CloseStreams(ctx.out.output_extended);
    CloseStreams(ctx.out.output_visual);
    CloseStream(ctx.out.output_pointcloud);
#ifdef MIP_Lichstein
    CloseStream(ctx.out.output_MIP_eco);
    CloseStream(ctx.out.output_MIP_ind);
#endif
#ifdef WATER
    CloseStreams(ctx.out.output);
#endif
#ifdef TRACK_INDIVIDUALS
    CloseStreams(ctx.out.output_track);
#endif
}
