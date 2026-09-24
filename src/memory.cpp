#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

//! Global function: initialise intraspecific variables
//! Global function: initialise output streams
void InitialiseOutputStreams(Context &ctx)
{
    char nnn[200];
    if (!mpi_rank)
    {
        sprintf(nnn, "%s_%i_sumstats.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_basic[0].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_initial_pattern.txt", ctx.buffers.buf, easympi_rank); // previously "state" output, but not used anymore, overwritten for initial pattern
        ctx.out.output_basic[1].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_final_pattern.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_basic[2].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_final_SWC3D.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_basic[3].open(nnn, ios::out);

        // write headers for files
        ctx.out.output_basic[0] << "iter\tsum1\tsum10\tsum30\tba\tba10\tagb\tgpp\tnpp\trday\trnight\trstem\tlitterfall" << endl;
        // headers for initial and final patterns are written automatically

        sprintf(nnn, "%s_%i_sumstats_species.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_extended[0].open(nnn, ios::out);
        ctx.out.output_extended[0] << "iter\tspecies\tsum1\tsum10\tsum30\tba\tba10\tagb\tgpp\tnpp\trday\trnight\trstem\tlitterfall" << endl;

#ifdef MIP_Lichstein
        sprintf(nnn, "%s_%i_MIP_eco.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_MIP_eco.open(nnn, ios::out);
        sprintf(nnn, "%s_%i_MIP_ind.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_MIP_ind.open(nnn, ios::out);

        ctx.out.output_MIP_eco << "YEAR\tMONTH\tDAY\tGPP\tNPP\tET\tLAI\tLFLIT\tSW1\tSW2\tSW3\tSW4" << endl;
        ctx.out.output_MIP_ind << "YEAR\tSP\tID\tNLIVE\tNDEAD\tDBH\tHT\tTB\tAGB\tWD\tSLA\tNMASS\tPMASS\tDBHMAX\tTLP\tLA" << endl;
#endif

        if (ctx.opt._OUTPUT_extended)
        {
            sprintf(nnn, "%s_%i_sumstats_species.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[0].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_ppfd0.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[1].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_death.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[2].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_death_snapshots.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[3].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_deathrate.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[4].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_sdd.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[5].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_vertd.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[6].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_LAI.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[7].open(nnn, ios::out);
            sprintf(nnn, "%s_%i_CHM.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_extended[8].open(nnn, ios::out);

            // write headers
            ctx.out.output_extended[0] << "iter\tspecies\tsum1\tsum10\tsum30\tba\tba10\tagb\tgpp\tnpp\trday\trnight\trstem\tlitterfall" << endl;
            ctx.out.output_extended[1] << "iter\tvariable\tvalue\tresidual" << endl;
            if (ctx.opt._BASICTREEFALL)
                ctx.out.output_extended[2] << "iter\tnbdead_n1\tnbdead_n10\tnbTreefall1\tnbTreefall10" << endl;
            else
                ctx.out.output_extended[2] << "iter\tnbdead_n1\tnbdead_n10" << endl;
            ctx.out.output_extended[3] << "iter\tspecies\tage\tdbh\theight" << endl;
            ctx.out.output_extended[4] << "iter\twsg\tdbh\tbasal\tdr" << endl;
            ctx.out.output_extended[5] << "iter\td\tfreq" << endl;
            ctx.out.output_extended[6] << "iter\th\tfreq" << endl;

            if (ctx.crown.extent_visual > 0)
            {
                sprintf(nnn, "%s_%i_visual_field.txt", ctx.buffers.buf, easympi_rank);
                ctx.out.output_visual[0].open(nnn, ios::out);
#ifdef CHM_SPIKEFREE
                ctx.out.output_visual[0] << "iter" << "\t" << "row" << "\t" << "col" << "\t" << "height" << "\t" << "height_spikefree" << "\t" << "LAI" << endl; // header
#else
                ctx.out.output_visual[0] << "iter" << "\t" << "row" << "\t" << "col" << "\t" << "height" << "\t" << "LAI" << endl; // header
#endif

                sprintf(nnn, "%s_%i_visual_slice.txt", ctx.buffers.buf, easympi_rank);
                ctx.out.output_visual[1].open(nnn, ios::out);
                ctx.out.output_visual[1] << "iter" << "\t" << "row" << "\t" << "col" << "\t" << "height" << "\t" << "sp_lab" << "\t" << "ratio_height_Ct" << "\t" << "ratio_NPP_GPP" << endl; // header
            }
        }

        // v.3.1.6 output for point cloud
        if (ctx.opt._OUTPUT_pointcloud)
        {
            sprintf(nnn, "%s_%i.las", ctx.buffers.buf, easympi_rank);
            ctx.out.output_pointcloud.open(nnn, ios::out | ios::binary);
            ctx.out.output_pointcloud.imbue(locale::classic()); // justification here: https://stackoverflow.com/questions/14750496/sending-integer-to-fstream-as-little-endian; locale regulates how streams print and read values (i.e. commas vs. points for decimals, etc.); setting it to classic to ensure portability, but not entirely sure how important this is in practice for binary files
        }

#ifdef Output_ABC
        sprintf(nnn, "%s_%i_abc_traitconservation.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[0].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_ground.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[1].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_chm.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[2].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_chmALS.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[3].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_transmittance.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[4].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_transmittanceALS.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[5].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_species.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[6].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_species10.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[7].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_traits.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[8].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_traits10.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[9].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_abc_biomass.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output_abc[10].open(nnn, ios::out);
#endif

#ifdef WATER
        // sprintf(nnn,"%s_%i_water_balance.txt",ctx.buffers.buf, easympi_rank);
        // output_water[0].open(nnn, ios::out);

        sprintf(nnn, "%s_%i_water_balance.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[11].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site1.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[12].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site2.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[13].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site3.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[14].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site4.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[15].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site5.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[16].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_site6.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[17].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_state_begin.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[28].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_state_mid.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[29].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_state_end.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[30].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWC_begin.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[1].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWC_mid.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[2].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWC_end.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[3].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWP_begin.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[4].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWP_mid.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[5].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_SWP_end.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[6].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_litterfall.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[7].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_waterfluxes_begin.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[18].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_waterfluxes_mid.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[19].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_waterfluxes_end.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[20].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_LAIdynamics.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[21].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_LAIyoung.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[22].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_LAImature.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[23].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_LAIold.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[24].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_phi_root.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[31].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_LAIprofile.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[32].open(nnn, ios::out);
        sprintf(nnn, "%s_%i_soilproperties.txt", ctx.buffers.buf, easympi_rank);
        ctx.out.output[33].open(nnn, ios::out);

        // write headers
        ctx.out.output[21] << "iter\tLAI" << endl;
        ctx.out.output[22] << "iter";
        ctx.out.output[23] << "iter";
        ctx.out.output[24] << "iter";
        for (int l = 0; l < ctx.grid.HEIGHT + 1; l++)
        {
            ctx.out.output[22] << "\t" << "h_" << l;
            ctx.out.output[23] << "\t" << "h_" << l;
            ctx.out.output[24] << "\t" << "h_" << l;
        }
        ctx.out.output[22] << endl;
        ctx.out.output[23] << endl;
        ctx.out.output[24] << endl;
        ctx.out.output[31] << "iter\tabundance_weighted_phi_root\tabundance_weighted_phi_root_n10\tagb_weighted_phi_root" << endl;
        ctx.out.output[32] << "iter\theight\tLAI" << endl;
        ctx.out.output[1] << "layer";
        ctx.out.output[2] << "layer";
        ctx.out.output[3] << "layer";
        ctx.out.output[4] << "layer";
        ctx.out.output[5] << "layer";
        ctx.out.output[6] << "layer";
        ctx.out.output[18] << "variable";
        ctx.out.output[19] << "variable";
        ctx.out.output[20] << "variable";
        for (int d = 0; d < ctx.grid.nbdcells; d++)
        {
            ctx.out.output[1] << "\t" << "dcell_" << d;
            ctx.out.output[2] << "\t" << "dcell_" << d;
            ctx.out.output[3] << "\t" << "dcell_" << d;
            ctx.out.output[4] << "\t" << "dcell_" << d;
            ctx.out.output[5] << "\t" << "dcell_" << d;
            ctx.out.output[6] << "\t" << "dcell_" << d;
            ctx.out.output[18] << "\t" << "dcell_" << d;
            ctx.out.output[19] << "\t" << "dcell_" << d;
            ctx.out.output[20] << "\t" << "dcell_" << d;
        }
        ctx.out.output[1] << endl;
        ctx.out.output[2] << endl;
        ctx.out.output[3] << endl;
        ctx.out.output[4] << endl;
        ctx.out.output[5] << endl;
        ctx.out.output[6] << endl;
        ctx.out.output[18] << endl;
        ctx.out.output[19] << endl;
        ctx.out.output[20] << endl;
        ctx.out.output[11] << "iter\tprecipitation\tinterception\tthroughfall\trunoff\tleak\tevaporation";
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            ctx.out.output[11] << "\t" << "transpitation_" << l;
        ctx.out.output[11] << "\t" << "transpiration1016";
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            ctx.out.output[11] << "\t" << "SWC_" << l;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            ctx.out.output[11] << "\t" << "SWP_" << l;
        ctx.out.output[11] << endl;
        string sites_headers = "iter\tsite\tsp_lab\theight\tdbh\tlitter\tage\tLA\tyougLA\tmatureLA\toldLA\tCR\tCD\tGPP\tNPP\tRstem\tRnight\tLAI3D_canopy\tLAI3D_base_crown\troot_depth\tphi_root\tWSFs\tWSFns\ttranspiration\tLAImax\tLAmax";
        ctx.out.output[12] << sites_headers;
        ctx.out.output[13] << sites_headers;
        ctx.out.output[14] << sites_headers;
        ctx.out.output[15] << sites_headers;
        ctx.out.output[16] << sites_headers;
        ctx.out.output[17] << sites_headers;
        ctx.out.output[28] << sites_headers;
        ctx.out.output[29] << sites_headers;
        ctx.out.output[30] << sites_headers;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            ctx.out.output[12] << "\t" << "root_biomass_" << l;
            ctx.out.output[13] << "\t" << "root_biomass_" << l;
            ctx.out.output[14] << "\t" << "root_biomass_" << l;
            ctx.out.output[15] << "\t" << "root_biomass_" << l;
            ctx.out.output[16] << "\t" << "root_biomass_" << l;
            ctx.out.output[17] << "\t" << "root_biomass_" << l;
            ctx.out.output[28] << "\t" << "root_biomass_" << l;
            ctx.out.output[29] << "\t" << "root_biomass_" << l;
            ctx.out.output[30] << "\t" << "root_biomass_" << l;
        }
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            ctx.out.output[12] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[13] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[14] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[15] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[16] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[17] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[28] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[29] << "\t" << "soil_layer_weight_" << l;
            ctx.out.output[30] << "\t" << "soil_layer_weight_" << l;
        }
        ctx.out.output[12] << endl;
        ctx.out.output[13] << endl;
        ctx.out.output[14] << endl;
        ctx.out.output[15] << endl;
        ctx.out.output[16] << endl;
        ctx.out.output[17] << endl;
        ctx.out.output[28] << endl;
        ctx.out.output[29] << endl;
        ctx.out.output[30] << endl;

#endif

#ifdef WATER

        if (ctx.opt._WATER_RETENTION_CURVE == 1)
        {
            ctx.out.output[33] << "layer" << "\t" << "depth" << "\t" << "sat" << "\t" << "max" << "\t" << "fc" << "\t" << "res" << "\t" << "min" << "\t" << "Ksat" << "\t" << "a_vgm" << "\t" << "m_vgm" << endl;
            for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            {
                ctx.out.output[33] << l << "\t" << ctx.soil.layer_depth[l] << "\t" << ctx.soil.Sat_SWC[l] << "\t" << ctx.soil.Max_SWC[l] << "\t" << ctx.soil.FC_SWC[l] << "\t" << ctx.soil.Res_SWC[l] << "\t" << ctx.soil.Min_SWC[l] << "\t" << ctx.soil.Ksat[l] << "\t" << ctx.soil.a_vgm[l] << "\t" << ctx.soil.m_vgm[l] << endl;
            }
        }
        else if (ctx.opt._WATER_RETENTION_CURVE == 0)
        {

            ctx.out.output[33] << "layer" << "\t" << "depth" << "\t" << "sat" << "\t" << "max" << "\t" << "fc" << "\t" << "res" << "\t" << "min" << "\t" << "Ksat" << "\t" << "phi_e" << "\t" << "b" << endl;
            for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            {
                ctx.out.output[33] << l << "\t" << ctx.soil.layer_depth[l] << "\t" << ctx.soil.Sat_SWC[l] << "\t" << ctx.soil.Max_SWC[l] << "\t" << ctx.soil.FC_SWC[l] << "\t" << ctx.soil.Res_SWC[l] << "\t" << ctx.soil.Min_SWC[l] << "\t" << ctx.soil.Ksat[l] << "\t" << ctx.soil.phi_e[l] << "\t" << ctx.soil.b[l] << endl;
            }
        }

#endif

#ifdef TRACK_INDIVIDUALS
        if (ctx.opt._OUTPUT_extended)
        {
            // these are the "cases" of trees that are followed
            sprintf(nnn, "%s_%i_trees_fortracking.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_track[0].open(nnn, ios::out);
            ctx.out.output_track[0] << "site" << "\t" << "timeofyear_born" << "\t" << "col" << "\t" << "row" << "\t" << "species" << "\t" << "dbh" << "\t" << "cr" << "\t" << "height" << "\t" << "agb" << "\t" << "multiplier_cr" << "\t" << "multiplier_height" << "\t" << "wsg" << "\t" << "Nmass" << "\t" << "Pmass" << "\t" << "LMA" << "\t" << "deviation_wsg" << "\t" << "multiplier_Nmass" << "\t" << "multiplier_Pmass" << "\t" << "multiplier_lma" << "\t" << "Vcmax" << "\t" << "Jmax" << "\t" << "Rdark" << "\t" << "LAImax" << "\t" << "leaflifespan" << endl;

            // these are the "activities" that are recorded for each tree
            sprintf(nnn, "%s_%i_trees_tracked.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_track[1].open(nnn, ios::out);

            ctx.out.output_track[1] << "site" << "\t" << "timeofyear_born" << "\t" << "Iter" << "\t" << "age" << "\t" << "seeds" << "\t" << "seedstotal" << "\t" << "carbstarv" << "\t" << "carbstarvtotal" << "\t" << "dbh" << "\t" << "dbhgrowth" << "\t" << "height" << "\t" << "heightgrowth" << "\t" << "cr" << "\t" << "crgrowth" << "\t" << "agb" << "\t" << "agbgrowth" << "\t" << "GPP" << "\t" << "GPPsq" << "\t" << "NPP" << "\t" << "NPPsq" << "\t" << "Rday" << "\t" << "Rnight" << "\t" << "Rstem" << "\t" << "LAIabove_avg" << "\t" << "LAIabove_effavg" << "\t" << "carbstore_avg" << endl;

            // this is to get information on the dead trees > 10cm
            sprintf(nnn, "%s_%i_trees_aftertracking.txt", ctx.buffers.buf, easympi_rank);
            ctx.out.output_track[2].open(nnn, ios::out);

            ctx.out.output_track[2] << "site" << "\t" << "timeofyear_born" << "\t" << "Iter" << "\t" << "age" << "\t" << "seedstotal" << "\t" << "carbstarvtotal" << "\t" << "dbh" << "\t" << "height" << "\t" << "cr" << "\t" << "agb" << "\t" << "GPP" << "\t" << "NPP" << "\t" << "LAIabove_avg" << "\t" << "LAIabove_effavg" << "\t" << "GPPsquared" << "\t" << "NPPsquared" << "\t" << "LAIabovesquared_avg" << "\t" << "LAIabovesquared_effavg" << endl;
        }
#endif
    }
}

// ######################################
//  Global function: initialisation with bare ground conditions
// ######################################
void Initialise(Context &ctx)
{
    //** Initialization of the simulation parameters **
    //*************************************************

    ctx.time.iter = -1; // changed in v.3.0.1, previously undefined; new function GetTimeofyear accepts also negative iterations

    ctx.diag.nblivetrees = 0;

    ReadInputSpecies(ctx);
#ifdef FULL_CLIMATE
    ReadInputClimate(ctx);
    ReadInputDailyvar(ctx);
#else
    ReadInputDailyvar(ctx);
    ReadInputClimate(ctx);
#endif

#ifdef WATER
    ReadInputSoil(ctx);
#endif

    //** Initialization of trees **
    //*****************************
    ctx.T.reserve(ctx.grid.sites);
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        Tree T_site(ctx);
        ctx.T.push_back(T_site);
#ifdef WATER
        // FF: not sure this check is necessary anymore (it's a check for memory allocation problems, I presume?), but I kept it just in case
        if (&T_site.t_soil_layer_weight[0] == &T_site.t_root_biomass[4])
        {
            cout << "Warning mem_alloc root biomass and soil ctx.diag.layer weight at site " << site << ": " << endl;
            cout << "t_soil_layer_weight adresses: " << &T_site.t_soil_layer_weight[0] << "\t" << &T_site.t_soil_layer_weight[1] << "\t" << &T_site.t_soil_layer_weight[2] << "\t" << &T_site.t_soil_layer_weight[3] << "\t" << &T_site.t_soil_layer_weight[4] << endl;
            cout << "t_root biomass adresses: " << &T_site.t_root_biomass[0] << "\t" << &T_site.t_root_biomass[1] << "\t" << &T_site.t_root_biomass[2] << "\t" << &T_site.t_root_biomass[3] << "\t" << &T_site.t_root_biomass[4] << endl;
        }
#endif
    }
    InitialiseIntraspecific(ctx);
#ifdef CACHE_LUT
    if (!LoadLookUpTablesFromCache(ctx))
    {
        std::cout << "LUT cache missing → computing LUTs..." << std::endl;
        InitialiseLookUpTables(ctx);
        SaveLookUpTablesToCache(ctx);
    }
    else
    {
        std::cout << "Loaded LUTs from cache." << std::endl;
    }
#else  // CACHE_LUT
    InitialiseLookUpTables(ctx);
#endif // CACHE_LUT
}

#ifdef Output_ABC
// ######################################
//  Global ABC function: initialise ABC conditions
// ######################################
//! - determine area to be taken into consideration for ABC and margin not to be taken into consideration (artefacts along border)
void InitialiseABC(Context &ctx)
{
    margin = 0;
    row_start = margin;
    col_start = margin;
    row_end = ctx.grid.rows - margin;
    col_end = ctx.grid.cols - margin;
    sites_abc = (row_end - row_start) * (col_end - col_start);
    isites_abc = 1.0 / float(sites_abc);
    cout << "row start: " << row_start << " | row end: " << row_end << " | sites_abc: " << sites_abc << endl;
    // counter for recursive function (gap patches)
    nbvisited = 0;
    // distributions for simulated and empirical CHM
    // precomputed, not directly deduced from empirical fields, as a separate algorithm is used (from LAStools)
    if (NULL == (chm_field_previous = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (chm_field_current = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (chm_field_previous_ALS = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (chm_field_current_ALS = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (chm_field_changes = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (chm_field_changes_ALS = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    for (int s = 0; s < ctx.grid.sites; s++)
    {
        chm_field_previous[s] = 0;
        chm_field_current[s] = 0;
        chm_field_previous_ALS[s] = 0;
        chm_field_current_ALS[s] = 0;
        chm_field_changes[s] = 0;
        chm_field_changes_ALS[s] = 0;
    }
    // field for simulated transmittance
    if (NULL == (transmittance_simulatedALS = new float *[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (transmittance_direct = new float *[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (transmittance_simulatedALS_sampling = new int *[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        if (NULL == (transmittance_simulatedALS[h] = new float[ctx.grid.sites]))
            cerr << "!!! Mem_Alloc\n";
        if (NULL == (transmittance_direct[h] = new float[ctx.grid.sites]))
            cerr << "!!! Mem_Alloc\n";
        if (NULL == (transmittance_simulatedALS_sampling[h] = new int[ctx.grid.sites]))
            cerr << "!!! Mem_Alloc\n";
    }
    // set it to default 1 (full visibility)
    for (int row = 0; row < ctx.grid.rows; row++)
    {
        for (int col = 0; col < ctx.grid.cols; col++)
        {
            for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            {
                transmittance_simulatedALS[h][col + row * ctx.grid.cols] = 1.0;
                transmittance_direct[h][col + row * ctx.grid.cols] = 1.0;
                transmittance_simulatedALS_sampling[h][col + row * ctx.grid.cols] = 0;
            }
        }
    }
    // initialise output streams
    OutputABCWriteHeaders(ctx.out.output[11], ctx.out.output[12], ctx.out.output[13], ctx.out.output[14], ctx.out.output[15], ctx.out.output[16], ctx.out.output[17], ctx.out.output[18], ctx.out.output[19], ctx.out.output[23], ctx.out.output[24], ctx.out.output[25], ctx.out.output[26], ctx.out.output[27]);
}
#endif

// ######################################
//  Global function: point cloud generation
// ######################################
//  added v.3.1.6

void ReadInputPointcloud(Context &ctx)
{
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_pointcloud << endl;

    fstream InPointcloud(ctx.fileio.inputfile_pointcloud, ios::in);

    if (InPointcloud)
    {
        string parameter_names[5] = {"ctx.pc.mean_beam_pc", "ctx.pc.sd_beam_pc", "ctx.pc.klaser_pc", "ctx.pc.transmittance_laser", "ctx.pc.iter_pointcloud_generation"};
        int nb_parameters = 5;
        vector<string> parameter_values(nb_parameters, "");

        cout << endl
             << "Reading in file: " << ctx.fileio.inputfile << endl;
        InPointcloud.getline(ctx.fileio.buffer, 256, '\n');
        string parameter_name, parameter_value;

        while (InPointcloud >> parameter_name >> parameter_value)
        {
            InPointcloud.getline(ctx.fileio.buffer, 256, '\n');
            for (int i = 0; i < nb_parameters; i++)
            {
                if (parameter_name == parameter_names[i])
                    parameter_values[i] = parameter_value;
            }
        }
        // now we assign values
        for (int i = 0; i < nb_parameters; i++)
        {
            AssignValuePointcloud(ctx, parameter_names[i], parameter_values[i]);
        }
    }
}

// ######################################
//  Global function: initialisation from inventories
// ######################################

// completely rewritten in v.3.1
void ReadInputInventory(Context &ctx)
{

    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_inventory << endl;

    fstream InInventory(ctx.fileio.inputfile_inventory, ios::in);

    if (InInventory)
    {
        // possible parameters to initialise
        //        vector<string> parameter_names{"col","row","s_name","CrownDisplacement","Pmass","Nmass","LMA","wsg","Rdark","Vcmax","Jmax","leaflifespan","lambda_young","lambda_mature","lambda_old","dbhmature","dbhmax","hmax","ah","Ct","LAImax","fraction_filled","mult_height","mult_CR","mult_CD","mult_P","mult_N","mult_LMA","mult_dbhmax","dev_wsg","age","dbh","sapwood_area","height","CD","CR","GPP","NPP","Rday","Rnight","Rstem","LAmax","LA","youngLA","matureLA","oldLA","LAI","litter","carbon_storage","carbon_biometry","multiplier_seed","hurt","NPPneg"};
        //        int nb_parameters = int(parameter_names.size()); // only works from C++11 onwards

#ifdef WATER
        string parameter_names_hardcoded[67] = {"col", "row", "s_name", "CrownDisplacement",
                                                "Pmass", "Nmass", "LMA", "leafarea", "tlp", "wsg",
                                                "Rdark", "Vcmax", "Jmax", "leaflifespan",
                                                "lambda_young", "lambda_mature", "lambda_old",
                                                "dbhmature", "dbhmax", "hmax", "ah", "Ct",
                                                "LAImax", "fraction_filled",
                                                "mult_height", "mult_CR", "mult_CD", "mult_P", "mult_N", "mult_LMA",
                                                "mult_dbhmax", "mult_leafarea", "mult_tlp",
                                                "dev_wsg", "age", "dbh", "sapwood_area", "height", "CD", "CR",
                                                "GPP", "NPP", "Rday", "Rnight", "Rstem",
                                                "LAmax", "LA", "youngLA", "matureLA", "oldLA", "LAI", "litter",
                                                "carbon_storage", "carbon_biometry", "multiplier_seed", "hurt", "NPPneg",
                                                "root_depth", "phi_root", "WSF", "WSF_A", "transpiration", "g1_0", "g1",
                                                "Ndays_dry", "Ndays_wet", "pheno_factor"};

        int nb_parameters = 67;

#else

        string parameter_names_hardcoded[53] = {"col", "row", "s_name", "CrownDisplacement", "Pmass", "Nmass", "LMA", "wsg", "Rdark", "Vcmax", "Jmax", "leaflifespan", "lambda_young", "lambda_mature", "lambda_old", "dbhmature", "dbhmax", "hmax", "ah", "Ct", "LAImax", "fraction_filled", "mult_height", "mult_CR", "mult_CD", "mult_P", "mult_N", "mult_LMA", "mult_dbhmax", "dev_wsg", "age", "dbh", "sapwood_area", "height", "CD", "CR", "GPP", "NPP", "Rday", "Rnight", "Rstem", "LAmax", "LA", "youngLA", "matureLA", "oldLA", "LAI", "litter", "carbon_storage", "carbon_biometry", "multiplier_seed", "hurt", "NPPneg"};

        int nb_parameters = 53;

#endif

        // conversion to vector for future use (e.g. water module needs flexible number of input layers)
        vector<string> parameter_names(parameter_names_hardcoded, parameter_names_hardcoded + nb_parameters);

#ifdef WATER

        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            char par_name[30];
            sprintf(par_name, "root_biomass%d", l);
            parameter_names.push_back(par_name);
        }
        nb_parameters += ctx.soil.nblayers_soil;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            char par_name[30];
            sprintf(par_name, "soil_layer_weight%d", l);
            parameter_names.push_back(par_name);
        }
        nb_parameters += ctx.soil.nblayers_soil;

#endif

        // define position of col and row in string vector
        int pcol = 0, prow = 0;
        while (pcol < nb_parameters && parameter_names[pcol] != "col")
            pcol++;
        while (prow < nb_parameters && parameter_names[prow] != "row")
            prow++;

        // first get parameter names
        string line;
        getline(InInventory, line);
        istringstream firstlinestream(line);

        string parameter_name;
        vector<int> parameter_index;
        int nb_parameters_fromfile = 0;

        int flag_col = 1, flag_row = 1, flag_dbh = 1, flag_species = 1;
        while (firstlinestream >> parameter_name)
        {
            int index = -1;
            for (int i = 0; i < nb_parameters; i++)
            {
                if (parameter_name == parameter_names[i])
                {
                    index = i;
                    nb_parameters_fromfile++;
                }
            }
            if (index == -1)
            {
                cout << "Ignoring unknown parameter: " << parameter_name << ". Parameters must be one of:";
                for (int i = 0; i < nb_parameters; i++)
                {
                    cout << "\t" << parameter_names[i];
                }
                cout << endl;
            }
            parameter_index.push_back(index);

            if (parameter_name == "col")
                flag_col = 0;
            if (parameter_name == "row")
                flag_row = 0;
            if (parameter_name == "dbh")
                flag_dbh = 0;
            if (parameter_name == "s_name")
                flag_species = 0;
        }

        // Warnings
        if (flag_dbh == 1)
        {
            cout << "WARNING! No diameter column provided, no initialization from file will be carried out." << endl;
        }
        else
        {
            if ((flag_col == 1) | (flag_row == 1))
                cout << "WARNING! At least one coordinate column (col/row) is missing, random coordinates chosen" << endl;
            if (flag_species == 1)
                cout << "WARNING! No species column provided, random species chosen" << endl;

            int nb_parameterlines = 0;
            int nb_individuals = 0;
            int nb_random = 0;
            int nb_moved = 0;
            int nb_speciesrandom = 0;

            // define a radius to look for alternative tree positions in case site is already initialized
            int radius_max = 5;
            int area_max = int(radius_max * radius_max * PI);

            // define an array of randomly permuted ctx.grid.sites so that random number generator does not have to be called too often
            vector<int> sites_shuffled;
            sites_shuffled.reserve(ctx.grid.sites);
            for (int s = 0; s < ctx.grid.sites; s++)
            {
                sites_shuffled.push_back(s);
            }

            shuffle(sites_shuffled.begin(), sites_shuffled.end(), default_random_engine());

            int sites_shuffled_index = 0;

            // we go through all lines in the input file */
            while (getline(InInventory, line))
            {
                istringstream linestream(line);

                // first we initialise the potential parameter arrays with empty parameter values
                vector<string> parameter_values(nb_parameters, "");

                // we update the values from the file
                int s = 0;
                string parameter_value;
                while (linestream >> parameter_value)
                {
                    int index = parameter_index[s];
                    if (index >= 0)
                    {
                        parameter_values[index] = parameter_value;
                    }
                    s++;
                }

                // additional checks on location of tree
                int col = -1, row = -1;
                bool quiet = 1;
                if (flag_col == 0 && flag_row == 0)
                {
                    SetParameter(parameter_names[pcol], parameter_values[pcol], col, 0, ctx.grid.cols - 1, -1, quiet);
                    SetParameter(parameter_names[prow], parameter_values[prow], row, 0, ctx.grid.rows - 1, -1, quiet);
                }

                if (col >= 0 && row >= 0)
                {
                    int site = col + row * ctx.grid.cols;
                    int i = 0;
                    while (ctx.T[site].t_age != 0.0 && i < area_max)
                    {
                        i++;
                        int site_relative = ctx.lookup.LookUp_Crown_site[i];
                        int row_new = row + site_relative / 51 - 25;
                        int col_new = col + site_relative % 51 - 25;
                        if (row_new >= 0 && row_new < ctx.grid.rows && col_new >= 0 && col_new < ctx.grid.cols)
                        {
                            site = col_new + row_new * ctx.grid.cols;
                        }
                    }

                    if (ctx.T[site].t_age == 0.0)
                    {
                        int success = ctx.T[site].BirthFromInventory(ctx, site, parameter_names, parameter_values, nb_speciesrandom);
                        if (success == 1)
                        {
                            nb_individuals++;
                            if (i > 0)
                                nb_moved++;
                        }
                    }
                }
                else
                {
                    // find a random free site
                    while (sites_shuffled_index < ctx.grid.sites && ctx.T[sites_shuffled[sites_shuffled_index]].t_age != 0.0)
                        sites_shuffled_index++;

                    if (sites_shuffled_index < ctx.grid.sites)
                    {
                        // as long as the search has stopped and the index has not run outside the range, the tree can be initialized
                        int site = sites_shuffled[sites_shuffled_index];
                        int success = ctx.T[site].BirthFromInventory(ctx, site, parameter_names, parameter_values, nb_speciesrandom);
                        if (success == 1)
                        {
                            nb_individuals++;
                            nb_random++;
                        }
                    }
                }

                nb_parameterlines++;
            }

            if (nb_parameterlines > 0)
            {
                cout << "Successfully initialised " << nb_individuals << " out of " << nb_parameterlines << " trees from file." << endl;
                cout << "Coordinates: " << nb_moved << " trees were moved due to overlapping coordinates, and " << nb_random << " were placed randomly on the grid as coordinates were incomplete or could not be read in." << endl;
                cout << "Species: " << nb_speciesrandom << " trees were assigned a random species." << endl;
                if (sites_shuffled_index == ctx.grid.sites)
                    cout << "WARNING: shuffle index: " << sites_shuffled_index << " is equal site number. This means that random placement of trees has stopped early as the algorithm has run out of grid cells to place trees. This may indicate that there were more trees than grid cells." << endl;
            }
            else
            {
                cout << "WARNING! Inventory file was empty. No trees were initialised." << endl;
            }
        }
    }
    else
    {
        cout << "ERROR with the inventory file" << endl;
    }

    InInventory.close();

    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_SWC << endl;

    fstream InSWC(ctx.fileio.inputfile_SWC, ios::in);

    if (InSWC)
    {

        string line;
        int l = 0;
        int dcell = 0;
        while (getline(InSWC, line))
        {
            istringstream linestream(line);
            linestream >> dcell;
            float SWC_value;
            while (linestream >> SWC_value)
            {
                ctx.soil.SWC3D[l][dcell] = SWC_value;
                l++;
            }
            if (l != ctx.soil.nblayers_soil)
                cout << "Error in SWC input file (number of layers for dcell " << dcell << endl;
            l = 0;
            dcell++;
        }

        if (dcell != ctx.grid.nbdcells)
            cout << "Error in SWC input file (number of dcells) " << endl;
    }
    else
    {
        cout << "ERROR with the input SWC file" << endl;
    }

    InSWC.close();

    // compute ctx.soil.soil_phi3D for RecruitTree function
    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            float theta_w = (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) / (ctx.soil.Max_SWC[l] - ctx.soil.Min_SWC[l]);
            if (ctx.opt._WATER_RETENTION_CURVE == 1)
            {
                if (theta_w == 0)
                {
                    theta_w = 0.001; // SS addition for limit value
                    cout << "Warning theta_w = 0 " << endl;
                }
                ctx.soil.soil_phi3D[l][d] = ctx.soil.a_vgm[l] * pow((pow(theta_w, -ctx.soil.b_vgm[l]) - 1), ctx.soil.c_vgm[l]); // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                float inter = 1 - pow((1 - pow(theta_w, ctx.soil.b_vgm[l])), ctx.soil.m_vgm[l]);
                ctx.soil.Ks[l][d] = ctx.soil.Ksat[l] * pow(theta_w, 0.5) * inter * inter;                           // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                if (isnan(ctx.soil.soil_phi3D[l][d]) || isnan(ctx.soil.Ks[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0) //|| ctx.soil.KsPhi[l][d]==0.0 || ctx.soil.Ks[l][d]==0.0 || ctx.soil.soil_phi3D[l][d]==0.0)
                    cout << "In bucket model, ctx.diag.layer " << l << " dcell " << d << " theta_w=" << theta_w << " ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) << " ctx.soil.soil_phi3D[l][d]=" << ctx.soil.soil_phi3D[l][d] << " ctx.soil.Ksat=" << ctx.soil.Ksat[l] << " ctx.soil.Ks[l][d]=" << ctx.soil.Ks[l][d] << endl;
            }
            else if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {
                ctx.soil.soil_phi3D[l][d] = ctx.soil.phi_e[l] * pow(theta_w, -ctx.soil.b[l]);                                                      // this is the soil water characteristic of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                ctx.soil.Ks[l][d] = ctx.soil.Ksat[l] * pow(theta_w, 2.5 + 2 * ctx.soil.b[l]);                                                      // this is the hydraulic conductivity curve of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                ctx.soil.KsPhi[l][d] = ctx.soil.Ksat[l] * ctx.soil.phi_e[l] * pow(theta_w, 2.5 + ctx.soil.b[l]);                                            // ctx.soil.Ks times ctx.soil.soil_phi3D, computed directly as the exact power of ctx.params.theta.
                if (isnan(ctx.soil.soil_phi3D[l][d]) || isnan(ctx.soil.Ks[l][d]) || isnan(ctx.soil.KsPhi[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0) //|| ctx.soil.KsPhi[l][d]==0.0 || ctx.soil.Ks[l][d]==0.0 || ctx.soil.soil_phi3D[l][d]==0.0)
                    cout << "In bucket model, ctx.diag.layer " << l << " dcell " << d << " theta_w=" << theta_w << " ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) << " ctx.soil.soil_phi3D[l][d]=" << ctx.soil.soil_phi3D[l][d] << " ctx.soil.Ksat=" << ctx.soil.Ksat[l] << " ctx.soil.phi_e=" << ctx.soil.phi_e[l] << " ctx.soil.b[l]=" << ctx.soil.b[l] << " ctx.soil.KsPhi[l][d]=" << ctx.soil.KsPhi[l][d] << " ctx.soil.Ks[l][d]=" << ctx.soil.Ks[l][d] << endl;
            }
        }
    }

    // compute LAID for RecruitTree function
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        for (int sbsite = 0; sbsite < ctx.grid.sites + 2 * ctx.grid.SBORD; sbsite++)
            ctx.field.LAI3D[h][sbsite] = 0.0;
    for (int site = 0; site < ctx.grid.sites; site++)
        CalcLAI(ctx, ctx.T[site]); // Each tree contribues to ctx.field.LAI3D
    for (int h = ctx.grid.HEIGHT; h > 0; h--)
    { // LAI is computed by summing LAI from the canopy top to the ground
        for (int site = 0; site < ctx.grid.sites; site++)
        {
            int sbsite = site + ctx.grid.SBORD;
            ctx.field.LAI3D[h - 1][sbsite] += ctx.field.LAI3D[h][sbsite];
        }
    }

    // added in v.3.1 to reduce transient behavior of forest on reinitialization when leaf parameters and tree dimensions are not known, but not used yet: does not massively improve the initial configuration of trees, as the main problem is in overlapping crowns and not how they allocate their leaf area; could be solved by a preceding Canopy Constructor reshuffle (Fischer et al. 2020, RSE), or by
    //*#################################################################*/
    //*### Update Leaf Area for more realistic initial configuration ###*/
    //*#################################################################*/
    // To do so, we order by height first and then allocate leaf area from top to bottom
    //     vector<float> heights_trees;
    //     vector<int> sites_trees;
    //
    //     heights_trees.reserve(ctx.grid.sites);
    //     sites_trees.reserve(ctx.grid.sites);
    //
    //     for(int site = 0; site < ctx.grid.sites; site++){
    //        // Only consider non-initialized trees (i.e. t_LA < 0.0)
    //        if(ctx.T[site].t_age > 0.0 & ctx.T[site].t_LA < 0.0){
    //            float height = ctx.T[site].t_height;
    //            heights_trees.push_back(height);
    //            sites_trees.push_back(site);
    //
    //            int index_tree_current = int(sites_trees.size()) - 1;
    //
    //            while(index_tree_current > 0 && height > heights_trees[index_tree_current - 1]){
    //                heights_trees[index_tree_current] = heights_trees[index_tree_current - 1];
    //                sites_trees[index_tree_current] = sites_trees[index_tree_current - 1];
    //
    //                heights_trees[index_tree_current - 1] = height;
    //                sites_trees[index_tree_current - 1] = site;
    //                index_tree_current--;
    //            }
    //        }
    //     }
    //
    //     // clear voxel field
    //     for(int h=0;h<(ctx.grid.HEIGHT+1);h++)
    //         for(int sbsite=0;sbsite<ctx.grid.sites+2*ctx.grid.SBORD;sbsite++)
    //             ctx.field.LAI3D[h][sbsite] = 0.0;
    //
    //     // allocate and compute leaf area
    //     for(int index_site = 0; index_site < sites_trees.size(); index_site++){
    //        int site = sites_trees[index_site];
    //        //cout << site << " Site of tree: " << ctx.T[site].t_site << " Height: " << ctx.T[site].t_height << " Height from index: " << heights_trees[index_site] << endl;
    //        ctx.T[site].CalcLAinitial();
    //     }
}

// ######################################
//  Global function: Field dynamic memory allocation
// ######################################
void AllocMem(Context &ctx)
{
    // this needs better commenting and probably a rethink (mainly has to do with MPI)
    float d = 0.0; // maximum diameter possible
    for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
    {
        d = fmaxf(d, ctx.S[spp].s_dbhmax * 1.5);
    }
    float r = 25.0; // simply set to maximum crown radius possible in simulations

    ctx.grid.RMAX = int(r + ctx.params.p_nonvert * ctx.grid.NH * ctx.grid.LV * ctx.grid.HEIGHT);
    //  ctx.grid.RMAX = int(r);
    ctx.grid.SBORD = ctx.grid.cols * ctx.grid.RMAX;
    ctx.grid.dbhmaxincm = int(100. * d);
    cout << "ctx.grid.SBORD: " << ctx.grid.SBORD << endl;
    if (!mpi_rank)
    {
        // cout << "ctx.grid.HEIGHT : " << ctx.grid.HEIGHT << " ctx.grid.RMAX : " << ctx.grid.RMAX << " DBH : " << DBH <<"\n"; cout.flush();
        if (ctx.grid.RMAX > ctx.grid.rows)
        {
            // Consistency tests
            cerr << "Error : ctx.grid.RMAX > ctx.grid.rows \n";
            exit(-1);
        }
        if (ctx.grid.HEIGHT > ctx.grid.rows)
        {
            cerr << "Error : ctx.grid.HEIGHT > ctx.grid.rows \n";
            exit(-1);
        }
    }

    //** Initialization of dynamic Fields **
    //**************************************
    if (NULL == (ctx.diag.nbdbh = new int[ctx.grid.dbhmaxincm]))
        cerr << "!!! Mem_Alloc\n"; // Field for DBH histogram
    if (NULL == (ctx.diag.layer = new float[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n"; // Field for variables averaged by vertical ctx.diag.layer
#ifdef Output_ABC
    if (NULL == (abundances_species = new int[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // vector to save species abundances every recorded ctx.time.timestep
    if (NULL == (abundances_species10 = new int[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // vector to save species abundances every recorded ctx.time.timestep (dbh > 10cm)
    if (NULL == (biomass_species = new float[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // vector to save species abundances every recorded ctx.time.timestep (dbh > 10cm)

    if (NULL == (traits_species = new float *[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // vector to save species traits every recorded ctx.time.timestep

    for (int spp = 0; spp < (ctx.grid.nbspp + 1); spp++)
        if (NULL == (traits_species[spp] = new float[10]))
            cerr << "!!! Mem_Alloc\n";
    if (NULL == (traits_species10 = new float *[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // vector to save species traits every recorded ctx.time.timestep  (dbh > 10cm)
    for (int spp = 0; spp < (ctx.grid.nbspp + 1); spp++)
        if (NULL == (traits_species10[spp] = new float[10]))
            cerr << "!!! Mem_Alloc\n";
#endif
    if (NULL == (ctx.species.SPECIES_GERM = new int[ctx.grid.nbspp + 1]))
        cerr << "!!! Mem_Alloc\n"; // Field for democratic seed germination
    if (NULL == (ctx.species.SPECIES_SEEDS = new int *[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";                              // Field of seeds
    for (int site = 0; site < ctx.grid.sites; site++)                    // For each processor, we define a stripe above (labelled 0) and a stripe below (1). Each stripe is ctx.grid.SBORD in width.
        if (NULL == (ctx.species.SPECIES_SEEDS[site] = new int[ctx.grid.nbspp + 1])) // ALL the ctx.grid.sites need to be updated.
            cerr << "!!! Mem_Alloc\n";
    for (int site = 0; site < ctx.grid.sites; site++)
        for (int spp = 0; spp <= ctx.grid.nbspp; spp++)
            ctx.species.SPECIES_SEEDS[site][spp] = 0;
    if (NULL == (ctx.species.p_seed = new double[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.species.n_seed = new unsigned int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";

    double prob_seed = 1.0 / double(ctx.grid.sites);
    for (int i = 0; i < ctx.grid.sites; i++)
        ctx.species.p_seed[i] = prob_seed;
    for (int i = 0; i < ctx.grid.sites; i++)
        ctx.species.n_seed[i] = 0;

    if (NULL == (ctx.species.p_species = new double[ctx.grid.nbspp]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.species.n_species = new unsigned int[ctx.grid.nbspp]))
        cerr << "!!! Mem_Alloc\n";

    for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
    {
        double prob_species = double(ctx.S[spp].s_nbext);
        // cout << "prob_species: " << prob_species << endl;
        ctx.species.p_species[spp - 1] = prob_species;
    }

    for (int i = 0; i < ctx.grid.nbspp; i++)
        ctx.species.n_species[i] = 0;

    if (ctx.opt._SEEDTRADEOFF)
        if (NULL == (ctx.species.PROB_S = new float[ctx.grid.nbspp + 1]))
            cerr << "!!! Mem_Alloc\n";
    if (ctx.opt._NDD)
        if (NULL == (ctx.species.PROB_S = new float[ctx.grid.nbspp + 1]))
            cerr << "!!! Mem_Alloc\n";
    //  if (NULL==(persist=new long int[ctx.time.nbiter])) cerr<<"!!! Mem_Alloc\n";                  // Field for persistence
    //  if (NULL==(distr=new int[ctx.grid.cols])) cerr<<"!!! Mem_Alloc\n";

    if (NULL == (ctx.field.LAI3D = new float *[ctx.grid.HEIGHT + 1]))             // Field 3D
        cerr << "!!! Mem_Alloc\n";                             // Trees at the border of the simulated forest need to know the canopy occupancy by trees in the neighboring processor.
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)                     // For each processor, we define a stripe above (labelled 0) and a stripe below (1). Each stripe is ctx.grid.SBORD in width.
        if (NULL == (ctx.field.LAI3D[h] = new float[ctx.grid.sites + 2 * ctx.grid.SBORD])) // ALL the ctx.grid.sites need to be updated.
            cerr << "!!! Mem_Alloc\n";
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        for (int site = 0; site < ctx.grid.sites + 2 * ctx.grid.SBORD; site++)
            ctx.field.LAI3D[h][site] = 0.0;
    if (NULL == (ctx.field.Thurt[0] = new unsigned short[3 * ctx.grid.sites])) // Field for treefall impacts
        cerr << "!!! Mem_Alloc\n";
    for (int i = 1; i < 3; i++)
        if (NULL == (ctx.field.Thurt[i] = new unsigned short[ctx.grid.sites]))
            cerr << "!!! Mem_Alloc\n";

#ifdef WATER
    if (NULL == (ctx.soil.SWC3D = new float *[ctx.soil.nblayers_soil]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.soil_phi3D = new float *[ctx.soil.nblayers_soil]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Ks = new float *[ctx.soil.nblayers_soil]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.KsPhi = new float *[ctx.soil.nblayers_soil]))
        cerr << "!!! Mem_Alloc\n";
    // if (NULL==(KsPhi2=new float*[ctx.soil.nblayers_soil])) cerr<<"!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Transpiration = new float *[ctx.soil.nblayers_soil]))
        cerr << "!!! Mem_Alloc\n";
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        if (NULL == (ctx.soil.SWC3D[l] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
        if (NULL == (ctx.soil.soil_phi3D[l] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
        if (NULL == (ctx.soil.Ks[l] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
        if (NULL == (ctx.soil.KsPhi[l] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
        // if (NULL==(KsPhi2[l]=new float[ctx.grid.nbdcells])) cerr<<"!!! Mem_Alloc\n";
        if (NULL == (ctx.soil.Transpiration[l] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
        for (int dcell = 0; dcell < ctx.grid.nbdcells; dcell++)
        {
            // ctx.soil.SWC3D[l][dcell]=ctx.soil.Max_SWC[l];
            ctx.soil.SWC3D[l][dcell] = ctx.soil.FC_SWC[l];
            ctx.soil.soil_phi3D[l][dcell] = 0.0;
            ctx.soil.Ks[l][dcell] = 0.0;
            ctx.soil.KsPhi[l][dcell] = 0.0;
            // KsPhi2[l][dcell]=0.0;
            ctx.soil.Transpiration[l][dcell] = 0.0;
            if (ctx.soil.SWC3D[l][dcell] <= 0.0)
            {
                cout << ctx.soil.SWC3D[l][dcell] << "\t" << ctx.soil.Max_SWC[l] << "\n";
            }
        }
    }
    if (NULL == (ctx.soil.LAI_DCELL = new float *[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        if (NULL == (ctx.soil.LAI_DCELL[h] = new float[ctx.grid.nbdcells]))
            cerr << "!!! Mem_Alloc\n";
    }
    if (NULL == (ctx.soil.LAI_young = new float[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.LAI_mature = new float[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.LAI_old = new float[ctx.grid.HEIGHT + 1]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Canopy_height_DCELL = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.TopWindSpeed_DCELL = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.HSum_DCELL = new int[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Interception = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Throughfall = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Runoff = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Leakage = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    if (NULL == (ctx.soil.Evaporation = new float[ctx.grid.nbdcells]))
        cerr << "!!! Mem_Alloc\n";
    for (int dcell = 0; dcell < ctx.grid.nbdcells; dcell++)
    {
        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            ctx.soil.LAI_DCELL[h][dcell] = 0.0;
        ctx.soil.Canopy_height_DCELL[dcell] = 0.0;
        ctx.soil.HSum_DCELL[dcell] = 0;
        ctx.soil.TopWindSpeed_DCELL[dcell] = 0.0;
        ctx.soil.Interception[dcell] = 0.0;
        ctx.soil.Throughfall[dcell] = 0.0;
        ctx.soil.Runoff[dcell] = 0.0;
        ctx.soil.Leakage[dcell] = 0.0;
        ctx.soil.Evaporation[dcell] = 0.0;
    }
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        ctx.soil.LAI_young[h] = 0.0;
        ctx.soil.LAI_mature[h] = 0.0;
        ctx.soil.LAI_old[h] = 0.0;
    }

    if (NULL == (ctx.grid.site_DCELL = new int[ctx.grid.sites]))
        cerr << "!!! Mem_Alloc\n";
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        int x = site % ctx.grid.cols;
        int y = site / ctx.grid.cols;
        int dcol = x / ctx.grid.length_dcell;
        int drow = y / ctx.grid.length_dcell;
        ctx.grid.site_DCELL[site] = dcol + ctx.grid.linear_nb_dcells * drow;
    }
#endif

#ifdef MPI // Fields for MPI operations
    for (i = 0; i < 2; i++)
    {                                                             //  Two fields: one for the CL north (0) one for the CL south (1)
        if (NULL == (LAIc[i] = new unsigned short *[ctx.grid.HEIGHT + 1])) // These fields contain the light info in the neighboring procs (2*ctx.grid.SBORD in width, not ctx.grid.SBORD !). They are used to update local fields
            cerr << "!!! Mem_Alloc\n";
        for (h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            if (NULL == (LAIc[i][h] = new unsigned short[2 * ctx.grid.SBORD]))
                cerr << "!!! Mem_Alloc\n";
    }
#endif
}

//!  Free dynamic memory
void FreeMem(Context &ctx)
{
    delete[] ctx.diag.nbdbh;
    delete[] ctx.diag.layer;
    delete[] ctx.species.SPECIES_GERM;
#ifdef WATER
    delete[] ctx.grid.site_DCELL;
#endif
    for (int site = 0; site < ctx.grid.sites; site++)
        delete[] ctx.species.SPECIES_SEEDS[site];
    delete[] ctx.species.SPECIES_SEEDS;
    delete[] ctx.species.p_seed;
    delete[] ctx.species.n_seed;
    delete[] ctx.species.p_species;
    delete[] ctx.species.n_species;

    if (ctx.opt._SEEDTRADEOFF || ctx.opt._NDD)
        delete[] ctx.species.PROB_S;
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        delete[] ctx.field.LAI3D[h];
    delete[] ctx.field.LAI3D;

    for (int i = 0; i < 3; i++)
        delete[] ctx.field.Thurt[i];

#ifdef Output_ABC
    delete[] chm_field_previous;
    delete[] chm_field_previous_ALS;
    delete[] chm_field_current;
    delete[] chm_field_current_ALS;
    delete[] chm_field_changes;

    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        delete[] transmittance_simulatedALS[h];
        delete[] transmittance_direct[h];
        delete[] transmittance_simulatedALS_sampling[h];
    }

    delete[] transmittance_simulatedALS;
    delete[] transmittance_direct;
    delete[] transmittance_simulatedALS_sampling;

    delete[] abundances_species;
    delete[] abundances_species10;
    delete[] biomass_species;

    for (int spp = 0; spp < (ctx.grid.nbspp + 1); spp++)
    {
        delete[] traits_species[spp];
        delete[] traits_species10[spp];
    }
    delete[] traits_species;
    delete[] traits_species10;

#endif
    delete[] ctx.lookup.LookUp_T;
    delete[] ctx.lookup.LookUp_KmT;
    delete[] ctx.lookup.LookUp_VPD;
    delete[] ctx.lookup.LookUp_flux_absorption;
    delete[] ctx.lookup.LookUp_flux;
    delete[] ctx.lookup.LookUp_ExtinctLW;
    // delete [] LookUp_Rday; //newIM: no redundancy anymore between LookUp_Rday and ctx.lookup.LookUp_Rleaf
    delete[] ctx.lookup.LookUp_JmaxT;
    delete[] ctx.lookup.LookUp_Rstem;
    delete[] ctx.lookup.LookUp_Rleaf; // newIM: no redundancy anymore between LookUp_Rday and ctx.lookup.LookUp_Rleaf
    delete[] ctx.lookup.LookUp_GammaT;
    // delete [] LookUp_Rnight; //newIM: no redundancy anymore between LookUp_Rday and ctx.lookup.LookUp_Rleaf
    delete[] ctx.lookup.LookUp_VcmaxT;
#ifdef WATER
    delete[] ctx.lookup.LookUp_INLR;
    delete[] ctx.lookup.LookUp_SLOPE;
    delete[] ctx.lookup.LookUp_GRADN;
#endif
}
