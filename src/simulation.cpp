#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"


// ######################################
//  Global function: Evolution at each ctx.time.timestep
// ######################################
void Evolution(Context &ctx)
{
#ifdef CHECK_CARBON
    if (ctx.time.iter == 0)
    {
        carbon_assimilated_total = 0.0;
        carbon_net_total = 0.0;
    }
#endif

#ifdef WATER
    ctx.soil.transpiration_1016 = 0.0;
#endif

    UpdateField(ctx); // Update light fields and seed banks
    ctx.diag.nbtrees_n10 = ctx.diag.nbtrees_n30 = ctx.diag.nbdead_n1 = ctx.diag.nbdead_n10 = ctx.diag.nbdead_n30 = 0;
    ctx.diag.nbtrees_carbstarv_n1 = ctx.diag.nbtrees_carbstarv_n10 = ctx.diag.nbtrees_carbstarv_n30 = 0;

#ifdef Output_ABC
    nbdead_n10_abc = 0;
#endif
    if (ctx.opt._BASICTREEFALL)
    {
        // secondary treefalls are triggered first, since they have been caused in the previous iteration
        TriggerTreefallSecondary(ctx); // Compute and distribute Treefall events, caused by treefalls in the previous iteration
        TriggerTreefall(ctx);          // Compute and distribute Treefall events, caused by wind drag
    }

#ifdef WATER
    // Parallel phase: water availability and leaf fluxes of every living tree. These
    // dominate run time, read only fields that are fixed during the tree loop below
    // (updated in UpdateField), write only the tree's own state and ctx.resp_gpp, and draw
    // no random numbers. Mortality draws, deaths, growth and the plot-level sums stay in
    // the serial loop, in site order, so results do not depend on the number of threads.
    ctx.resp_gpp.resize(ctx.grid.sites);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        Tree &tree = ctx.T[site];
        if (tree.t_age)
        {
            tree.Water_availability(ctx);
            ctx.resp_gpp[site] = tree.ComputeRespGPP(ctx);
        }
    }
#endif

    for (int site = 0; site < ctx.grid.sites; site++)
    {
        //**** Tree evolution: Growth or death ****
#ifdef WATER
        ctx.T[site].Update(ctx, true);
#else
        ctx.T[site].Update(ctx);
#endif
    }

    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            ctx.soil.Transpiration[l][d] = 0.0;
        }
    }
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        ctx.T[site].Water_uptake(ctx);
    }

    // Update trees
    Average(ctx); //! Compute averages for outputs
    if (ctx.opt._OUTPUT_extended)
        OutputField(ctx); //! Output the statistics
}

// #################################
//  Global function: Compute field Seed
// #################################
void UpdateSeeds(Context &ctx)
{
    // With MPI option: Pass seeds across processors => two more fields to be communicated between n.n. (nearest neighbor) processors. NB: dispersal distance is bounded by the value of 'ctx.grid.rows'. At least 99 % of the seeds should be dispersed within the stripe or on the n.n. stripe. Hence ctx.grid.rows > 4.7*max(dist_moy_dissemination),for an exponential dispersal kernel.
    // dispersal only once a year
    if (ctx.time.iter % ctx.time.iterperyear == 0)
    {
        // acceleration, using the multinomial distribution
        // Plot area in hectares. This was an integer division, which gave zero seed rain
        // for plots smaller than 1 ha and truncated non-integer areas (e.g. 1.5 ha -> 1 ha).
        // For whole-hectare plots the result is unchanged.
        double ha = double(ctx.grid.sites) / 10000.0;
        unsigned int nb_seedrain = (unsigned int)(ctx.params.Cseedrain * ha);
        gsl_ran_multinomial(ctx.rng.gslrand, ctx.grid.sites, nb_seedrain, ctx.species.p_seed.data(), ctx.species.n_seed.data());
        cout << ctx.grid.sites << " Seedrain: " << nb_seedrain << endl;
        int seedsadded = 0;
        for (int s = 0; s < ctx.grid.sites; s++)
        {
            // if(ctx.T[s].t_age == 0){
            int nbseeds = ctx.species.n_seed[s];
            // cout << "Site: " << s << " nbseeds: " << nbseeds << " ctx.grid.nbspp: " << ctx.grid.nbspp << endl;
            gsl_ran_multinomial(ctx.rng.gslrand, ctx.grid.nbspp, nbseeds, ctx.species.p_species.data(), ctx.species.n_species.data());
            for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
            {
                int nbseeds_species = ctx.species.n_species[spp - 1];
                // cout << "Site: " << s << " Species: " << spp << " nbseeds: " << nbseeds_species << endl;
                if (nbseeds_species > 0)
                {
                    ctx.species.SPECIES_SEEDS[s][spp] = 1;
                    seedsadded++;
                }
                else
                    ctx.species.SPECIES_SEEDS[s][spp] = 0;
            }
            //}
        }

        // now disperse seeds from the trees on site
        // dispersion comes after seedrain calculation, because seedrain automatically removes seeds from places where there was no incoming seed
        int trees_mature = 0;
        for (int site = 0; site < ctx.grid.sites; site++)
        { // disperse seeds produced by mature trees
            if (ctx.T[site].t_age)
            {
                if (ctx.T[site].t_dbh >= ctx.T[site].t_dbhmature)
                    trees_mature++;
                ctx.T[site].DisperseSeed(ctx);
            }
        }

        int nbspecies_affected = 0;
        int seedsadded_effective = 0;
        for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
        {
            int seedsadded_species = 0;
            for (int s = 0; s < ctx.grid.sites; s++)
            {
                seedsadded_species += ctx.species.SPECIES_SEEDS[s][spp];
            }
            if (seedsadded_species > 0)
                nbspecies_affected++;
            seedsadded_effective += seedsadded_species;
        }
        cout << "Trees_mature: " << trees_mature << " Nbseedsadded: " << seedsadded << " effective: " << seedsadded_effective << " nbspeciesaffected: " << nbspecies_affected << endl;
    }
}

// #################################
//  Global function: Update all fields
// #################################
//! - This is an important function for TROLL -- Includes many of the operations
//! - set the iteration environment -- nb: the current structure of code suppose that environment is periodic (a period = a year), if one wants to input a variable climate, with interannual variation and climate change along the simulation, a full climatic input needs to be input (ie number of columns=ctx.time.iter and not ctx.time.iterperyear) and change ctx.time.iterperyear by ctx.time.nbiter here.
void UpdateField(Context &ctx)
{

#ifdef FULL_CLIMATE

    ctx.climate.tnight = ctx.climate.NightTemperature[ctx.time.iter % ctx.time.nbdays];
    ctx.climate.precip = ctx.climate.Rainfall[ctx.time.iter % ctx.time.nbdays];
    ctx.climate.WSDailyMean = ctx.climate.DailyMeanWindSpeed[ctx.time.iter % ctx.time.nbdays];
    ctx.climate.WDailyMean = ctx.climate.DailyMeanIrradiance[ctx.time.iter % ctx.time.nbdays] * ctx.params.SWtoPPFD;
    ctx.climate.tDailyMean = ctx.climate.DailyMeanTemperature[ctx.time.iter % ctx.time.nbdays];
    ctx.climate.VPDDailyMean = ctx.climate.DailyMeanVapourPressureDeficit[ctx.time.iter % ctx.time.nbdays];

#else
    ctx.climate.tnight = ctx.climate.NightTemperature[ctx.time.iter % ctx.time.iterperyear];
    ctx.climate.precip = ctx.climate.Rainfall[ctx.time.iter % ctx.time.iterperyear];
    ctx.climate.WSDailyMean = ctx.climate.DailyMeanWindSpeed[ctx.time.iter % ctx.time.iterperyear];
    ctx.climate.WDailyMean = ctx.climate.DailyMeanIrradiance[ctx.time.iter % ctx.time.iterperyear] * ctx.params.SWtoPPFD;
    ctx.climate.tDailyMean = ctx.climate.DailyMeanTemperature[ctx.time.iter % ctx.time.iterperyear];
    ctx.climate.VPDDailyMean = ctx.climate.DailyMeanVapourPressureDeficit[ctx.time.iter % ctx.time.iterperyear];

#endif // FULL_CLIMATE

    UpdateSeeds(ctx);

    if (ctx.opt._NDD)
    {
        // Evolution of the field NDDfield

        float normBA = 10000.0 / (0.001 + PI * ctx.crown.Rndd * ctx.crown.Rndd * ctx.crown.BAtot);
        for (int site = 0; site < ctx.grid.sites; site++)
        {

            for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
            {
                // if ((ctx.time.iter == int(ctx.time.nbiter-1))&&(site>80000)&&(site<85000))  { sor[142]<< ctx.T[site].t_NDDfield[spp] << "\t" ;}
                ctx.T[site].t_NDDfield[spp] = 0;
            }
            // if (ctx.time.iter == int(ctx.time.nbiter-1))  sor[142]<< "\n";

            int row0 = ctx.T[site].t_site / ctx.grid.cols;
            int col0 = ctx.T[site].t_site % ctx.grid.cols;
            for (int col = max(0, int(col0 - ctx.crown.Rndd)); col <= min(ctx.grid.cols - 1, int(col0 + ctx.crown.Rndd)); col++)
            {
                for (int row = max(0, int(row0 - ctx.crown.Rndd)); row <= min(ctx.grid.rows - 1, int(row0 + ctx.crown.Rndd)); row++)
                { // loop over the neighbourhood
                    int xx = col0 - col;
                    int yy = row0 - row;
                    float d = sqrt(xx * xx + yy * yy);
                    if ((d <= ctx.crown.Rndd) && (d > 0))
                    { // is the voxel within the neighbourhood?
                        int j = ctx.grid.cols * row + col;
                        if (ctx.T[j].t_age)
                            ctx.T[site].t_NDDfield[ctx.T[j].t_sp_lab] += PI * ctx.T[j].t_dbh * ctx.T[j].t_dbh * 0.25 * normBA;
                    }
                }
            }
        }
    }

    RecruitTree(ctx);

    //  Compute Field ctx.field.LAI3D

    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        for (int sbsite = 0; sbsite < ctx.grid.sites + 2 * ctx.grid.SBORD; sbsite++)
            ctx.field.LAI3D[h][sbsite] = 0.0;
    for (int site = 0; site < ctx.grid.sites; site++)
        ctx.T[site].CalcLAI(ctx); // Each tree contribues to ctx.field.LAI3D

    for (int h = ctx.grid.HEIGHT; h > 0; h--)
    { // LAI is computed by summing LAI from the canopy top to the ground
        for (int site = 0; site < ctx.grid.sites; site++)
        {
            int sbsite = site + ctx.grid.SBORD;
            ctx.field.LAI3D[h - 1][sbsite] += ctx.field.LAI3D[h][sbsite];
        }
    }

#ifdef WATER
    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        ctx.soil.Runoff[d] = 0.0;
        ctx.soil.Interception[d] = 0.0;
        ctx.soil.Throughfall[d] = 0.0;
        ctx.soil.Evaporation[d] = 0.0;
        ctx.soil.Leakage[d] = 0.0;
        // for (int l=0;l<ctx.soil.nblayers_soil;l++) {
        //     ctx.soil.Transpiration[l][d]=0.0;
        // }
        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        {
            ctx.soil.LAI_DCELL[h][d] = 0.0;
        }
        ctx.soil.Canopy_height_DCELL[d] = 0.0;
        ctx.soil.HSum_DCELL[d] = 0;
        ctx.soil.TopWindSpeed_DCELL[d] = 0.0;
    }

    // ctx.soil.LAI_DCELL[h][dcell] provide the average LAI at height h in dcell, to estimate water interception and evaporation in each dcell at each ctx.time.timestep, as well as wind speed for a given height (this latter is new and was added simultaneously to BOUNDARY LAYER_ITERATIVE_SCHEME (IM June 2021).
    for (int site = 0; site < ctx.grid.sites; site++)
    {

        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
        {
            ctx.soil.LAI_DCELL[h][ctx.grid.site_DCELL[site]] += ctx.field.LAI3D[h][site + ctx.grid.SBORD];
        }

        float Htop = 0.0;
        int H = ctx.grid.HEIGHT;
        while (Htop == 0.0 && H > 0)
        {
            if (ctx.field.LAI3D[H][site + ctx.grid.SBORD] > 0.0)
            {
                Htop = H;
            }
            H--;
        }

        if (Htop > 0)
        {
            ctx.soil.Canopy_height_DCELL[ctx.grid.site_DCELL[site]] += Htop;
            ctx.soil.HSum_DCELL[ctx.grid.site_DCELL[site]]++;
        }
    }

    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {

        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            ctx.soil.LAI_DCELL[h][d] *= ctx.grid.i_sites_per_dcell;
        if (ctx.soil.HSum_DCELL[d] > 0)
        {
            ctx.soil.Canopy_height_DCELL[d] *= 1.0 / float(ctx.soil.HSum_DCELL[d]);
        }
        else
        {
            ctx.soil.Canopy_height_DCELL[d] = 0.0;
        }

#ifdef FULL_CLIMATE
        if (ctx.soil.Canopy_height_DCELL[d] <= MeteoStation_Height)
        {
            ctx.soil.TopWindSpeed_DCELL[d] = 1.204 / log(16.67 * ((MeteoStation_Height / ctx.soil.Canopy_height_DCELL[d]) - 0.8)); // WS is the ctx.time.timestep windspeed at a height=MeteoStation_Height, and ctx.soil.TopWindSpeed_DCELL is the wind speed computed at a height=ctx.soil.Canopy_height_DCELL[d], according to the model of Monteith & Unsworth 2008 (see Rau et al's TROLL manuscript), with d=0.8H and z0=0.06H; 16.67~1/0.06, 1.204=log(0.2/0.06).
        }
        else
            ctx.soil.TopWindSpeed_DCELL[d] = exp(alphaInoue * (1 - MeteoStation_Height / ctx.soil.Canopy_height_DCELL[d]));
        if (ctx.soil.Canopy_height_DCELL[d] == 0)
        {
            cout << "in UpdateField: d=" << d << "; Canopyheight_DCELL[d]=" << ctx.soil.Canopy_height_DCELL[d] << "; ctx.soil.HSum_DCELL[d]=" << ctx.soil.HSum_DCELL[d] << "; ctx.soil.TopWindSpeed_DCELL[d]=" << ctx.soil.TopWindSpeed_DCELL[d] << endl;
        }
#else
        if (ctx.soil.Canopy_height_DCELL[d] <= MeteoStation_Height)
        {
            ctx.soil.TopWindSpeed_DCELL[d] = ctx.climate.WSDailyMean * 1.204 / log(16.67 * ((MeteoStation_Height / ctx.soil.Canopy_height_DCELL[d]) - 0.8)); // WS is the ctx.time.timestep windspeed at a height=MeteoStation_Height, and ctx.soil.TopWindSpeed_DCELL is the wind speed computed at a height=ctx.soil.Canopy_height_DCELL[d], according to the model of Monteith & Unsworth 2008 (see Rau et al's TROLL manuscript), with d=0.8H and z0=0.06H; 16.67~1/0.06, 1.204=log(0.2/0.06).
        }
        else
            ctx.soil.TopWindSpeed_DCELL[d] = ctx.climate.WSDailyMean * exp(alphaInoue * (1 - MeteoStation_Height / ctx.soil.Canopy_height_DCELL[d]));
#endif

        // cout << "iter=" << ctx.time.iter << " ctx.climate.WSDailyMean=" << ctx.climate.WSDailyMean << " d=" << d <<  " canopy_height_DCELL[d]=" << ctx.soil.Canopy_height_DCELL[d] << " ctx.soil.HSum_DCELL[d]=" << ctx.soil.HSum_DCELL[d] <<" ctx.soil.TopWindSpeed_DCELL[d]=" << ctx.soil.TopWindSpeed_DCELL[d];
        // if (d==225) cout << " d=" << d <<  " canopy_height_DCELL[d]=" << ctx.soil.Canopy_height_DCELL[d] ;
        // cout << endl;
    }
#endif


#ifdef WATER
    //**  Evolution of belowground hydraulic fields: Soil bucket model

    // for(int site=0;site<ctx.grid.sites;site++) ctx.T[site].Water_uptake(); // Update of ctx.soil.Transpiration: tree water uptake, each tree will deplete soil water content through its transpiration. Now made ate the end of the evolution loop so that the outputs for water uptake match the others (otherwise lag of one ctx.time.timestep)

    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        //****   BUCKET MODEL in each dcell   ****
        // the unit used for water volume throughout the bucket model is m3.
        // NOTE: under the assumption of a flat terrain and no lateral fluxes, as it is assumed here for a first implementation, the order with which dcells are visited during the loop does not matter. With topography, we will need to visit the soil voxels (ie. dcells*ctx.diag.layer) from highest to lowest elevation so that run-off from highest voxels contribute to the water flux entering the lowest ones.
        //  to be investigated: does the order in which transpiration and evaporation are retrieved from the soil affect the overall outcome? which one should be retrieved first?

        // Water uptake through tree transpiration
        float w_uptake = 0.0;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            w_uptake = fminf(ctx.soil.Transpiration[l][d], (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]));
            if (ctx.soil.Transpiration[l][d] < 0.0 || isnan(ctx.soil.Transpiration[l][d]) || isnan(w_uptake) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0)
            {
                cout << "l=" << l << " d= " << d << " transpiration=" << ctx.soil.Transpiration[l][d] << " and ctx.soil.SWC3D[l][d]=" << ctx.soil.SWC3D[l][d] << " and ctx.soil.Min_SWC[l]=" << ctx.soil.Min_SWC[l] << " and ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l] << endl;
            }
            ctx.soil.SWC3D[l][d] -= w_uptake;
            ctx.soil.SWC3D[l][d] = fmaxf(ctx.soil.SWC3D[l][d], ctx.soil.Min_SWC[l]);
            if (ctx.soil.Transpiration[l][d] < 0.0 || isnan(ctx.soil.Transpiration[l][d]) || isnan(w_uptake) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0)
            {
                cout << "After w_uptake, l=" << l << " d= " << d << " transpiration=" << ctx.soil.Transpiration[l][d] << " and ctx.soil.SWC3D[l][d]=" << ctx.soil.SWC3D[l][d] << " and ctx.soil.Min_SWC[l]=" << ctx.soil.Min_SWC[l] << " and ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l] << endl;
            }
        }

        // ctx.soil.Evaporation from soil

        // if this should be negligible in dense forest understory, it should have a more important effect in open areas, especially through species filtering at germination stage, at the beginning of a succession or in gaps in drier conditions; see Marthews et al. 2008 Ecological Modelling
        // However it is sometimes neglected and not represented in models, eg. Laio et al. 2001, Guterriez et al. 2014, Fischer et al. 2014.
        // note that, in this version, evaporation only depletes the most superficial soil ctx.diag.layer. This could be changed, especially if the superficial soil ctx.diag.layer is particularly thin and the energy reaching the soil high.

        // here, we use a phenomenological approach, following Granier et al. 1999 Ecological Modelling and Wagner et al. 2011 AFM, which assumed that evaporation is proportional to the energy reaching the soil.[this is an approximation as as the soil gets drier, more energy would be needed to remove the same amount of water from the soil as water molecules should be more tighly bound to soil particules and cavitation also occur in the soil...] ==> see if a model under which evaporation also depends on the soil water potential would not be better -- I guess so.
        // parameter values are not so clear, so TO BE CHECKED.
        // float e_factor=ctx.params.PPFDtoSW * 3600*0.000001*ctx.time.nbhours_covered* 0.1 * ctx.grid.sites_per_dcell*ctx.grid.LH*ctx.grid.LH*0.001; // to be moved outside of the loop to avoid repeating calculation.
        // float e_Granier = e_factor* ctx.climate.WDailyMean * exp(-ctx.params.klight*ctx.soil.LAI_DCELL[0][d]);
        // 3600*0.000001*ctx.time.nbhours_covered to convert Wmax in micromol of PAR /s /m2 into  Joule, and 10^-6 to MJoule as in Wagner et al. 2011 (however the value provided by Wagner et al. 2011 seems really weird -too high-, and the values we obtained here are in agreement with the ones reported in Marthews et al. 2014.
        // the value 0.1 is drawn from Wagner et al. 2011, but not really explained... to be checked!
        // ctx.grid.sites_per_dcell*ctx.grid.LH*ctx.grid.LH*0.001 is to convert the amount of water in mm, ie. in 10-3 m3/m2, to the amount of water evaporated for the focal dcell in m3

        // in this newer version, we used the framework provided by Sellers et al. 1992, which is better mechanistically grounded: depends on the soil ctx.diag.layer resistance, which varies with its water potential, and the aerodynamic resistance in series and the differences of vapour pressure between the top soil ctx.diag.layer and air just above
        float absorb_prev = ctx.soil.LAI_DCELL[1][d];
        float absorb_current = ctx.field.LAI3D[0][d];
        float absorb_delta = absorb_current - absorb_prev;
        if (absorb_delta < 0.0)
            absorb_delta = 0.0; //! eliminate rounding errors
        int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);
        float VPDground = ctx.climate.VPDDailyMean * ctx.lookup.LookUp_VPD[intabsorb] * 1000; // in Pa
        float Tsoil = ctx.climate.tDailyMean - ctx.lookup.LookUp_T[intabsorb];
        float esat_ground = 611.21 * exp((18.678 - (Tsoil / 234.5)) * (Tsoil / (257.14 + Tsoil))); // Buck equation; in Pa (see Jones p. 348)
        float esoil = esat_ground * exp(2.17 * ctx.soil.soil_phi3D[0][d] / (Tsoil - ABSZERO));              // esoil variation with the top soil ctx.diag.layer water potential, following Duursma & Medlyn 2012 equ. 17, Cochard et al. 2021 equ. 36., see equ. 5.14 in Jones (p. 102), in Pa
        float eair = esat_ground - VPDground;                                                      // in Pa
        // float r_soil = exp(8.206 - 4.255*ctx.soil.SWC3D[0][d]/ctx.soil.Max_SWC[0]) ; // soil surface resistance in s ctx.params.m-1, following Sellers et al. 1992 equ. 19, see also equ 12 in Merlin et al. 2016 (also used in CLM, Oleson et al. 2007).
        float r_soil = exp(8.206 - 4.255 * ctx.soil.SWC3D[0][d] / ctx.soil.FC_SWC[0]); // soil surface resistance in s ctx.params.m-1, following Sellers et al. 1992 equ. 19, see also equ 12 in Merlin et al. 2016 (also used in CLM, Oleson et al. 2007).
#ifdef FULL_CLIMATE
        float r_aero = 43.17347 * exp(alphaInoue * (1 - 1 / ctx.soil.Canopy_height_DCELL[d])) / (ctx.climate.WSDailyMean * ctx.soil.TopWindSpeed_DCELL[d]); // aerodynamic resistance to hear transfer (boundary ctx.diag.layer just above the soil surface), in s ctx.params.m-1 (see equ. 7 and 14 in Duursma & Medlyn 2012; and equ. B10 in Merlin et al. 2016). 43.17347= log(1/0.001)/(0.40*0.40), where 1= the reference height where the wind speed is measured, in ctx.params.m, 0.001=the momentum soil roughness in ctx.params.m (set to 0.001 following Yang et al. 2008 and Stefan et al 2015 in Merlin et al. 2016 equ B10), and 0.40=the von Karman constant.
#else
        float r_aero = 43.17347 * exp(alphaInoue * (1 - 1 / ctx.soil.Canopy_height_DCELL[d])) / ctx.soil.TopWindSpeed_DCELL[d]; // aerodynamic resistance to hear transfer (boundary ctx.diag.layer just above the soil surface), in s ctx.params.m-1 (see equ. 7 and 14 in Duursma & Medlyn 2012; and equ. B10 in Merlin et al. 2016). 43.17347= log(1/0.001)/(0.40*0.40), where 1= the reference height where the wind speed is measured, in ctx.params.m, 0.001=the momentum soil roughness in ctx.params.m (set to 0.001 following Yang et al. 2008 and Stefan et al 2015 in Merlin et al. 2016 equ B10), and 0.40=the von Karman constant.
#endif
        float Rtot = r_soil + r_aero;                                                                                 // in s ctx.params.m-1
        float e = ctx.time.nbhours_covered * ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * 0.0078 * (esoil - eair) / ((Tsoil - ABSZERO) * Rtot); // 0.0078 = 0.001*3600*18e-3/8.31 with 18e-3 = the molar mass of water vapor in kg/mol and 8.31 the ideal gas constant in J/mol/K; 0.001*3600*ctx.time.nbhours_covered*ctx.grid.sites_per_dcell*ctx.grid.LH*ctx.grid.LH is used to convert evaporation in kg ctx.params.m-2 s-1 to m3 per day per dcell.

        // if (ctx.soil.soil_phi3D[0][d] < -1) {
        //  cout << "r_soil=" << r_soil << " r_soil_sellers=" << r_soil_sellers <<" r_aero=" <<r_aero << " VPDground=" << VPDground << " esat_ground=" << esat_ground << " esoil=" << esoil << " eair=" << eair << " ctx.soil.soil_phi3D[0][d]=" << ctx.soil.soil_phi3D[0][d] << " Tsoil=" << Tsoil << " ctx.soil.TopWindSpeed_DCELL[d]=" << ctx.soil.TopWindSpeed_DCELL[d] << " Wind ground level=" << exp(-alphaInoue*(1-1/ctx.soil.Canopy_height_DCELL[d]))*ctx.soil.TopWindSpeed_DCELL[d] << " evaporation S92=" << e  << " evaporation S92_sellers=" << e_sellers << " e_granier=" << e_Granier << " ctx.soil.SWC3D[0][d]-ctx.soil.Min_SWC[0]=" << ctx.soil.SWC3D[0][d]-ctx.soil.Min_SWC[0] << endl;
        //  }

        ctx.soil.Evaporation[d] = fmaxf(0.0, fminf(e, (ctx.soil.SWC3D[0][d] - ctx.soil.Min_SWC[0]))); // the amount of water evaporated from the soil cannot result in a water content below the residual water content. A model depending on soil matric potential would not need this.
        if (ctx.soil.Evaporation[d] < 0 || isnan(ctx.soil.Evaporation[d]) || (ctx.soil.SWC3D[0][d] - ctx.soil.Min_SWC[0]) < 0)
        {
            cout << "evaporation=" << ctx.soil.Evaporation[d] << " and e=" << e << " and ctx.soil.SWC3D[0][d]=" << ctx.soil.SWC3D[0][d] << " and ctx.soil.Min_SWC[0]=" << ctx.soil.Min_SWC[0] << "and ctx.soil.SWC3D[0][d]-ctx.soil.Min_SWC[0]=" << ctx.soil.SWC3D[0][d] - ctx.soil.Min_SWC[0] << "; raero=" << r_aero << "; ctx.soil.Canopy_height_DCELL[d]=" << ctx.soil.Canopy_height_DCELL[d] << "; ctx.soil.TopWindSpeed_DCELL[d]=" << ctx.soil.TopWindSpeed_DCELL[d] << endl;
        }

        ctx.soil.SWC3D[0][d] -= ctx.soil.Evaporation[d];

        // Refilling by rainfall

        ctx.soil.Interception[d] = fminf(ctx.climate.precip, 0.2 * ctx.soil.LAI_DCELL[0][d]); // This is the amount of rainfall - in mm, as rainfall -, intercepted by vegetation cover, following the approach used in Liang et al. 1994 Journal of Geophysical Reserach, and also used by Laio et al. 2001 Advances in Water Resources and Fischer et al. 2014 Environmental Modelling & Software (FORMIX3, Madagascar). More complex approach can be used however - see eg. Gutierrez et al. 2014 Plos One (FORMIND, Chili), or Wagner et al. 2011 AFM (Paracou)
        ctx.soil.Throughfall[d] = ctx.climate.precip - ctx.soil.Interception[d];
        ctx.soil.Throughfall[d] *= ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * 0.001; // to convert in absolute amount of water entering the soil voxel in m3

        if (isnan(ctx.soil.Throughfall[d]) || (ctx.soil.Throughfall[d]) < 0)
        {
            cout << "Incorrect throughfall" << endl;
            cout << ctx.climate.precip << "\t" << ctx.soil.Interception[d] << "\t" << ctx.soil.LAI_DCELL[0][d] << endl;
        }

        float in = ctx.soil.Throughfall[d];

        /*if(ctx.soil.SWC3D[0][d]<ctx.soil.Max_SWC[0]) {
            int l=0;
            while((l<ctx.soil.nblayers_soil) && (in>0.0)) {
                if(in>(ctx.soil.Max_SWC[l]-ctx.soil.SWC3D[l][d])) {
                    in-=(ctx.soil.Max_SWC[l]-ctx.soil.SWC3D[l][d]);
                    ctx.soil.SWC3D[l][d]=ctx.soil.Max_SWC[l];
                    if(isnan(ctx.soil.SWC3D[l][d]) || (ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l])<=0) {
                        cout << "incorrect ctx.soil.SWC3D, Min/ctx.soil.Max_SWC" << endl;
                        cout <<ctx.soil.Max_SWC[l] << endl;
                    }
                }
                else{
                    ctx.soil.SWC3D[l][d]+=in;
                    if (isnan(ctx.soil.SWC3D[l][d]) || (ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l])<0) {
                        cout << "incorrect ctx.soil.SWC3D, Min/ctx.soil.Max_SWC" << endl;
                        cout << ctx.soil.Throughfall[d] << "\t" <<in <<"\t" <<  ctx.climate.precip << "\t" << ctx.soil.Interception[d] << "\t" << ctx.soil.LAI_DCELL[0][d] << endl;
                    }
                    in=0.0;
                }
                l++;
            }
        }*/
        if (ctx.soil.SWC3D[0][d] < ctx.soil.Max_SWC[0])
        {
            int l = 0;
            while ((l < ctx.soil.nblayers_soil) && (in > 0.0))
            {
                if (in > (ctx.soil.FC_SWC[l] - ctx.soil.SWC3D[l][d]))
                {
                    in -= (ctx.soil.FC_SWC[l] - ctx.soil.SWC3D[l][d]);
                    ctx.soil.SWC3D[l][d] = ctx.soil.FC_SWC[l];
                    if (isnan(ctx.soil.SWC3D[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) <= 0)
                    {
                        cout << "incorrect ctx.soil.SWC3D, Min/ctx.soil.Max_SWC" << endl;
                        cout << ctx.soil.Max_SWC[l] << endl;
                    }
                }
                else
                {
                    ctx.soil.SWC3D[l][d] += in;
                    if (isnan(ctx.soil.SWC3D[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0)
                    {
                        cout << "incorrect ctx.soil.SWC3D, Min/ctx.soil.Max_SWC" << endl;
                        cout << ctx.soil.Throughfall[d] << "\t" << in << "\t" << ctx.climate.precip << "\t" << ctx.soil.Interception[d] << "\t" << ctx.soil.LAI_DCELL[0][d] << endl;
                    }
                    in = 0.0;
                }
                l++;
            }
        }
        else
        { // if the top soil ctx.diag.layer is already saturated (eg. inundated forest), throughfall -> runoff
            ctx.soil.Runoff[d] = ctx.soil.Throughfall[d];
        }
        // ctx.soil.Leakage
        ctx.soil.Leakage[d] = in;
    }
    // END of the BUCKET MODEL.

    // Update of soil water potential field
    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            // ctx.soil.soil_phi3D[l][d]=ctx.soil.phi_e[l]*pow((ctx.soil.SWC3D[l][d]/ctx.soil.Max_SWC[l]), -ctx.soil.b[l]);

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
                ctx.soil.Ks[l][d] = ctx.soil.Ksat[l] * pow(theta_w, 0.5) * inter * inter; // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)

                if (isnan(ctx.soil.soil_phi3D[l][d]) || isnan(ctx.soil.Ks[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0) //|| ctx.soil.KsPhi[l][d]==0.0 || ctx.soil.Ks[l][d]==0.0 || ctx.soil.soil_phi3D[l][d]==0.0)
                    cout << "In bucket model, ctx.diag.layer " << l << " dcell " << d << " theta_w=" << theta_w << " ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) << " ctx.soil.soil_phi3D[l][d]=" << ctx.soil.soil_phi3D[l][d] << " ctx.soil.Ksat=" << ctx.soil.Ksat[l] << " ctx.soil.Ks[l][d]=" << ctx.soil.Ks[l][d] << endl;
            }
            else if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {
                ctx.soil.soil_phi3D[l][d] = ctx.soil.phi_e[l] * pow(theta_w, -ctx.soil.b[l]);           // this is the soil water characteristic of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                ctx.soil.Ks[l][d] = ctx.soil.Ksat[l] * pow(theta_w, 2.5 + 2 * ctx.soil.b[l]);           // this is the hydraulic conductivity curve of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                ctx.soil.KsPhi[l][d] = ctx.soil.Ksat[l] * ctx.soil.phi_e[l] * pow(theta_w, 2.5 + ctx.soil.b[l]); // ctx.soil.Ks times ctx.soil.soil_phi3D, computed directly as the exact power of ctx.params.theta.

                if (isnan(ctx.soil.soil_phi3D[l][d]) || isnan(ctx.soil.Ks[l][d]) || isnan(ctx.soil.KsPhi[l][d]) || (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) < 0) //|| ctx.soil.KsPhi[l][d]==0.0 || ctx.soil.Ks[l][d]==0.0 || ctx.soil.soil_phi3D[l][d]==0.0)
                    cout << "In bucket model, ctx.diag.layer " << l << " dcell " << d << " theta_w=" << theta_w << " ctx.soil.SWC3D[l][d]-ctx.soil.Min_SWC[l]=" << (ctx.soil.SWC3D[l][d] - ctx.soil.Min_SWC[l]) << " ctx.soil.soil_phi3D[l][d]=" << ctx.soil.soil_phi3D[l][d] << " ctx.soil.Ksat=" << ctx.soil.Ksat[l] << " ctx.soil.phi_e=" << ctx.soil.phi_e[l] << " ctx.soil.b[l]=" << ctx.soil.b[l] << " ctx.soil.KsPhi[l][d]=" << ctx.soil.KsPhi[l][d] << " ctx.soil.Ks[l][d]=" << ctx.soil.Ks[l][d] << endl;
                // KsPhi2[l][d]=ctx.soil.Ksat[l]*ctx.soil.phi_e[l]*pow(theta_w, 2.5);
                //  we may want to shift to the van Genuchten-Mualem expressions of ctx.soil.soil_phi3D and ctx.soil.Ks, as the van genuchten-Mualem model is currently defacto the more standard soil hydraulic model (see ref in Table 1 in Marthews et al. 2014). To do so, see if we have data of soil pH, cation exchange capacity, organic carbon content, to explicitly compute the parameters with Hodnett & Tomasella 2002 (as recommended by Marthews et al. 2014 -- Table 2; or instead directly use the parameter provided by the map in Marthews et al. 2014.
            }
        }
    }
#endif
}

// #############################
//  Global function: update ctx.species.SPECIES_SEEDS field
// #############################
void FillSeed(Context &ctx, int col, int row, int spp)
{
    if ((col >= 0) && (col < ctx.grid.cols))
    {
        if ((row >= 0) && (row < ctx.grid.rows))
        {
            int site = col + ctx.grid.cols * row;
            // if(ctx.T[site].t_age == 0){
            if (ctx.opt._SEEDTRADEOFF)
                ctx.species.SPECIES_SEEDS[site][spp]++; // ifdef SEEDTRADEOFF, ctx.species.SPECIES_SEEDS[site][spp] is the number of seeds of this species at that site
            else
                ctx.species.SPECIES_SEEDS[site][spp] = 1; // If s_Seed[site] = 0, site is not occupied, if s_Seed[site] > 1, s_Seed[site] is the presence of a seed
            // cout << "site: " << site << " spp: " << spp << " Seed added!!! " << endl;
            // }
        }
    }
}

// #############################
//  Global function: tree germination module
// #############################
void RecruitTree(Context &ctx)
{
    for (int site = 0; site < ctx.grid.sites; site++)
    { //**** Local germination ****
        if (ctx.T[site].t_age == 0)
        {
            int spp_withseeds = 0;
            for (int spp = 1; spp <= ctx.grid.nbspp; spp++)
            { // lists all the species with a seed present at given site...
                if (ctx.species.SPECIES_SEEDS[site][spp] > 0)
                {
                    // write species that are present to an extra array
                    ctx.species.SPECIES_GERM[spp_withseeds] = spp;
                    spp_withseeds++;
                }
            }
            if (spp_withseeds > 0)
            { // ... and then randomly select one of these species

                // new in v.2.4.1: for consistency use genrand2() instead of rand(), since v.2.5: use gsl RNG
                int spp_index = int(gsl_rng_uniform_int(ctx.rng.gslrand, spp_withseeds));
                int spp = ctx.species.SPECIES_GERM[spp_index];
                // otherwise all species with seeds present are equiprobable

#ifdef LCP_alternative

#ifdef WATER
                if (ctx.soil.soil_phi3D[0][ctx.grid.site_DCELL[site]] > 0.5 * ctx.S[spp].s_tlp)
                {
                    ctx.T[site].Birth(ctx, spp, site); // in this version, the light environment is checked within Birth() function
                }

#else
                ctx.T[site].Birth(ctx, spp, site); // in this version, the light environment is checked within Birth() function

#endif

#else // LCP_alternative

                float flux = ctx.climate.WDailyMean * exp(-fmaxf(ctx.field.LAI3D[0][site + ctx.grid.SBORD], 0.0) * ctx.params.kpar);
#ifdef WATER
                if (flux > (ctx.S[spp].s_LCP) && ctx.soil.soil_phi3D[0][ctx.grid.site_DCELL[site]] > 0.5 * ctx.S[spp].s_tlp)
                {
                    ctx.T[site].Birth(ctx, spp, site);
                    ctx.species.SPECIES_SEEDS[site][spp] = 0; // newIM nov2021, to adjust to the yearly update of Species_seeds
                }
                // in addition to a condition of light availability- hence light demanding species may not be able to grow in deep sahde conditions in understorey -, a condition on water availability is added - hence drought-intolerant species may not be recruited in water-stressd conditions
#else
                // If enough light, germination, initialization of NPP (LCP is the species light compensation point
                // here, light is the sole environmental resources tested as a limiting factor for germination, but we should think about adding nutrients (N,P) and water conditions...
                if (flux > (ctx.S[spp].s_LCP))
                {
                    ctx.T[site].Birth(ctx, spp, site);
                    ctx.species.SPECIES_SEEDS[site][spp] = 0; // newIM nov2021, to adjust to the yearly update of Species_seeds
                }
#endif // WATER
#endif // LCP_alternative
            }
        }
    }
}

// #############################
//  Global function: Treefall gap formation
// #############################
//! change in v.2.4: resetting ctx.field.Thurt[0] field is done in TriggerSecondaryTreefall() at the beginning of each iteration. Further changes: rewriting of Tree::FallTree() which is now Tree::Treefall(angle). t_hurt can now persist longer, so new treefall events are added to older damages (that, in turn are decaying)
void TriggerTreefall(Context &ctx)
{
    for (int site = 0; site < ctx.grid.sites; site++)
        if (ctx.T[site].t_age)
        {
            // treefall is triggered given a certain flexural force
            // ctx.opt._BASICTREEFALL: just dependent on height threshold + random uniform distribution
            float angle = 0.0, c_forceflex = 0.0;
            if (ctx.opt._BASICTREEFALL)
            {
                c_forceflex = (1 - (1 - gsl_rng_uniform(ctx.rng.gslrand)) / (12 * ctx.time.timestep)) * ctx.T[site].t_height; // probability of treefall per month = 1-t_Ct/t_height , compare to genrand2(), if ctx.time.timestep=1/12: genrand2() < 1 - t_Ct/t_height, or: genrand2() > t_Ct/t_height
                angle = float(twoPi * gsl_rng_uniform(ctx.rng.gslrand));                                         // random angle
            }
            // above a given stress threshold the tree falls
            if (c_forceflex > ctx.T[site].t_Ct)
            {
                ctx.T[site].Treefall(ctx, angle);
            }
        }
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        // Update of Field hurt
        if (ctx.T[site].t_age)
        {
            ctx.T[site].t_hurt = max(ctx.field.Thurt[0][site + ctx.grid.sites], ctx.T[site].t_hurt); // NEW in v.2.4: addition of damages, alternative: max()
        }
    }
}

// #############################
//  Global function: Secondary treefall gap formation
// #############################
//! - NEW in v.2.4: TriggerSecondaryTreefall(), called at the beginning of each iteration
//! - translates damages from previous round into tree deaths, partly treefalls, partly removing them only (e.g. splintering)
//! - in the limit of ctx.params.p_tfsecondary = 0.0, this is equivalent to the previous computation
void TriggerTreefallSecondary(Context &ctx)
{
    ctx.diag.nbTreefall1 = 0;
    ctx.diag.nbTreefall10 = 0;
    ctx.diag.nbTreefall30 = 0;
#ifdef Output_ABC
    nbTreefall10_abc = 0;
#endif
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        ctx.field.Thurt[0][site] = ctx.field.Thurt[0][site + 2 * ctx.grid.sites] = 0;
        ctx.field.Thurt[0][site + ctx.grid.sites] = 0;
    }
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        if (ctx.T[site].t_age)
        {
            float height_threshold = ctx.T[site].t_height / ctx.T[site].t_mult_height; // since 2.5: a tree's stability is defined by its species' average height, i.e. we divide by the intraspecific height multiplier to account for lower stability in quickly growing trees; otherwise slender, faster growing trees would be treated preferentially and experience less secondary treefall than more heavily built trees
            if (2.0 * ctx.T[site].t_hurt * (1 - (1 - gsl_rng_uniform(ctx.rng.gslrand)) / (12 * ctx.time.timestep)) > height_threshold)
            { // check whether tree dies: probability of death per month is 1.0-0.5*t_height/t_hurt, so, when ctx.time.timestep=1/12, ctx.rng.gslrand <= 1.0 - 0.5 * t_height/t_hurt, or ctx.rng.gslrand > 0.5 * t_height/t_hurt; modified in v.2.5: probability of death is 1.0 - 0.5*t_height/(t_mult_height * t_hurt), so the larger the height deviation (more slender), the higher the risk of being thrown by another tree
                if (ctx.params.p_tfsecondary > gsl_rng_uniform(ctx.rng.gslrand))
                {                                                          // check whether tree falls or dies otherwise
                    float angle = float(twoPi * gsl_rng_uniform(ctx.rng.gslrand)); // random angle
                    ctx.T[site].Treefall(ctx, angle);
                }
                else
                {
                    ctx.T[site].Death(ctx);
                }
            }
            else
            {
                ctx.T[site].t_hurt = short(ctx.params.hurt_decay * float(ctx.T[site].t_hurt)); // reduction of t_hurt according to ctx.params.hurt_decay, could be moved to Tree::Growth() function and made dependent on the tree's carbon gain
            }
        }
    }

}

// Helper function
int GetTimeofyear(Context &ctx)
{
    // new function to derive time of year, extended to negative iterations (-1 would be treated as last iteration of previous year)
    int timeofyear;
    if (ctx.time.iter < 0)
        timeofyear = ctx.time.iterperyear - abs(ctx.time.iter) % ctx.time.iterperyear;
    else
        timeofyear = ctx.time.iter % ctx.time.iterperyear;
    return (timeofyear);
};
float CalcHeightBaseline(float &ah, float &hmax, float &dbh)
{
    // height allometry
    float height = hmax * dbh / (dbh + ah);
    return (height);
}

float CalcCRBaseline(Context &ctx, float &dbh)
{
    // crown radius allometry
    float CR;
    if (!ctx.opt._CROWN_MM)
    {
        CR = exp(ctx.params.CR_a + ctx.params.CR_b * log(dbh)); // power law, the default
    }
    else
        CR = ctx.params.CR_b * dbh / (dbh + ctx.params.CR_a); // Michaelis Menten type allometry !!!: requires ctx.params.CR_b to be the CR_max parameter and ctx.params.CR_a the initial increase */
    // for reference, two crown allometries that are reasonable in French Guiana
    // t_CR = t_mult_CR * exp(1.9472 + 0.5925*log(t_dbh)); // crown allometry deduced from Piste Saint-Elie */
    // t_CR = t_mult_CR * exp(1.8814 + 0.5869*log(t_dbh)); // this is crown allometry derived from data set compiled by Jucker et al. 2016 (Global Change Biology)
    return (CR);
}

float CalcCDBaseline(Context &ctx, float &height)
{
    // crown depth allometry
    // since v.2.5, simplification of the computation of the crown depth, in accordance with the Canopy Constructor algorithm
    float CD = (ctx.params.CD_a + ctx.params.CD_b * height);
    return (CD);
}

//! - upper bound on LAI within one voxel and above voxel (beyond 9.95 and 19.95, none of the environmental variables should change), needed for LookUp tables
int CalcIntabsorb(float absorb_prev, float absorb_delta)
{
    absorb_delta = fminf(absorb_delta, 9.95);
    absorb_prev = fminf(absorb_prev, 19.95);
    int intabsorb = int(absorb_prev * 20.0) + 400 * int(absorb_delta * 20.0);
    return (intabsorb);
}
// Helper function
//! - upper bound on LAI within one voxel and above voxel, needed for LookUp tables, for cases where there won't be any absorption within the voxel (absorb_delta = 0.0)
int CalcIntabsorb(float absorb_prev)
{
    absorb_prev = fminf(absorb_prev, 19.95);
    int intabsorb = int(absorb_prev * 20.0);
    return (intabsorb);
}
