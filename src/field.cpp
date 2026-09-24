#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

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
        int ha = ctx.grid.sites / 10000;
        gsl_ran_multinomial(ctx.rng.gslrand, ctx.grid.sites, ctx.params.Cseedrain * ha, ctx.species.p_seed, ctx.species.n_seed);
        cout << ctx.grid.sites << " Seedrain: " << ctx.params.Cseedrain * ha << endl;
        int seedsadded = 0;
        for (int s = 0; s < ctx.grid.sites; s++)
        {
            // if(ctx.T[s].t_age == 0){
            int nbseeds = ctx.species.n_seed[s];
            // cout << "Site: " << s << " nbseeds: " << nbseeds << " ctx.grid.nbspp: " << ctx.grid.nbspp << endl;
            gsl_ran_multinomial(ctx.rng.gslrand, ctx.grid.nbspp, nbseeds, ctx.species.p_species, ctx.species.n_species);
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
#ifdef MPI
    // Reinitialize field ctx.field.LAI3D
    for (int i = 0; i < 2; i++)
        for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
            for (int site = 0; site < 2 * ctx.grid.SBORD; site++)
                LAIc[i][h][site] = 0;
#endif

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

#ifdef MPI
    // Communicate border of field
    // MPI_ShareField(ctx.field.LAI3D,LAIc,2*ctx.grid.SBORD);
    This MPI command no longer exists in openMPI
            Action 20 /
        01 / 2016 TODO : FIX THIS MPI_ShareField(ctx.field.LAI3D, LAIc, 2 * ctx.grid.SBORD);
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        //! Add border effects in local fields
        if (mpi_rank)
            for (site = 0; site < 2 * ctx.grid.SBORD; site++)
                ctx.field.LAI3D[h][site] += LAIc[0][h][site];
        if (mpi_rank < mpi_size - 1)
            for (int site = 0; site < 2 * ctx.grid.SBORD; site++)
                ctx.field.LAI3D[h][site + ctx.grid.sites] += LAIc[1][h][site];
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

// Global function: a modifying function that converts LAI to the density of a specific ctx.diag.layer, using the GetDensity functions
//! - modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results
//! - LAI is the input, dens_layer the output
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop)
{
    int crown_top = int(height);
    int crown_base = int(height - CD);
    float dens_top, dens_belowtop, dens_lowerlayer; // RENAMED: dens → dens_lowerlayer (base-layer crown density from GetDensitiesGradient; distinct from output param dens_layer)
#ifdef LAI_gradient
    GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens_lowerlayer);
#else
    GetDensityUniform(LAI, CD, dens_lowerlayer);
    dens_top = dens_belowtop = dens_lowerlayer;
#endif

    if (CD < 3.0 && crown_top == crown_base)
    {
        dens_layer = LAI; /* full LAI allocation */
    }
    else if (CD < 3.0 && (crown_top - layer_fromtop == crown_base))
    {
        dens_layer = LAI;
    }
    else
    {
        float fraction_layer = height - floor(height); /* this is the fraction that each ctx.diag.layer apart from the topmost ctx.diag.layer will extend into the voxel above */
        if (layer_fromtop == 0)
            dens_layer = dens_top * fraction_layer;
        else if (layer_fromtop == 1)
            dens_layer = dens_top + dens_belowtop * fraction_layer;
        else if (layer_fromtop == 2)
            dens_layer = dens_top + dens_belowtop + dens_lowerlayer * fraction_layer;
        else
            dens_layer = LAI;
    }
}

// Global function: a modifying function that converts LAI to the density of a specific ctx.diag.layer, using the GetDensity functions
//! - modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results
//! - LAI is the input, dens_layer the output
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop)
{

    int crown_top = int(height);
    int crown_base = int(height - CD);
    float dens_top, dens_belowtop, dens_lowerlayer; // RENAMED: dens → dens_lowerlayer (base-layer crown density from GetDensitiesGradient; distinct from output param dens_layer)
#ifdef LAI_gradient
    GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens_lowerlayer);
#else
    GetDensityUniform(LAI, CD, dens_lowerlayer);
    dens_top = dens_belowtop = dens_lowerlayer;
#endif

    if (CD < 3.0 && crown_top == crown_base)
    {
        dens_layer = dens_top * CD;
    }
    else if (CD < 3.0 && (crown_top - layer_fromtop == crown_base))
    {
        float fraction_belowbase = float(crown_base + 1) - (height - CD);
        dens_layer = dens_lowerlayer * fraction_belowbase;
    }
    else
    {
        float fraction_layer = height - floor(height);         // this is the fraction that each ctx.diag.layer apart from the topmost ctx.diag.layer will extend into the voxel above
        float fraction_layer_fromabove = 1.0 - fraction_layer; // the inverse of the fraction above

        if (layer_fromtop == 0)
            dens_layer = dens_top * fraction_layer;
        else if (layer_fromtop == 1)
            dens_layer = dens_top * fraction_layer_fromabove + dens_belowtop * fraction_layer;
        else if (layer_fromtop == 2)
            dens_layer = dens_belowtop * fraction_layer_fromabove + dens_lowerlayer * fraction_layer;
        else
            dens_layer = dens_lowerlayer * fraction_layer_fromabove;
    }
}

// Global function: update of ctx.field.LAI3D field, called by CalcLAI()
void UpdateLAI3D(Context &ctx, int height, int site, float dens_layer, float &LA_cumulated) // RENAMED: ctx.params.dens → dens_layer
{
    ctx.field.LAI3D[height][site + ctx.grid.SBORD] += dens_layer;
    LA_cumulated += dens_layer;
}

// Global function: remove outliers in canopy height model (CHM); vector option
#ifdef CHM_SPIKEFREE
void UpdateCHMvector(int height, int site, float noinput, vector<int> &chm)
{
    if (chm[site] < height)
        chm[site] = height;
}
// Global function: remove outliers in canopy height model (CHM)
void UpdateCHM(int height, int site, float noinput, int *chm)
{
    if (chm[site] < height)
        chm[site] = height;
}
#endif

// Global function: PPFD retrieval for function leafarea_max()
void GetPPFDabove(Context &ctx, int height, int site, float noinput, float (&ppfd_CA)[2])
{
    // First get voxel field densities
    float absorb_prev = ctx.field.LAI3D[height + 1][site + ctx.grid.SBORD];
    int intabsorb = CalcIntabsorb(absorb_prev);

    // Obtain PPFD for the voxel, and also record the circled area
    ppfd_CA[0] += ctx.climate.WDailyMean * ctx.lookup.LookUp_flux[intabsorb];
    ppfd_CA[1] += 1.0; // add area
}

// Global function: calculates the canopy environment
//! - this function adds to the environmental variables provided in canopy_environment_cumulated
//! - the PPFD, VPD, Tmp and leafarea_layer retrieval function for Fluxh()
#ifdef WATER
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[6]) // RENAMED: ctx.params.dens → dens_layer
{
#else
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[4]) // RENAMED: ctx.params.dens → dens_layer
{
#endif
    // first get voxel field densities
    float absorb_prev = ctx.field.LAI3D[height + 1][site + ctx.grid.SBORD];
    float absorb_curr = ctx.field.LAI3D[height][site + ctx.grid.SBORD];
    float absorb_delta = absorb_curr - absorb_prev;
    if (absorb_delta < 0.0)
        absorb_delta = 0.0; // eliminate rounding errors
    int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);
#ifdef WATER
    int intincident = CalcIntabsorb(absorb_prev);
#endif

#ifdef FULL_CLIMATE
    // Obtain PPFD, VPD and T for the voxel
    float PPFD_voxel = ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD_voxel = ctx.lookup.LookUp_VPD[intabsorb];
    float T_voxel = ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float PPFD_voxel_incident = ctx.lookup.LookUp_flux[intincident];
    float ExtinctLW_voxel = ctx.lookup.LookUp_ExtinctLW[intincident];
#endif // WATER

#else // FULL_CLIMATE
      // Obtain PPFD, VPD and T for the voxel
    float PPFD_voxel = ctx.climate.WDailyMean * ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD_voxel = ctx.climate.VPDDailyMean * ctx.lookup.LookUp_VPD[intabsorb];
    float T_voxel = ctx.climate.tDailyMean - ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float PPFD_voxel_incident = ctx.climate.WDailyMean * ctx.lookup.LookUp_flux[intincident];
    float ExtinctLW_voxel = ctx.lookup.LookUp_ExtinctLW[intincident];
#endif // WATER

#endif // FULL_CLIMATE

    // Add the three variables up, weighted by leaf density inside voxel
    canopy_environment_cumulated[0] += dens_layer;
    canopy_environment_cumulated[1] += PPFD_voxel * dens_layer;
    canopy_environment_cumulated[2] += VPD_voxel * dens_layer;
    canopy_environment_cumulated[3] += T_voxel * dens_layer;
#ifdef WATER
    canopy_environment_cumulated[4] += PPFD_voxel_incident * dens_layer;
    canopy_environment_cumulated[5] += ExtinctLW_voxel * dens_layer;
#endif
}

// Global function: calculates packing densities
void AddCrownVolumeLayer(Context &ctx, int row_center, int col_center, float height, float CR, float CD, int crownvolume[70])
{
    int crown_top = int(height);
    int crown_base = int(height - CD);
    if (CD <= 3.0)
    {
        // for the smallest crowns it is simply cylinder rings being filled up
        int crown_intarea = GetCrownIntarea(CR);
        for (int h = crown_top; h >= crown_base; h--)
            crownvolume[h] += crown_intarea;
    }
    else
    {
        // For the rest of the crown, we go through different crown shells. We separate out the innermost sector (a slowly increasing cylinder), and the surrounding parts of the crown
        // first the metrics with respect to the internal crown structure (i.e. z coordinate with respect to crown base)
        float crownshell_base = height - CD + 2.0;                   // lower reference point for the crown slope function is two layers up from the crown base
        float crownshell_extent = height - crownshell_base;          // this is the extent from the "base ctx.diag.layer" to the top
        float crownshell_extent_toplayer = floor(crownshell_extent); // this is the extent to the lower limit of the toplayer
        // then we translate the crown coordinates into discretised variables with respect to the absolute location in the voxel field, as needed for location in the voxel field, with layers defined from top to bottom
        int shell_fromtop = 0;
        int height_innermost = crown_top - shell_fromtop;
        int height_toplayer = int(crownshell_base + crownshell_extent_toplayer) - shell_fromtop;
        int height_baselayer = int(crownshell_base + 1.0) - shell_fromtop;

        // now do calculations
        // first the inner crown shell section that grows dynamically
        float radius_innermost = GetRadiusSlope(ctx, CR, crownshell_extent, crownshell_extent_toplayer);
        int crown_intarea_innermost = GetCrownIntarea(radius_innermost);
        for (int h = height_innermost; h >= crown_base; h--)
        {
            crownvolume[h] += crown_intarea_innermost;
        }
        // now loop through the outer crown shell cylinders
        for (int h_outer = height_toplayer; h_outer >= crown_base; h_outer--)
        {
            // calculating the radius of the current ctx.diag.layer depending on the respective slopes, to be replaced by function
            // float radius_height = CR - crown_slope * (h_outer - height_baselayer);    // for the lowest ctx.diag.layer, i.e. h == height_baselayer, radius = t_CR
            int extent_layerouter = max(h_outer - height_baselayer, 0); // we also fill up underneath the baselayer
            float radius_height = GetRadiusSlope(ctx, CR, crownshell_extent, extent_layerouter);
            int crown_intarea = GetCrownIntarea(radius_height);

            crownvolume[h_outer] += (crown_intarea - crown_intarea_innermost);
        }
    }
}

// ###############################################
//  Update of the ctx.field.LAI3D field
//! called by UpdateField
// #################################################
//! - modified in v.2.3: additional contribution to voxels that are not fully occupied by the tree crown. !!!: this does not calculate ctx.field.LAI3D directly, this only calculates the density in each voxel belonging to a tree. The final LAI field is calculated outside of the class Tree
//! - modified in v.2.4 and v.2.5: introducing an alternative crown shape, "umbrella"-like, inspired by previous shell models and similar to the crown shapes in the PPA. If activated, crowns contain three layers of vegetation that, once the crown goes beyond 3m in depth, will bend downwards on the edges with a linear slope. Since v.2.5 all loops (CalcLAI, Fluxh, leafarea_max) are executed through the same template. This allows to implement other crown shapes in the future and ensures that modifications are carried through across the code
#ifdef CROWN_UMBRELLA
void CalcLAI(Context &ctx, Tree &tree)
{
    if (tree.t_age > 0)
    {
        int site_crowncenter = tree.t_site + tree.t_CrownDisplacement;
        int row_crowncenter = site_crowncenter / ctx.grid.cols;
        int col_crowncenter = site_crowncenter % ctx.grid.cols;

        float LA_cumulated = 0.0; // Currently, an output variable is required by LoopLayerUpdateCrownStatistic_template, we here use LA_cumulated as control variable

        int crown_top = int(tree.t_height);
        int crown_base = int(tree.t_height - tree.t_CD);
        int max_shells = min(crown_top - crown_base + 1, 4); // since the new crown shapes

        for (int shell_fromtop = 0; shell_fromtop < max_shells; shell_fromtop++)
        {
            LoopLayerUpdateCrownStatistic_template(
                ctx,
                row_crowncenter, col_crowncenter,
                tree.t_height, tree.t_CR, tree.t_CD,
                tree.t_fraction_filled,
                shell_fromtop,
                [&ctx](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); },
                tree.t_LAI,
                LA_cumulated,
                LAI2dens,
                [&ctx](int h, int s, float d, float &la){ UpdateLAI3D(ctx, h, s, d, la); });
        }
    }
}
#else
void CalcLAI(Context &ctx, Tree &tree)
{
    if (tree.t_age > 0)
    {
        int crown_base = int(tree.t_height - tree.t_CD),
            crown_top = int(tree.t_height);
        int site_crowncenter = tree.t_site + tree.t_CrownDisplacement;
        int row_crowncenter = site_crowncenter / ctx.grid.cols;
        int col_crowncenter = site_crowncenter % ctx.grid.cols;

        float fraction_abovetop = tree.t_height - float(crown_top);
        float fraction_belowbase = float(crown_base + 1) - (tree.t_height - tree.t_CD);

        float crown_area = PI * tree.t_CR * tree.t_CR;      // floor of crown_area to bound area accumulation
        int crown_intarea = int(crown_area);      // floor of crown_area to bound area accumulation
        crown_intarea = max(crown_intarea, 1);    // minimum area of crown (1)
        crown_intarea = min(crown_intarea, 1963); // maximum area of crown (radius 25), int(3.14*25*25), int(3.14*25*25)

        float dens_avg = tree.t_LAI / tree.t_CD;
        float fraction_filled_target = tree.t_fraction_filled;

        for (int h = crown_base; h <= crown_top; h++)
        {
            float fraction_filled_actual = 0.0;
            float dens_layer;

            dens_layer = dens_avg;
            if (crown_top == crown_base)
                dens_layer *= tree.t_CD;
            else if (h == crown_top)
                dens_layer *= fraction_abovetop;
            else if (h == crown_base)
                dens_layer *= fraction_belowbase;

            for (int i = 0; i < crown_intarea; i++)
            {
                if (fraction_filled_actual > fraction_filled_target)
                {
                    fraction_filled_actual = (fraction_filled_actual * float(i)) / (float(i) + 1.0);
                }
                else
                {
                    fraction_filled_actual = (fraction_filled_actual * float(i) + 1.0) / (float(i) + 1.0);
                    int site_relative = ctx.lookup.LookUp_Crown_site[i];
                    int row = row_crowncenter + site_relative / 51 - 25;
                    int col = col_crowncenter + site_relative % 51 - 25;

                    if (row >= 0 && row < ctx.grid.rows && col >= 0 && col < ctx.grid.cols)
                    {
                        int site = col + ctx.grid.cols * row + ctx.grid.SBORD;
                        ctx.field.LAI3D[h][site] += dens_layer;
                    }
                }
            }
        }
    }
}
#endif

#ifdef MPI
// MPI Routines

// ##############################################
//  Global MPI function: Communication of border fields in the parallel version of the code
// ##############################################
//! - Only if the MPI option has been enabled
void MPI_ShareSeed(unsigned char **c, int n)
{

    MPI_Status status;

    if (p_rank == size - 1)
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_CHAR, size - 2, 0, c[3], n, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    if (p_rank == 0)
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_CHAR, size - 1, 0, c[3], n, MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD, &status);
    if ((p_rank) && (p_rank < size - 1))
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_CHAR, p_rank - 1, 0, c[3], n, MPI_UNSIGNED_CHAR, p_rank + 1, 0, MPI_COMM_WORLD, &status);

    if (p_rank == 0)
        MPI_Sendrecv(c[1], n, MPI_UNSIGNED_CHAR, 1, 1, c[2], n, MPI_UNSIGNED_CHAR, size - 1, 1, MPI_COMM_WORLD, &status);
    if (p_rank == size - 1)
        MPI_Sendrecv(c[1], n, MPI_UNSIGNED_CHAR, 0, 1, c[2], n, MPI_UNSIGNED_CHAR, size - 2, 1, MPI_COMM_WORLD, &status);
    if ((p_rank) && (p_rank < size - 1))
        MPI_Sendrecv(c[1], n, MPI_UNSIGNED_CHAR, p_rank + 1, 1, c[2], n, MPI_UNSIGNED_CHAR, p_rank - 1, 1, MPI_COMM_WORLD, &status);
}

// ##############################################
//  Global MPI function: Communication of fields
// ##############################################
void MPI_ShareField(Context &ctx, unsigned short **cl, unsigned short ***cp, int n)
{

    MPI_Status status;
    for (int h = 0; h < (ctx.grid.HEIGHT + 1); h++)
    {
        if (p_rank == 0)
            MPI_Sendrecv(cl[h], n, MPI_UNSIGNED_SHORT, size - 1, h, cp[1][h], n, MPI_UNSIGNED_SHORT, 1, h, MPI_COMM_WORLD, &status);
        if (p_rank == size - 1)
            MPI_Sendrecv(cl[h], n, MPI_UNSIGNED_SHORT, size - 2, h, cp[1][h], n, MPI_UNSIGNED_SHORT, 0, h, MPI_COMM_WORLD, &status);
        if ((p_rank) && (p_rank < size - 1))
            MPI_Sendrecv(cl[h], n, MPI_UNSIGNED_SHORT, p_rank - 1, h, cp[1][h], n, MPI_UNSIGNED_SHORT, p_rank + 1, h, MPI_COMM_WORLD, &status);

        if (p_rank == 0)
            MPI_Sendrecv(cl[h] + ctx.grid.sites, n, MPI_UNSIGNED_SHORT, 1, h + ctx.grid.HEIGHT, cp[0][h], n, MPI_UNSIGNED_SHORT, size - 1, h + ctx.grid.HEIGHT, MPI_COMM_WORLD, &status);
        if (p_rank == size - 1)
            MPI_Sendrecv(cl[h] + ctx.grid.sites, n, MPI_UNSIGNED_SHORT, 0, h + ctx.grid.HEIGHT, cp[0][h], n, MPI_UNSIGNED_SHORT, size - 2, h + ctx.grid.HEIGHT, MPI_COMM_WORLD, &status);
        if ((p_rank) && (p_rank < size - 1))
            MPI_Sendrecv(cl[h] + ctx.grid.sites, n, MPI_UNSIGNED_SHORT, p_rank + 1, h + ctx.grid.HEIGHT, cp[0][h], n, MPI_UNSIGNED_SHORT, p_rank - 1, h + ctx.grid.HEIGHT, MPI_COMM_WORLD, &status);
    }
}

// ##############################################
//  Global MPI function: Communication of treefalls
// ##############################################
void MPI_ShareTreefall(unsigned short **c, int n)
{

    MPI_Status status;
    if (p_rank == 0)
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_SHORT, size - 1, 0, c[2], n, MPI_UNSIGNED_SHORT, 1, 0, MPI_COMM_WORLD, &status);
    if (p_rank == size - 1)
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_SHORT, size - 2, 0, c[2], n, MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);
    if ((p_rank) && (p_rank < size - 1))
        MPI_Sendrecv(c[0], n, MPI_UNSIGNED_SHORT, p_rank - 1, 0, c[2], n, MPI_UNSIGNED_SHORT, p_rank + 1, 0, MPI_COMM_WORLD, &status);

    if (p_rank == 0)
        MPI_Sendrecv(c[0] + 2 * n, n, MPI_UNSIGNED_SHORT, 1, 1, c[1], n, MPI_UNSIGNED_SHORT, size - 1, 1, MPI_COMM_WORLD, &status);
    if (p_rank == size - 1)
        MPI_Sendrecv(c[0] + 2 * n, n, MPI_UNSIGNED_SHORT, 0, 1, c[1], n, MPI_UNSIGNED_SHORT, size - 2, 1, MPI_COMM_WORLD, &status);
    if ((p_rank) && (p_rank < size - 1))
        MPI_Sendrecv(c[0] + 2 * n, n, MPI_UNSIGNED_SHORT, p_rank + 1, 1, c[1], n, MPI_UNSIGNED_SHORT, p_rank - 1, 1, MPI_COMM_WORLD, &status);
}
#endif

