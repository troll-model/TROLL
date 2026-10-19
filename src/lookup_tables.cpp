#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

void InitialiseIntraspecific(Context &ctx)
{
    float max_intraspecific_height = 0.0, min_intraspecific_height = 1000.0,
          max_intraspecific_CR = 0.0, min_intraspecific_CR = 1000.0,
          max_intraspecific_CD = 0.0, min_intraspecific_CD = 1000.0,
          max_intraspecific_P = 0.0, min_intraspecific_P = 1000.0,
          max_intraspecific_N = 0.0, min_intraspecific_N = 1000.0,
          max_intraspecific_LMA = 0.0, min_intraspecific_LMA = 1000.0,
          max_intraspecific_wsg = 0.0, min_intraspecific_wsg = 1000.0,
          max_intraspecific_dbhmax = 0.0, min_intraspecific_dbhmax = 1000.0;
#ifdef WATER
    float max_intraspecific_leafarea = 0.0, min_intraspecific_leafarea = 1000.0,
          max_intraspecific_tlp = 0.0, min_intraspecific_tlp = 1000.0;
    double variation_leafarea = 0.0, variation_tlp = 0.0;
#endif

    double variation_height = 0.0, variation_CR = 0.0, variation_CD = 0.0, variation_P = 0.0, variation_N = 0.0, variation_LMA = 0.0, variation_wsg = 0.0, variation_dbhmax = 0.0;
    for (int i = 0; i < 10000; i++)
    { // modified FF v.3.1.5 (reduced from 100000 to 10000)
        if (ctx.rng.covariance_status == 0)
        {
            variation_N = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_N);
            variation_P = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_P);
            variation_LMA = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_LMA);
        }
        else
        {
            gsl_ran_multivariate_gaussian(ctx.rng.gslrand, ctx.rng.mu_N_P_LMA, ctx.rng.mcov_N_P_LMA, ctx.rng.variation_N_P_LMA);
            variation_N = gsl_vector_get(ctx.rng.variation_N_P_LMA, 0);
            variation_P = gsl_vector_get(ctx.rng.variation_N_P_LMA, 1);
            variation_LMA = gsl_vector_get(ctx.rng.variation_N_P_LMA, 2);
        }
        gsl_ran_bivariate_gaussian(ctx.rng.gslrand, ctx.intra.sigma_height, ctx.intra.sigma_CR, ctx.intra.corr_CR_height, &variation_height, &variation_CR);
        // variation_height = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_height);
        // variation_CR = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_CR);
        variation_CD = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_CD);
        variation_wsg = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_wsg);
        variation_dbhmax = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_dbhmax);
#ifdef WATER
        variation_leafarea = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_leafarea);
        variation_tlp = gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_tlp);
#endif

        // limit extent of variation, some coarse biological limits for variation around allometries
        if (variation_N > 0.5)
            variation_N = 0.5;
        if (variation_N < -0.5)
            variation_N = -0.5;
        if (variation_P > 0.5)
            variation_P = 0.5;
        if (variation_P < -0.5)
            variation_P = -0.5;
        if (variation_LMA > 0.5)
            variation_LMA = 0.5;
        if (variation_LMA < -0.5)
            variation_LMA = -0.5;
        if (variation_CD > 0.5)
            variation_CD = 0.5;
        if (variation_CD < -0.5)
            variation_CD = -0.5;
        if (variation_CR > 0.5)
            variation_CR = 0.5;
        if (variation_CR < -0.5)
            variation_CR = -0.5;
        if (variation_height > 0.5)
            variation_height = 0.5;
        if (variation_height < -0.5)
            variation_height = -0.5;
#ifdef WATER
        if (variation_leafarea > 0.5)
            variation_leafarea = 0.5;
        if (variation_leafarea < -0.5)
            variation_leafarea = -0.5;
        if (variation_tlp > 0.5)
            variation_tlp = 0.5;
        if (variation_tlp < -0.5)
            variation_tlp = -0.5;
        ctx.intra.d_intraspecific_leafarea[i] = float(exp(variation_leafarea));
        ctx.intra.d_intraspecific_tlp[i] = float(exp(variation_tlp));
#endif

        ctx.intra.d_intraspecific_height[i] = float(exp(variation_height));
        ctx.intra.d_intraspecific_CR[i] = float(exp(variation_CR));
        ctx.intra.d_intraspecific_N[i] = float(exp(variation_N));
        ctx.intra.d_intraspecific_P[i] = float(exp(variation_P));
        ctx.intra.d_intraspecific_LMA[i] = float(exp(variation_LMA));
        ctx.intra.d_intraspecific_CD[i] = float(exp(variation_CD));
        ctx.intra.d_intraspecific_wsg[i] = float(variation_wsg); // normal, not log-normal
        ctx.intra.d_intraspecific_dbhmax[i] = float(exp(variation_dbhmax));
        // ctx.intra.d_intraspecific_height[i] = exp(float(gsl_ran_gaussian(ctx.rng.gslrand, ctx.intra.sigma_height)));
        max_intraspecific_height = fmaxf(max_intraspecific_height, ctx.intra.d_intraspecific_height[i]);
        min_intraspecific_height = fminf(min_intraspecific_height, ctx.intra.d_intraspecific_height[i]);
        max_intraspecific_CR = fmaxf(max_intraspecific_CR, ctx.intra.d_intraspecific_CR[i]);
        min_intraspecific_CR = fminf(min_intraspecific_CR, ctx.intra.d_intraspecific_CR[i]);
        max_intraspecific_N = fmaxf(max_intraspecific_N, ctx.intra.d_intraspecific_N[i]);
        min_intraspecific_N = fminf(min_intraspecific_N, ctx.intra.d_intraspecific_N[i]);
        max_intraspecific_P = fmaxf(max_intraspecific_P, ctx.intra.d_intraspecific_P[i]);
        min_intraspecific_P = fminf(min_intraspecific_P, ctx.intra.d_intraspecific_P[i]);
        max_intraspecific_LMA = fmaxf(max_intraspecific_LMA, ctx.intra.d_intraspecific_LMA[i]);
        min_intraspecific_LMA = fminf(min_intraspecific_LMA, ctx.intra.d_intraspecific_LMA[i]);
        max_intraspecific_CD = fmaxf(max_intraspecific_CD, ctx.intra.d_intraspecific_CD[i]);
        min_intraspecific_CD = fminf(min_intraspecific_CD, ctx.intra.d_intraspecific_CD[i]);
        max_intraspecific_wsg = fmaxf(max_intraspecific_wsg, ctx.intra.d_intraspecific_wsg[i]);
        min_intraspecific_wsg = fminf(min_intraspecific_wsg, ctx.intra.d_intraspecific_wsg[i]);
        max_intraspecific_dbhmax = fmaxf(max_intraspecific_dbhmax, ctx.intra.d_intraspecific_dbhmax[i]);
        min_intraspecific_dbhmax = fminf(min_intraspecific_dbhmax, ctx.intra.d_intraspecific_dbhmax[i]);
#ifdef WATER
        max_intraspecific_leafarea = fmaxf(max_intraspecific_leafarea, ctx.intra.d_intraspecific_leafarea[i]);
        min_intraspecific_leafarea = fminf(min_intraspecific_leafarea, ctx.intra.d_intraspecific_leafarea[i]);
        max_intraspecific_tlp = fmaxf(max_intraspecific_tlp, ctx.intra.d_intraspecific_tlp[i]);
        min_intraspecific_tlp = fminf(min_intraspecific_tlp, ctx.intra.d_intraspecific_tlp[i]);
#endif
    }
    cout << endl
         << "Intraspecific variation initialisation: " << endl;
    cout << "Max and min allometry deviation, lognormal (height): " << max_intraspecific_height << " | " << min_intraspecific_height << endl;
    cout << "Max and min allometry deviation, lognormal (crown radius): " << max_intraspecific_CR << " | " << min_intraspecific_CR << endl;
    cout << "Max and min trait deviation, lognormal (N): " << max_intraspecific_N << " | " << min_intraspecific_N << endl;
    cout << "Max and min trait deviation, lognormal (P): " << max_intraspecific_P << " | " << min_intraspecific_P << endl;
    cout << "Max and min trait deviation, lognormal (LMA): " << max_intraspecific_LMA << " | " << min_intraspecific_LMA << endl;
    cout << "Max and min allometry deviation, normal (crown depth): " << max_intraspecific_CD << " | " << min_intraspecific_CD << endl;
    cout << "Max and min trait deviation, normal (wsg): " << max_intraspecific_wsg << " | " << min_intraspecific_wsg << endl;
    cout << "Max and min trait deviation, lognormal (dmax): " << max_intraspecific_dbhmax << " | " << min_intraspecific_dbhmax << endl;
}

#ifdef LCP_alternative
// v.3.1.5: create a LookUp table based on species identity and intraspecific deviation
//! LAImax of one pseudo-tree (species spp, intraspecific variant dev) under the current
//! ctx conditions; callers make sure these are the initialisation conditions.
static float ComputeLAImaxEntry(Context &ctx, int spp, int dev)
{
    Tree pseudotree(ctx);
    pseudotree.t_sp_lab = spp;

    pseudotree.t_Pmass = ctx.S[spp].s_Pmass * ctx.intra.d_intraspecific_P[dev];
    pseudotree.t_Nmass = ctx.S[spp].s_Nmass * ctx.intra.d_intraspecific_N[dev];
    pseudotree.t_LMA = ctx.S[spp].s_LMA * ctx.intra.d_intraspecific_LMA[dev];

#ifdef WATER
    pseudotree.t_leafarea = ctx.S[spp].s_leafarea * ctx.intra.d_intraspecific_leafarea[dev];
    pseudotree.t_wsg = fmaxf(ctx.S[spp].s_wsg + ctx.intra.d_intraspecific_wsg[dev], 0.05);
    pseudotree.t_wleaf = sqrt(pseudotree.t_leafarea * 0.0001);
    pseudotree.t_WSF = 1;
    pseudotree.t_WSF_A = 1;
    pseudotree.t_g1 = (-3.97 * pseudotree.t_wsg + 6.53); // this is the relationship provided by Lin et al. 2015
#endif

    pseudotree.t_Vcmax = CalcVcmaxm(pseudotree.t_LMA, pseudotree.t_Nmass, pseudotree.t_Pmass) * pseudotree.t_LMA;
    pseudotree.t_Jmax = pseudotree.CalcJmaxm(ctx) * pseudotree.t_LMA;
    pseudotree.t_Rdark = CalcRdark(pseudotree.t_LMA, pseudotree.t_Nmass, pseudotree.t_Pmass, pseudotree.t_Vcmax);

    pseudotree.CalcLAImax(ctx); // reads ctx only: safe to call from several threads
    return pseudotree.t_LAImax;
}

//! LAImax for species spp (1-based) and intraspecific variant dev, computed on first use.
//! Computing all 10000 * nbspp entries upfront cost ~80 us each (about 11 s for 14
//! species, minutes for species-rich parameterisations), while short runs use few of them.
float GetLAImaxPrecomputed(Context &ctx, int spp, int dev)
{
    const int index = dev + (spp - 1) * 10000;
    if (!ctx.lookup.LookUpLAImax_ready[index])
    {
        // The table is defined for the conditions at initialisation (iteration -1, which makes
        // dailyFluxesLeaf use yearly mean climate, and the initial night temperature). An entry
        // computed later during the run temporarily restores them, so it is bit-identical to
        // the former eager computation.
        const int iter_now = ctx.time.iter;
        const float tnight_now = ctx.climate.tnight;
        ctx.time.iter = ctx.lookup.LookUpLAImax_iter;
        ctx.climate.tnight = ctx.lookup.LookUpLAImax_tnight;
        ctx.lookup.LookUpLAImax[index] = ComputeLAImaxEntry(ctx, spp, dev);
        ctx.time.iter = iter_now;
        ctx.climate.tnight = tnight_now;
        ctx.lookup.LookUpLAImax_ready[index] = 1;
    }
    return ctx.lookup.LookUpLAImax[index];
}

//! Set up the LAImax table (10000 intraspecific variants per species).
//! With several OpenMP threads, all entries are computed now, in parallel: they are
//! independent and use no random numbers, so results do not depend on the thread count.
//! With a single thread, entries are computed on first use by GetLAImaxPrecomputed(),
//! which saves most of the start-up time of short runs.
void InitialiseLookUpLAImax(Context &ctx)
{
    const int nb_entries = 10000 * ctx.grid.nbspp;
    ctx.lookup.LookUpLAImax.assign(nb_entries, 0.0f);
    ctx.lookup.LookUpLAImax_ready.assign(nb_entries, 0);
    ctx.lookup.LookUpLAImax_iter = ctx.time.iter;
    ctx.lookup.LookUpLAImax_tnight = ctx.climate.tnight;

    if (TrollThreads() <= 1)
    {
        cout << "LookUp table for LAImax: entries computed on first use (single thread)." << endl;
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (int index = 0; index < nb_entries; index++)
    {
        ctx.lookup.LookUpLAImax[index] = ComputeLAImaxEntry(ctx, index / 10000 + 1, index % 10000);
        ctx.lookup.LookUpLAImax_ready[index] = 1;
    }

    float minLAImax = 10.0; // fixed upper limit of LAImax is currently 10.0 m2/m2
    float maxLAImax = 0.0;
    float avgLAImax = 0.0;
    for (int index = 0; index < nb_entries; index++)
    {
        minLAImax = fminf(minLAImax, ctx.lookup.LookUpLAImax[index]);
        maxLAImax = fmaxf(maxLAImax, ctx.lookup.LookUpLAImax[index]);
        avgLAImax += ctx.lookup.LookUpLAImax[index];
    }
    avgLAImax *= 1.0 / float(nb_entries);
    cout << "Calculated LookUp table for LAImax on " << TrollThreads() << " threads. Min LAImax is: " << minLAImax << " | max LAImax is: " << maxLAImax << " avg LAImax is: " << avgLAImax << endl;
}
#endif

//! Global function: initialise lookup tables
//! - this contains a number of empirical functions including:
//! -# temperature dependence on Farquhar model parameters (following von Caemmerer 2000 and Bernacchi et al. 2003 PCE)
//! -# flux averaging with the canopy
//! -# ctx.grid.sites within a crown in order of distance from the center
void InitialiseLookUpTables(Context &ctx)
{

    ctx.lookup.nbTbins = 500;
    float Taccuracy = 0.1;
    ctx.lookup.iTaccuracy = 1.0 / Taccuracy;
    cout << endl
         << "Built-in maximal temperature: " << float(ctx.lookup.nbTbins) * Taccuracy << endl;
    ctx.lookup.LookUp_KmT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_GammaT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_VcmaxT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_JmaxT.assign(ctx.lookup.nbTbins, float());
    // if(NULL==(LookUp_Rday=new float[ctx.lookup.nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rday" << endl; //newIM: no redundancy anymore between LookUp_Rday and LookUP_Rnight
    ctx.lookup.LookUp_Rleaf.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_Rstem.assign(ctx.lookup.nbTbins, float());
    // if(NULL==(LookUp_Rnight=new float[ctx.lookup.nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rnight" << endl; //newIM: no redundancy anymore between LookUp_Rday and LookUP_Rnight
#ifdef WATER
    ctx.lookup.nbVPDbins = 600;
    float VPDaccuracy = 0.01;
    ctx.lookup.iVPDaccuracy = 1.0 / VPDaccuracy;
    cout << "Built-in maximal VPD: " << float(ctx.lookup.nbVPDbins) * VPDaccuracy << endl;
    ctx.lookup.LookUp_INLR.resize(ctx.lookup.nbTbins);
    for (int i = 0; i < ctx.lookup.nbTbins; i++) // row 0 included: it was never allocated (wild pointer), which broke the LUT cache
    {
        ctx.lookup.LookUp_INLR[i].assign(ctx.lookup.nbVPDbins, float());
    }
    // if(NULL==(LookUp_ESAT=new float[ctx.lookup.nbTbins])) cerr<<"!!! Mem_Alloc LookUp_ESAT" << endl;
    ctx.lookup.LookUp_SLOPE.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_GRADN.assign(ctx.lookup.nbTbins, float());
    float ESAT0 = 0.61121;
#endif
    for (int i = 1; i < ctx.lookup.nbTbins; i++)
    { // loop over "T" in GPPleaf()
        float temper = float(i) * Taccuracy;
        // !!!UPDATE provide references for these equations

#ifdef WATER
        ctx.lookup.LookUp_KmT[i] = 404.0 * exp(((temper - 25.0) / (298 * 0.00831 * (273 + temper))) * 59.36) *
                        (1 + 210 * 1.0 / 248.0 * exp(-(temper - 25.0) / (298 * 0.00831 * (273 + temper)) * 35.94)); // taken from von Caemmerer 2000, as in Domingues et al. 2010, for consistency
        ctx.lookup.LookUp_GammaT[i] = 37.0 * exp(((temper - 25.0) / (298 * 0.00831 * (273 + temper))) * 23.4);                 // taken from von Caemmerer 2000, as in Domingues et al. 2010, for consistency
#else
        ctx.lookup.LookUp_KmT[i] = 404.0 * exp(((temper - 25.0) / (298 * 0.00831 * (273 + temper))) * 59.36) *
                        (1 + 210 * 1.0 / 248.0 * exp(-(temper - 25.0) / (298 * 0.00831 * (273 + temper)) * 35.94)) * ctx.params.iCair; // taken from von Caemmerer 2000, as in Domingues et al. 2010, for consistency
        ctx.lookup.LookUp_GammaT[i] = 37.0 * exp(((temper - 25.0) / (298 * 0.00831 * (273 + temper))) * 23.4) * ctx.params.iCair;                 // taken from von Caemmerer 2000, as in Domingues et al. 2010, for consistency
#endif
        ctx.lookup.LookUp_VcmaxT[i] = exp(26.35 - 65.33 / (0.00831 * (temper + 273.15))); // taken from Bernacchi et al. 2003 PCE, as in Domingues et al. 2010 for consistency
        ctx.lookup.LookUp_JmaxT[i] = exp(17.57 - 43.54 / (0.00831 * (temper + 273.15)));  // taken from Bernacchi et al. 2003 PCE, as in Domingues et al. 2010 for consistency
        // LookUp_Rday[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper))); //newIM: no redundancy anymore between LookUp_Rday and LookUP_Rnight
        ctx.lookup.LookUp_Rleaf[i] = exp((temper - 25.0) * 0.1 * log(3.09 - 0.0215 * (25.0 + temper))); // this is equ. 1 in Atkin et al. 2015 New phytologist //newIM: no redundancy anymore between LookUp_Rday and LookUP_Rnight
        ctx.lookup.LookUp_Rstem[i] = 39.6 * 378.7 * ctx.time.timestep * exp(((temper - 25.0) / 10.0) * log(2.0));
        // LookUp_Rnight[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper))); //newIM: no redundancy anymore between LookUp_Rday and LookUP_Rnight
        //  exp((temp-25)/10*log(2)) is the temperature dependency of Rstem, supposing a constant Q10=2, according to Ryan et al 1994 and Meir & Grace 2002 exp((ctx.climate.tnight-25)*0.1*log(3.09-0.0215*(25+ctx.climate.tnight))) is the temperature dependencies used by Atkin 2015 (equ1)
#ifdef WATER
        // for computation of isothermal net long-range radiation
        float ESAT = 0.61121 * exp((18.678 - temper / 234.5) * temper / (257.14 + temper)); // Saturation partial pressure of water vapour in kPa From Jones 2013, Eq (5.15) page 102 (in agreement with Cochard 2019 equ 2 -- Buck equation)
        // float INLR=SIGMA*EMLEAF*pow(temper-ABSZERO,4.0); // Longwave radiation Jones (2014) Eq. 5.5, p101. // as in Jerome's script, but does not provide sound values (pb with EMLEAF ==> to be checked in Jones)
        for (int v = 1; v < ctx.lookup.nbVPDbins; v++)
        {
            float vpd = float(v) * VPDaccuracy;
            // float ea=1000.0*(ESAT - vpd);
            // float emaB0=0.642*pow(ea/(temper-ABSZERO),0.1428571); // cf. equ. D4 in Leuning et al. 1995 (and in Duursma's R package code)
            //  float emaB=1.24*pow(ea/(temper-ABSZERO),0.1428571); // cf. Brutsaert 1975 equ 11 (see also Marthews et al. 2012 Theor Appl Climatol, Table 2). Weird that this expression provides values >1 !
            // float emaI=1-0.261*exp(-0.000777*temper*temper); // cf. Idso & Jackson 1969 (see also Marthews et al. 2012 Theor Appl Climatol, Table 2)
            float ema = 1 / (SIGMA * pow(temper - ABSZERO, 4.0)) * (59.38 + (113.7 * pow((ABSZERO - temper) / ABSZERO, 6)) + 96.96 * sqrt(186 * vpd / (temper - ABSZERO))); // cf. Dilley & O'Brien 1998 (see also Marthews et al. 2012 Theor Appl Climatol, Table 2); 186=4650/25
            float INLR = (1 - ema) * SIGMA * pow(temper - ABSZERO, 4.0);                                                                                                    // as in Duursma R package, following Leuning et al. 1995 PCE, equs D2-D5
            ctx.lookup.LookUp_INLR[i][v] = INLR;
            // cout << "T=" << temper << " vpd=" << vpd <<  " ema=" << ema <<  " INLR=" << INLR << endl;
        }
        float SLOPE = 1000.0 * (ESAT - ESAT0) / Taccuracy; // Slope of the ESAT function, but in Pa (not in kPa)
        ESAT0 = ESAT;
        float GRADN = 4. * EMLEAF * SIGMA * pow(temper - ABSZERO, 3.0) / (CPAIR * AIRMA); // Equivalent conductance to radiative heat transfer Jones 2013 (page 101)
        // LookUp_ESAT[i]=ESAT;

        ctx.lookup.LookUp_SLOPE[i] = SLOPE;
        ctx.lookup.LookUp_GRADN[i] = GRADN;
#endif
    }

#ifdef WATER

    ctx.lookup.nbHbins = 2000;
    float Haccuracy = 0.01;
    ctx.lookup.iHaccuracy = 1.0 / Haccuracy;
    ctx.lookup.LookUp_Wind.assign(ctx.lookup.nbHbins, float());

    for (int i = 1; i < ctx.lookup.nbHbins; i++)
    { // loop over "T" in GPPleaf()
        float h_ratio = float(i) * Haccuracy;

        if (h_ratio >= 1)
        {
            // cppcheck-suppress invalidFunctionArg ; false positive: h_ratio >= 1 here, so the argument is >= 3.33
            ctx.lookup.LookUp_Wind[i] = 0.831 * log(16.67 * ((h_ratio)-0.8)); // drawn from Monteith & Unsworth 2008, with d = 0.8 H and z0 = 0.06 H (cf. Rau et al's TROLL wind manuscript).
        }
        else
        {
            ctx.lookup.LookUp_Wind[i] = exp(-alphaInoue * (1 - h_ratio));
        }
    }

#endif

    // look up table for flux averaging/integration
    // division into absorption prior to current voxel (absorb_prev) and absorption in current voxel (absorb_delta)
    // prior absorption has a maximum of 20 m2/m3, while absorption within one voxel can maximally reach 10 m2/m3
    ctx.lookup.LookUp_flux_absorption.assign(80000, float());
    ctx.lookup.LookUp_flux.assign(80000, float());
    ctx.lookup.LookUp_ExtinctLW.assign(80000, float());
    ctx.lookup.LookUp_VPD.assign(80000, float());
    ctx.lookup.LookUp_T.assign(80000, float());
    for (int i = 0; i < 400; i++)
    { // loop over "absorb" in Fluxh()
        for (int j = 0; j < 200; j++)
        {
            float absorb_prev = float(i) / 20.0;
            float absorb_delta = float(j) / 20.0;
            if (absorb_delta == 0)
            {
                // flux is now computed simply as absorption on a per m2 plant matter basis (cf. explanation below) but since there is no plant matter in the case of absorb_delta == 0, PPFD should be zero here
                ctx.lookup.LookUp_flux_absorption[i + 400 * j] = 0.0; // in this case

                // if the voxel does not contain any plant matter, values are constant across the voxel, e.g. just the top value calculated from absorb_prev
                ctx.lookup.LookUp_flux[i + 400 * j] = exp(-ctx.params.kpar * absorb_prev);
                ctx.lookup.LookUp_ExtinctLW[i + 400 * j] = 0.8 * exp(-0.8 * absorb_prev); // extinction of thermal radiation, following Leuning et al. 1995 equ. D1
                ctx.lookup.LookUp_VPD[i + 400 * j] = 0.25 + sqrt(fmaxf(0.0, 0.08035714 * (7.0 - absorb_prev)));
                ctx.lookup.LookUp_T[i + 400 * j] = 0.4285714 * (fminf(7.0, absorb_prev));
            }
            else
            {
                // Flux is now computed simply as absorption on a per m2 plant matter basis, and not incident flux, since the FvCB model requires a transformation of incident into absorbed flux (cf. original Farquhar 1980 paper, or Medlyn et al. 2002, Plant, Cell & Environment). For example, Medlyn et al. use 0.093 quantum yield, 4 mol electron/mol photon and an absorptance of leaves of 0.8 to arrive at a factor of 0.3 to be multiplied with incident PPFD. Now, in a dense forest, ctx.params.kpar modifies the absorbed PPFD per leaf area. For simplicity, we separate the 0.093 and 4 mol electrion/mol photon from the absorptance and calculate the latter directly from an effective ctx.params.kpar which includes a general k (i.e. leaf angle distribution) and a leaf absorptance factor, e.g. 0.9. The main effect of this scheme is that lowering k (reflecting, for example, steeper leaf angles) may result in more incident light per m2 ground, but absorption also gets lower, since leaves are not perfectly illuminated
                // To calculate absorbed PPFD, a formula can either be derived through integration or be motivated as follows: 1) incoming flux is exp(-ctx.params.kpar * absorb_prev), equivalent to what was previously computed as flux, 2) the absorbed fraction of the incoming flux in a ctx.diag.layer of "absorb_delta" is (1.0 - exp(-ctx.params.kpar * absorb_delta), and 3) the amount of leaf area per ground area is absorb_delta, which is needed as divisor to convert to absorption per m2 leaf area
                // Since PPFD is a density (i.e. given relative to m2 leaf area), a lower absorb_delta results in slightly higher PPFD. A lower absorb_delta implies that leaves are less densely distributed in space, so there is a slight increase in absorbed photons per leaf area, even though overall absorbed photon numbers decrease. While an absorb_delta = 0 implies zero absorption, in the limit of very low absorb_delta (-> 0), the absorption approaches ctx.params.kpar (Taylor expansion: exp(x) ~ 1 + x, so (1 - exp(-ctx.params.kpar*x))/x ~ ctx.params.kpar*x/x ~ ctx.params.kpar). This is not realistic, since leaves cannot get infinitesimally small and the assumptions of Beer-Lambert breaks down beforehand. But since the linear approximation should be justified in low density layers, maybe this could be used to accelerate the computation?
                ctx.lookup.LookUp_flux_absorption[i + 400 * j] = exp(-ctx.params.kpar * absorb_prev) * (1.0 - exp(-ctx.params.kpar * absorb_delta)) / absorb_delta;

                // an alternative to calculating the absorbed flux density, is to calculate the average flux density
                // for voxels of 1 unit length depth, this corresponds just to the integral over LAI, which can be decomposed into a constant absorb_prev and a linearly increasing absorb_delta
                // once LAI reaches the critical value of 7.0, VPD And T do not decrease anymore, hence the distinction between two cases

                ctx.lookup.LookUp_flux[i + 400 * j] = exp(-ctx.params.kpar * absorb_prev) * (1.0 - exp(-ctx.params.kpar * absorb_delta)) / (ctx.params.kpar * absorb_delta);
                ctx.lookup.LookUp_ExtinctLW[i + 400 * j] = 0.8 * exp(-0.8 * absorb_prev); // extinction of thermal radiation, following Leuning et al. 1995 equ. D1

                if (absorb_prev + absorb_delta >= 7)
                { // This is related to the empirical description of VPD and temperature decrease through the canopy. These empirical description were derived from literature and data from HOBO on the COPAS tower, and assumed that above a LAI of 7, the VPD and temperature are pretty constant (e.g. Within the understory).
                    ctx.lookup.LookUp_VPD[i + 400 * j] = 0.25;
                    ctx.lookup.LookUp_T[i + 400 * j] = 3.0; // 0.4285714 * 7.0
                }
                else
                {
                    ctx.lookup.LookUp_VPD[i + 400 * j] = 0.25 + (0.188982 / absorb_delta) * (pow((7.0 - absorb_prev), 1.5) - pow((7.0 - absorb_prev - absorb_delta), 1.5));
                    ctx.lookup.LookUp_T[i + 400 * j] = 0.4285714 * (absorb_prev + 0.5 * absorb_delta);
                }
            }
        }
    }

    // new in v.2.4: LookUp table that gives the ctx.grid.sites within a crown in order of distance from the center
    // crowns can thus be assembled from inside out, in a radial fashion
    int Crown_dist[2601]; // this saves the distances from the center of the crown
    int extent = 25;      // maximum extent of crowns (after test: either enlarge or allocate dynamically)
    int extent_full = 2 * 25 + 1;
    int index_crown = 0, xx, yy, site_rel, dist;
    Crown_dist[index_crown] = 0;                                    // this is the distance of the center of the crown from the center (0)
    ctx.lookup.LookUp_Crown_site[index_crown] = extent + extent_full * extent; // this is the label of the site at the center of the crown ( x = extent, y = extent)

    // loop over crown
    for (int col = 0; col < extent_full; col++)
    {
        for (int row = 0; row < extent_full; row++)
        {
            xx = col - extent; // distance from center (x = extent) in x direction
            yy = row - extent; // distance from center (y = extent) in y direction
            if (!((xx == 0) && (yy == 0)))
            {
                site_rel = col + extent_full * row;
                dist = xx * xx + yy * yy;
                // now order the arrays according to distance from center
                // index_crown saves last filled position in array
                // for every voxel we run through the array from position zero to last filled position, and check where to put the new value
                for (int i = 0; i < index_crown + 1; i++)
                {
                    int temp = Crown_dist[i];
                    int site_rel_temp = ctx.lookup.LookUp_Crown_site[i];
                    if (dist <= temp)
                    { // if distance is smaller than at current array position, move everything up
                        Crown_dist[i] = dist;
                        ctx.lookup.LookUp_Crown_site[i] = site_rel;
                        dist = temp;
                        site_rel = site_rel_temp;
                    }
                }
                Crown_dist[index_crown + 1] = dist; // the last value that has been pushed into the dist variable fills a new position
                index_crown = index_crown + 1;
            }
        }
    }

    // new in v.3.1.6: reshuffle the allocation patterns slightly to avoid too regular-looking crowns (a cosmetic change)
    // in previous versions, a small minority of trees had fraction_filled values of very close to 0.75, 0.5, 0.25. Since we do not randomize allocation of filled/empty pixels per crown, but simply increment it step-wise (i.e. one empty gap in the crown every 4 pixels for a filled fraction of 0.75, or one empty gap every 2 pixels for 0.5, etc.), this meant that the determination of "holes" in the crown could accidentally align with the allocation of crown pixels in general (which usually have 4 or 8 pixels of the same distance following each other), and this could lead to weirdly regular shapes. This usually only occurred for a tiny fraction of trees and did not influence model behavior, but may create odd-looking outputs for people unfamiliar with models/TROLL, so it is desirable to remove it.
    // in the new version, we simply shuffle a few values that are close to each other - they have a similar or the same distance, so the overall outward-spiralling allocation remains the same, just with a slightly shifted order. To break the 4-based pattern, we use an increment of 7 and a distance of 3 for the exchanges. We start at the 6th pixel as the order will not matter for the innermost crown part.
    for (int i = 5; i < 2598; i += 7)
    {
        int site_current = ctx.lookup.LookUp_Crown_site[i];
        int site_exchange = ctx.lookup.LookUp_Crown_site[i + 3];
        ctx.lookup.LookUp_Crown_site[i] = site_exchange;
        ctx.lookup.LookUp_Crown_site[i + 3] = site_current;
    }

#ifdef LCP_alternative
    InitialiseLookUpLAImax(ctx);
#endif
}

