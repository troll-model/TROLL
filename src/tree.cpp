#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"


// #############################################
//  Tree constructor
// #############################################
Tree::Tree(Context &ctx)
{
    t_from_Data = 0;
    t_sp_lab = 0;
    t_age = 0;
    t_hurt = 0;
    t_NPP = t_GPP = t_Rday = t_Rnight = t_Rstem = 0.0; // new v.2.2
    t_dbh = t_height = t_CR = t_CD = 0.0;
    t_CrownDisplacement = 0;
#ifdef MIP_Lichstein
    t_inInventory = 0;
#endif

    if (ctx.opt._NDD)
    {
        t_NDDfield.reserve(ctx.grid.nbspp + 1);
        for (int sp = 0; sp < (ctx.grid.nbspp + 1); sp++)
            t_NDDfield.push_back(0.0);
    }
#if defined(WATER)
    t_transpiration = 0.0;
    t_root_biomass.reserve(ctx.soil.nblayers_soil);
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        t_root_biomass.push_back(0.0);
    }
    t_soil_layer_weight.reserve(ctx.soil.nblayers_soil);
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        t_soil_layer_weight.push_back(0.0);
    }
#endif

#ifdef Output_ABC
    t_dbh_previous = 0.0;
#endif
    if (ctx.opt._BASICTREEFALL)
        t_Ct = 0.0;
};

// #############################################
//  Tree birth
// #############################################
//! Actions needed at birth for a tree (precisely, at time when a stem enters the > 1 cm trunk diameter class)
void Tree::Birth(Context &ctx, int nume, int site0)
{

    int dev_rand = int(gsl_rng_uniform_int(ctx.rng.gslrand, 10000)); // modified FF v.3.1.5 (reduced to 10000)

#ifdef LCP_alternative

    // ######################
    // # first test LAImax ##
    // ######################
    int index_LAImax = dev_rand + (nume - 1) * 10000;
    float LAImax_precomputed = ctx.lookup.LookUpLAImax[index_LAImax];

    if (ctx.field.LAI3D[0][site0 + ctx.grid.SBORD] < LAImax_precomputed)
    {

#endif

        ctx.diag.nblivetrees++;
        t_site = site0;
        t_sp_lab = nume;
        ctx.S[t_sp_lab].s_nbind++;
        t_multiplier_seed = 1;
        t_age = 1;
        t_hurt = 0;
        t_NPPneg = 0;

#ifdef MIP_Lichstein
        t_inInventory = 0;
#endif

        // ######################
        // # intraspecific var ##
        // ######################

        t_mult_height = ctx.intra.d_intraspecific_height[dev_rand];
        t_mult_CR = ctx.intra.d_intraspecific_CR[dev_rand];
        t_mult_N = ctx.intra.d_intraspecific_N[dev_rand];
        t_mult_P = ctx.intra.d_intraspecific_P[dev_rand];
        t_mult_LMA = ctx.intra.d_intraspecific_LMA[dev_rand];
        t_mult_CD = ctx.intra.d_intraspecific_CD[dev_rand];
        t_dev_wsg = ctx.intra.d_intraspecific_wsg[dev_rand];
        t_mult_dbhmax = ctx.intra.d_intraspecific_dbhmax[dev_rand];
#ifdef WATER
        t_mult_leafarea = ctx.intra.d_intraspecific_leafarea[dev_rand];
        t_mult_tlp = ctx.intra.d_intraspecific_tlp[dev_rand];
#endif

        // #####################
        // # leaf/wood traits ##
        // #####################
        t_Pmass = ctx.S[t_sp_lab].s_Pmass * t_mult_P;
        t_Nmass = ctx.S[t_sp_lab].s_Nmass * t_mult_N;
        t_LMA = ctx.S[t_sp_lab].s_LMA * t_mult_LMA;
        t_wsg = fmaxf(ctx.S[t_sp_lab].s_wsg + t_dev_wsg, 0.05); // updated in v.3.1: Lower cutoff for normal distribution set to 0.05 instead of 0.1 (some wsg measurements in Global Wood Density Database go below 0.1)
#ifdef WATER
        t_leafarea = ctx.S[t_sp_lab].s_leafarea * t_mult_leafarea;
        t_tlp = ctx.S[t_sp_lab].s_tlp * t_mult_tlp;
        t_wleaf = sqrt(t_leafarea * 0.0001);
        t_phi_lethal = -0.9842 + 3.1795 * t_tlp; // Inferred from data provided in Bartlett et al. 2016 PNAS
        t_itlp = 1 / t_tlp;
        t_b = -2.23 * t_itlp; // this results from the following assumption: (i) ctx.params.phi(TLP)=0.97*ctx.params.phi(gs90), from Martin-StPaul et al. 2017 Ecology letters; (ii) WSF=0.1 at ctx.params.phi(gs90)
#endif

        t_Vcmax = CalcVcmaxm(t_LMA, t_Nmass, t_Pmass) * t_LMA; // in micromolC ctx.params.m-2 s-1
        t_Jmax = CalcJmaxm(ctx) * t_LMA;                          // in micromolC ctx.params.m-2 s-1
        t_Rdark = CalcRdark(t_LMA, t_Nmass, t_Pmass, t_Vcmax); // in micromolC ctx.params.m-2 s-1
#ifdef WATER
        // t_g1_0=-0.0224 * t_LMA + 4.8278; // this is the relationship provided by Wu et al. 2019 GCB
        t_g1_0 = (-3.97 * t_wsg + 6.53); // this is the relationship provided by Lin et al. 2015
        t_g1 = t_g1_0;
        if (t_g1 < 0)
        {
            cout << "Warning: t_g1 <0 in Tree::Birth ! t_LMA=" << t_LMA << "; s_LMA=" << ctx.S[t_sp_lab].s_LMA << "; multiplierLMA=" << t_mult_LMA << endl;
        }
#endif

        // #################$####
        // # biometry/allometry ##
        // ##################$###

        t_hmax = ctx.S[t_sp_lab].s_hmax;
        t_ah = ctx.S[t_sp_lab].s_ah;

        t_dbh = ctx.params.DBH0;
        t_dbhmax = ctx.S[t_sp_lab].s_dbhmax;
        t_dbhmax *= t_mult_dbhmax;
        t_dbhmax = fmaxf(t_dbhmax, ctx.params.DBH0 * 2.0);
        t_dbhmature = t_dbhmax * 0.5; // Mean threshold of tree size to maturity - see Visser et al. 2016 Functional Ecology (suited to both understory and top-canopy species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230

        UpdateHeight(ctx);
        UpdateCR(ctx);
        UpdateCD(ctx);

        t_CrownDisplacement = 0;

        if (ctx.opt._BASICTREEFALL)
        {
            t_Ct = CalcCt(ctx);
        }

        // ##################
        // # leaf dynamics ##
        // ##################

#ifdef WATER
        t_WSF = t_WSF_A = 1.0; // Tree leaf lifespan and LAImax are computed under no water stress and with average climatic conditions (radiation, temperature).
                               // This could/should be updated later, so that these two quantities change under water stress, and also seasonally.
                               // t_WSF and t_WSF_A are then updated later in Tree::Birth to account for the real water conditions at birth.
#endif
        float crown_area = PI * t_CR * t_CR;
        float fraction_filled_general = 1.0 - ctx.crown.crown_gap_fraction;
        t_fraction_filled = fminf(fraction_filled_general / (t_mult_CR * t_mult_CR), 1.0);
        float crown_area_nogaps = GetCrownAreaFilled(crown_area);

        if (ctx.opt._LA_regulation > 0)
        {

            // In this case, we determine both maximum LAI (theoretical if fully exposed to sunlight) and maximum leaf area (actual, i.e. also considering shading by neighboring trees)
            //  v.3.1.5, we can simply set t_LAImax to the precomputed value
            // CalcLAImax(ctx);
            // cout << "LAImax: " << t_LAImax << " precomputed: " << LAImax_precomputed << endl; // check consistency
#ifdef LCP_alternative
            t_LAImax = LAImax_precomputed;
#else
        CalcLAImax(ctx);
#endif
            float LAIexperienced_eff;
            CalcLAmax(ctx, LAIexperienced_eff, t_LAmax);
            // t_LA = 0.5 * t_LAmax;    //tets post_fusion 27/01/2023 IM
            t_LA = 0.25 * t_LAmax; // Initially trees are set to have a quarter of their maximum leaf area

            t_LAI = t_LA / crown_area_nogaps;
            t_carbon_storage = CalcCarbonStorageMax(ctx) * 0.5; // Initial value with half of the maximum storage
            t_carbon_biometry = 0.0;
        }
        else
        {
#ifdef CROWN_UMBRELLA
            t_LAI = ctx.params.dens * fminf(t_CD, 3.0);
#else
        t_LAI = ctx.params.dens * t_CD;
#endif
            t_LA = t_LAI * crown_area_nogaps;
        }
        CalcLeafLifespan(ctx);
        InitialiseLeafPools(ctx); // new in v.3.1, formerly part of CalcLeafLifespan()

#ifdef PHENO_DROUGHT
        t_Ndays_dry = 0;
        t_Ndays_wet = 0;
        t_pheno_factor = 1.0;
#endif

        // ############
        // # sapwood ##
        // ############
        float ddbh = t_dbh;   // Initially each stem can be treated as juvenile and stem cross section is mostly sapwood
        t_sapwood_area = 0.0; // the sapwood area is initialized to zero and then updated
        UpdateSapwoodArea(ctx, ddbh);

#ifdef WATER
        // ##########
        // # water ##
        // ##########
        // Water_availability(ctx);
        t_transpiration = 0.0;
        Water_availability(ctx); // Roots are not set here, but at the beginning of Tree::Update (however see comments within function Tree::Water_availability)
                              // UpdateRootDistribution();
#endif

        // ###############
        // # for output ##
        // ###############
#ifdef Output_ABC
        t_dbh_previous = t_dbh;
#endif

#ifdef TRACK_INDIVIDUALS
        StartTracking(ctx);
#endif

#ifdef LCP_alternative
        ctx.species.SPECIES_SEEDS[site0][nume] = 0; // newIM nov2021, to adjust to the yearly update of Species_seeds
    }
#endif
}

// #############################################
//  Tree initialisation from field data
// ##############################################
//! - Modelled for compatibility with Tree::Birth. For comments cf. Tree::Birth
//! - For comments regarding allometries and t_LA cf. Tree::Growth.
int Tree::BirthFromInventory(Context &ctx, int site, vector<string> &parameter_names, vector<string> &parameter_values, int &nb_speciesrandom)
{
    int success = 0;
    bool quiet = 1;
    int dev_rand = int(gsl_rng_uniform_int(ctx.rng.gslrand, 10000)); // in case traits have to be drawn at random, modified FF v.3.1.5

    //*############################################*/
    //*## First diameter (minimum data condition ##*/
    //*############################################*/

    string parameter_name = "dbh";
    string parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
    SetParameter(parameter_name, parameter_value, t_dbh, 0.001f, 15.0f, 0.0f, quiet);

    if (t_dbh > 0.0)
    {
        success = 1;

        //*####################*/
        //*## General traits ##*/
        //*####################*/

        ctx.diag.nblivetrees++;
        // assign site, by default initialized to site of tree
        t_site = site;

        // get fromData indicator, 1 by default
        t_from_Data = 1;

        // get species label first
        parameter_name = "s_name";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        int species_exists = 0;
        // species run from 1 to ctx.grid.nbspp!
        for (int s = 1; s <= ctx.grid.nbspp; s++)
        {
            if (parameter_value == ctx.S[s].s_name)
            {
                t_sp_lab = s;
                species_exists = 1;
            }
        }
        if (species_exists == 0)
        {
            t_sp_lab = int(gsl_rng_uniform_int(ctx.rng.gslrand, ctx.grid.nbspp)) + 1;
            nb_speciesrandom++;
            // cout << "Species: " << parameter_value << " not found. Initializing as random species: " << ctx.S[t_sp_lab].s_name << endl;
        }

        // update species counter
        ctx.S[t_sp_lab].s_nbind++;

        // CrownDisplacement
        parameter_name = "CrownDisplacement";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_CrownDisplacement, 0, 25, 0, quiet); // maximum displacement of 25m is unrealistic, but upper end where even largest possible crowns would still overlap with stem, defaults to 0

        //*#############################*/
        //*## Basic functional traits ##*/
        //*#############################*/

        // !!!: in future versions, this could be condensed by creating a template
        parameter_name = "Pmass";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Pmass, 0.0f, 0.01f, 0.0f, quiet); // maximum value of 0.01 based on typical values in the range of 0.0001-0.001

        if (t_Pmass == 0.0)
        {
            // draw random trait
            t_mult_P = ctx.intra.d_intraspecific_P[dev_rand];
            t_Pmass = ctx.S[t_sp_lab].s_Pmass * t_mult_P;
        }
        else
        {
            // infer multiplier trait and check against supplied values if available
            // t_mult_P = t_Pmass/ctx.S[t_sp_lab].s_Pmass;
            parameter_name = "mult_P";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_P, 0.0f, 2.0f, 1.0f, quiet);
            // example for how parameters could be checked
            //        string parameter_name_emp = "mult_P";
            //        CompareParameters(t_mult_P, 0.0f, 10.0f, 0.0f, parameter_name_emp, parameter_names, parameter_values, quiet);
        }

        parameter_name = "Nmass";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Nmass, 0.0f, 0.1f, 0.0f, quiet); // maximum value of 0.1 based on typical values in the range of 0.01-0.05

        if (t_Nmass == 0.0)
        {
            // draw random trait
            t_mult_N = ctx.intra.d_intraspecific_N[dev_rand];
            t_Nmass = ctx.S[t_sp_lab].s_Nmass * t_mult_N;
        }
        else
        {
            // infer multiplier trait and check against supplied values if available
            // t_mult_N = t_Nmass/ctx.S[t_sp_lab].s_Nmass;
            parameter_name = "mult_N";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_N, 0.0f, 2.0f, 1.0f, quiet);
        }

        parameter_name = "LMA";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_LMA, 0.0f, 1000.0f, 0.0f, quiet); // maximum value of 1000 is based on assumption that some trees may reach LMA as high as 500 (gymnosperms, some palm trees, etc.)

        if (t_LMA == 0.0)
        {
            // draw random trait
            t_mult_LMA = ctx.intra.d_intraspecific_LMA[dev_rand];
            t_LMA = ctx.S[t_sp_lab].s_LMA * t_mult_LMA;
        }
        else
        {
            // infer multiplier trait and check against supplied values if available
            // t_mult_LMA = t_LMA/ctx.S[t_sp_lab].s_LMA;
            parameter_name = "mult_LMA";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_LMA, 0.0f, 2.0f, 1.0f, quiet);
        }

#ifdef WATER
        parameter_name = "leafarea";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_leafarea, 0.0f, 1000.0f, 0.0f, quiet); // maximum value of 1000 is based on assumption that some trees may reach LMA as high as 500 (gymnosperms, some palm trees, etc.)

        if (t_leafarea == 0.0)
        {
            t_mult_leafarea = ctx.intra.d_intraspecific_leafarea[dev_rand];
            t_leafarea = ctx.S[t_sp_lab].s_leafarea * t_mult_leafarea;
        }
        else
        {
            parameter_name = "mult_leafarea";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_leafarea, 0.0f, 2.0f, 1.0f, quiet);
        }

        parameter_name = "tlp";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_tlp, -10.0f, 0.0f, 0.0f, quiet); // maximum value of 1000 is based on assumption that some trees may reach LMA as high as 500 (gymnosperms, some palm trees, etc.)

        if (t_leafarea == 0.0)
        {
            t_mult_tlp = ctx.intra.d_intraspecific_tlp[dev_rand];
            t_tlp = ctx.S[t_sp_lab].s_tlp * t_mult_tlp;
        }
        else
        {
            parameter_name = "mult_tlp";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_tlp, 0.0f, 2.0f, 1.0f, quiet);
        }

#endif

        parameter_name = "wsg";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_wsg, 0.0f, 1.5f, 0.0f, quiet); // hard upper limit of 1.5 based on density of woody substance

        if (t_wsg == 0.0)
        {
            // draw random trait
            t_dev_wsg = ctx.intra.d_intraspecific_wsg[dev_rand];
            t_wsg = fmaxf(ctx.S[t_sp_lab].s_wsg + t_dev_wsg, 0.05);
        }
        else
        {
            // infer deviation trait and check against supplied values if available
            // t_dev_wsg = t_wsg - ctx.S[t_sp_lab].s_wsg;
            parameter_name = "dev_wsg";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_dev_wsg, -1.0f, 1.0f, 0.0f, quiet); //
        }

        parameter_name = "dbhmax";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_dbhmax, 0.0f, 15.0f, 0.0f, quiet); // maximum value of 15m (Baobabs, Sequoias and Taxodium mucronatum can reach somewhere around 10m)

        if (t_dbhmax == 0.0)
        {
            t_dbhmax = ctx.S[t_sp_lab].s_dbhmax;
            t_mult_dbhmax = ctx.intra.d_intraspecific_dbhmax[dev_rand];
            t_dbhmax *= t_mult_dbhmax;
            if (t_dbhmax < t_dbh * 1.5)
            {
                t_dbhmax = t_dbh * 1.5;
                t_mult_dbhmax = t_dbhmax / ctx.S[t_sp_lab].s_dbhmax;
            }
        }
        else
        {
            // infer multiplier trait and check against supplied values if available
            // t_mult_dbhmax = t_dbhmax/ctx.S[t_sp_lab].s_dbhmax;
            parameter_name = "mult_dbhmax";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_dbhmax, 0.0f, 2.0f, 1.0f, quiet);
        }

        // t_dbhmature = t_dbhmax * 0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
        parameter_name = "dbhmature";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_dbhmature, 0.0f, 15.0f, 0.0f, quiet); // maximum value of 15m (Baobabs, Sequoias and Taxodium mucronatum can reach somewhere around 10m)

        //*##############################*/
        //*## Calculate derived traits ##*/
        //*##############################*/
        // LAImax and fraction_filled, although constant throughout tree's life will be calculated further down

        t_Vcmax = CalcVcmaxm(t_LMA, t_Nmass, t_Pmass) * t_LMA;
        parameter_name = "Vcmax";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Vcmax, 0.0f, 200.0f, 50.0f, quiet);

        // t_Jmax = CalcJmaxm() * t_LMA;
        parameter_name = "Jmax";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Jmax, 0.0f, 200.0f, 50.0f, quiet);

        t_Rdark = CalcRdark(t_LMA, t_Nmass, t_Pmass, t_Vcmax);

#ifdef WATER
        t_wleaf = sqrt(t_leafarea * 0.0001);
        t_phi_lethal = -0.9842 + 3.1795 * t_tlp; // Inferred from data provided in Bartlett et al. 2016 PNAS
        t_itlp = 1 / t_tlp;
        t_b = -2.23 * t_itlp; // this results from the following assumption: (i) ctx.params.phi(TLP)=0.97*ctx.params.phi(gs90), from Martin-StPaul et al. 2017 Ecology letters; (ii) WSF=0.1 at ctx.params.phi(gs90)
#endif

        // height allometry
        t_hmax = ctx.S[t_sp_lab].s_hmax;
        t_ah = ctx.S[t_sp_lab].s_ah;

        // leaf related traits

        // if(ctx.opt._LA_regulation > 0) CalcLAImax();

        parameter_name = "LAImax";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_LAImax, 0.0f, 10.0f, 2.0f, quiet);

        parameter_name = "leaflifespan";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_leaflifespan, 3.0f, 4000.0f, 0.0f, quiet);

        parameter_name = "lambda_young";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_lambda_young, 0.0f, 1.0f, 0.0f, quiet);

        parameter_name = "lambda_mature";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_lambda_mature, 0.0f, 1.0f, 0.0f, quiet);

        parameter_name = "lambda_old";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_lambda_old, 0.0f, 1.0f, 0.0f, quiet);

        // if((t_leaflifespan == 0.0) | (t_lambda_young == 0.0) | (t_lambda_mature == 0.0) | (t_lambda_old == 0.0)) CalcLeafLifespan(); // if the Kikuzawa model is used, the leaflifespan will be modified by a random error term, which may considerably affect tree performance if the value is recomputed

        //*#######################################*/
        //*## Traits that vary during tree life ##*/
        //*#######################################*/
        // note that this rough division encodes approximations/choices in the TROLL model to not assume changes in certain functional traits during different ontogenetic stages

        parameter_name = "age";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_age, ctx.time.timestep, 10000.0f, 1.0f, quiet);

        parameter_name = "height";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_height, 0.0f, 150.0f, 0.0f, quiet); // maximum value of 150m to leave a bit of room

        // !!!: assumption of random deviations from allometry, irrespective of neighborhood, is a problematic assumption for generic inventories
        if (t_height == 0.0)
        {
            t_mult_height = ctx.intra.d_intraspecific_height[dev_rand];
            UpdateHeight(ctx);
        }
        else
        {
            // t_mult_height = t_height/CalcHeightBaseline(t_ah, t_hmax, t_dbh);
            parameter_name = "mult_height";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_height, 0.0f, 2.0f, 1.0f, quiet);
        }

        parameter_name = "CD";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_CD, 0.0f, t_height * 0.5f, 0.0f, quiet);

        // !!!: assumption of random deviations from allometry, irrespective of neighborhood, is a problematic assumption for generic inventories
        if (t_CD == 0.0)
        {
            t_mult_CD = ctx.intra.d_intraspecific_CD[dev_rand];
            UpdateCD(ctx);
        }

        else
        {
            // t_mult_CD = t_CD/CalcCDBaseline(t_height);
            parameter_name = "mult_CD";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_CD, 0.0f, 2.0f, 1.0f, quiet);
        }

        parameter_name = "CR";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_CR, 0.0f, 25.0f, 0.0f, quiet); // maximum value of 15m (Baobabs, Sequoias and Taxodium mucronatum can reach somewhere around 10m)

        // !!!: assumption of random deviations from allometry, irrespective of neighborhood, is a problematic assumption for generic inventories
        if (t_CR == 0.0)
        {
            t_mult_CR = ctx.intra.d_intraspecific_CR[dev_rand];
            UpdateCR(ctx);
        }
        else
        {
            // t_mult_CR = t_CR/CalcCRBaseline(t_dbh);
            parameter_name = "mult_CR";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_mult_CR, 0.0f, 2.0f, 1.0f, quiet);
        }

        //*##############################*/
        //*## Calculate derived traits ##*/
        //*##############################*/

        if (ctx.opt._BASICTREEFALL)
        {
            parameter_name = "Ct";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_Ct, 0.0f, float(ctx.grid.HEIGHT), 0.0f, quiet); // the default is a height fall threshold of maximum height (so sth the tree never reaches)
            if (t_Ct == 0.0)
                t_Ct = CalcCt(ctx);
        }

        // float fraction_filled_general = 1.0 - ctx.crown.crown_gap_fraction;
        // t_fraction_filled = fminf(fraction_filled_general/(t_mult_CR * t_mult_CR),1.0);
        parameter_name = "fraction_filled";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_fraction_filled, 0.0f, 1.0f, 1.0f, quiet);

        //*######################################*/
        //*## Further photosynthetic variables ##*/
        //*######################################*/
        parameter_name = "LAmax";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_LAmax, 0.0f, 100000.0f, -1.0f, quiet); // assuming LAI cannot exceed 10 m2/m2 and maximum crown area in TROLL is 25^2 * pi, then upper limit is around 20000, 100000 to be safe

        parameter_name = "LA";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_LA, 0.0f, 100000.0f, -1.0f, quiet); // assuming LAI cannot exceed 10 m2/m2 and maximum crown area in TROLL is 25^2 * pi, then upper limit is around 20000, 100000 to be safe - note that LA can be larger than LAmax when tree has been recently shaded and leaf turnover is slow

        parameter_name = "youngLA";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_youngLA, 0.0f, 100000.0f, -1.0f, quiet);

        parameter_name = "matureLA";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_matureLA, 0.0f, 100000.0f, -1.0f, quiet);

        parameter_name = "oldLA";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_oldLA, 0.0f, 100000.0f, -1.0f, quiet);

        parameter_name = "LAI";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_LAI, 0.0f, 10.0f, 0.5f, quiet);

        parameter_name = "litter";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_litter, 0.0f, 100000.0f, -1.0f, quiet);

        float crown_area = PI * t_CR * t_CR;
        float crown_area_nogaps = GetCrownAreaFilled(crown_area);

        /* if(t_LAmax >= 0.0 && t_LA >= 0.0 && t_youngLA >= 0.0 && t_matureLA >= 0.0 && t_oldLA >= 0.0 && t_litter >= 0.0){
            // there could be more consistency checks here (e.g. t_youngLA + t_matureLA + t_oldLA should be approximately t_LA)
            t_LAI = t_LA/crown_area_nogaps;
        } else {
            if(ctx.opt._LA_regulation > 0){
                float LAIexperienced_eff;
                CalcLAmax(LAIexperienced_eff, t_LAmax);     // !!!: problematic calculation if only classic inventory data is provided (i.e. x/y/dbh/s_name): we calculate the maximum leaf area for an empty canopy, thus grossly overestimating it for trees in the understory --> first improvement would be implementing calculation of this variable from highest tree to smallest tree and successively allocating leaves. Cf. the currently not used CalcLAinitial(). Interestingly, when trying this, the effects were not very large, probably because they are overshadowed by the impact of random allometric deviations and crown overlap.
                t_LA = t_LAmax;   // assume half the maximum leaf area?
                t_LAI = t_LA/crown_area_nogaps;
            } else {
#ifdef CROWN_UMBRELLA
                t_LAI = ctx.params.dens * fminf(t_CD, 3.0);
#else
                t_LAI = ctx.params.dens * t_CD;
#endif
                t_LA = t_LAI * crown_area_nogaps;
            }
            InitialiseLeafPools();
        }*/

        parameter_name = "sapwood_area";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        float ba = t_dbh * t_dbh * 0.25 * PI;
        SetParameter(parameter_name, parameter_value, t_sapwood_area, 0.0f, ba, 0.0f, quiet);

        if (t_sapwood_area == 0.0)
        {
            float ddbh = fminf(t_dbh, 0.04f); // limit the approximate sapwood thickness to 5cm
            UpdateSapwoodArea(ctx, ddbh);
        }

        if (ctx.opt._LA_regulation > 0)
        {
            // float carbon_storage_max = CalcCarbonStorageMax();
            parameter_name = "carbon_storage";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_carbon_storage, 0.0f, 1000000.0f, -1.0f, quiet);
            // SetParameter(parameter_name, parameter_value, t_carbon_storage, 0.0f, carbon_storage_max, -1.0f, quiet);
            // if(t_carbon_storage == -1.0) t_carbon_storage = carbon_storage_max * 0.5; // We here initialize trees with half their maximum storage

            parameter_name = "carbon_biometry";
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            // SetParameter(parameter_name, parameter_value, t_carbon_biometry, 0.0f, carbon_storage_max * 10.0f, 0.0f, quiet);
            SetParameter(parameter_name, parameter_value, t_carbon_biometry, 0.0f, 1000000.0f * 10.0f, 0.0f, quiet);
        }

        parameter_name = "hurt";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        unsigned short hurt_min = 0, hurt_max = USHRT_MAX, hurt_default = 0;
        SetParameter(parameter_name, parameter_value, t_hurt, hurt_min, hurt_max, hurt_default, quiet);

        parameter_name = "NPPneg";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_NPPneg, 0, INT_MAX, 0, quiet);

        parameter_name = "multiplier_seed";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_multiplier_seed, 1, INT_MAX, 1, quiet);

        //*##############################*/
        //*## "Convenience" parameters ##*/
        //*##############################*/
        // these parameters can be set to 0, as they are recomputed each ctx.time.timestep in the actual light environment
        // will be computed for initial configuration after allocation of leaves to voxel field
        t_GPP = t_NPP = t_Rday = t_Rnight = t_Rstem = 0.0;

        parameter_name = "GPP";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_GPP, -10.0f, 5000.0f, 0.0f, quiet);

        parameter_name = "NPP";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_NPP, 0.0f, 1000.0f, 0.0f, quiet);

        parameter_name = "Rday";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Rday, 0.0f, 1000.0f, 0.0f, quiet);

        parameter_name = "Rnight";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Rnight, 0.0f, 1000.0f, 0.0f, quiet);

        parameter_name = "Rstem";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Rstem, 0.0f, 1000.0f, 0.0f, quiet);

        parameter_name = "Rdark";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Rdark, 0.0f, 1000.0f, 0.0f, quiet);

#ifdef WATER
        //*###########*/
        //*## water ##*/
        //*###########*/

        t_transpiration = 0.0;

        parameter_name = "transpiration";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_transpiration, 0.0f, 1.0f, 0.0f, quiet);

        parameter_name = "root_depth";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_root_depth, 0.0f, 100.0f, 2.0f, quiet);

        parameter_name = "phi_root";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_phi_root, -15.0f, 0.0f, 0.0f, quiet);

        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            char par_name[30];
            sprintf(par_name, "root_biomass%d", l);
            string par_name_s = par_name;
            parameter_name = par_name_s;
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_root_biomass[l], 0.0f, 100000.0f, 1.0f, quiet);
        }

        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            char par_name[30];
            sprintf(par_name, "soil_layer_weight%d", l); //
            string par_name_s = par_name;
            parameter_name = par_name_s;
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_soil_layer_weight[l], 0.0f, 1.0f, 1.0f, quiet);
        }

        parameter_name = "g1_0"; // not necessarily needed currently, since computable from other traits, but this can be changed given the low confindence in these relationships (and should at least be recomputed here anyway)
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_g1_0, 0.0f, 10.0f, 3.77f, quiet);

        parameter_name = "WSF";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_WSF, 0.0f, 1.0f, 1.0f, quiet);

        // t_g1=t_g1_0*t_WSF;
        parameter_name = "g1"; // not necessarily needed currently, since computable from g1_0 and WSF
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_g1, 0.0f, 10.0f, 3.77f, quiet);

        parameter_name = "WSF_A";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_WSF_A, 0.0f, 1.0f, 1.0f, quiet);

        parameter_name = "Ndays_dry";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Ndays_dry, 0, 1000, 0, quiet);

        parameter_name = "Ndays_wet";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_Ndays_wet, 0, INT_MAX, 1, quiet);

        parameter_name = "pheno_factor";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_pheno_factor, 0.0f, 1.0f, 1.0f, quiet);

#endif
        //*################*/
        //*## for output ##*/
        //*################*/
#ifdef Output_ABC
        // parameter_name = "dbh_previous";
        // parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        // SetParameter(parameter_name, parameter_value, t_dbh_previous, 0.0f, t_dbh, t_dbh, quiet);
        t_dbh_previous = t_dbh; // check whether proper initialisation is needed
#endif
#ifdef TRACK_INDIVIDUALS
        StartTracking(ctx);
#endif
    }
    return (success);
}

// v.3.1: new function to initialize trees more realistically from inventories, not used yet
// void Tree::CalcLAinitial(){
//    if(t_LA < 0.0){
//        // 1. determine leaf area and related variables, based on current environment
//        float crown_area = PI*t_CR*t_CR;
//        float crown_area_nogaps = GetCrownAreaFilled(crown_area);
//
//        if(ctx.opt._LA_regulation > 0){
//            float LAIexperienced_eff;
//            CalcLAmax(LAIexperienced_eff, t_LAmax);
//            t_LA = t_LAmax;   // assume half the maximum leaf area?
//            t_LAI = t_LA/crown_area_nogaps;
//        } else {
// #ifdef CROWN_UMBRELLA
//            t_LAI = ctx.params.dens * fminf(t_CD, 3.0);
// #else
//            t_LAI = ctx.params.dens * t_CD;
// #endif
//            t_LA = t_LAI * crown_area_nogaps;
//        }
//        InitialiseLeafPools();
//
//        if(t_sapwood_area == 0.0){
//            float ddbh = fminf(t_dbh, 0.04f);   // initial sapwood thickness of ca. 2.5cm
//            UpdateSapwoodArea(ddbh);
//        }
//
//        // 2. allocate the leaves to the ctx.field.LAI3D field
//        // nota bene: we here use the same function as CalcLAI3D, with one exception: LAI2dens_cumulated instead of LAI2dens; this means that we allocate the cumulated LAI in each ctx.diag.layer and do not require any summation afterwards
//        int site_crowncenter = t_site + t_CrownDisplacement;
//        int row_crowncenter = site_crowncenter/ctx.grid.cols;
//        int col_crowncenter = site_crowncenter%ctx.grid.cols;
//
//        float LA_cumulated = 0.0;   //Currently, an output variable is required by LoopLayerUpdateCrownStatistic_template, we here use LA_cumulated as control variable
//
//        int crown_top = int(t_height);
//        int crown_base = int(t_height - t_CD);
//        int max_shells = min(crown_top - crown_base + 1, 4);    //since the new crown shapes
//        float fraction_filled_target = t_fraction_filled;
//
//        for(int shell_fromtop = 0; shell_fromtop < max_shells; shell_fromtop++){
//            LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, t_height, t_CR, t_CD, fraction_filled_target, shell_fromtop, GetRadiusSlope, t_LAI, LA_cumulated, LAI2dens_cumulated, UpdateLAI3D);
//        }
//    }
//}

#ifdef WATER

// ##################################################
//  Computation of the tree water availability in the root zone
// ####################################################
//! - Tree::Water_availability could be divided into two functions:
//! -# Water_availability: the equivalent of Tree::Fluxh but for water, ie compute the tree water conditions, and which updates t_phi_root and t_WSF
//! -# another function that could be called UpdateRootDistribution, similar to UpdateLeafDynamics and UpdateTreeBiometry, that would update t_root_depth and t_soil_layer_weight and t_root_biomass, at the end of Tree::Growth (and in Tree::Birth)
//! The current use of one function was chosen for simplicity's sake, but we could shift to a version with 2 functions if we prefer similarity to what has been already done
void Tree::Water_availability(Context &ctx)
{

    // if(t_LA > 0.0){       // TO BE CHECKED: what happens to a tree with leafarea=0 in the following timesteps ?
    // Tree root biomass
    float total_root_biomass = t_LA * t_LMA;
    // In the absence of a more explicit allocation scheme, we assume that the total root biomass=the total leaf biomass. This may be discussed and changed. Eg. in Schippers et al. 2015 Functional Plant Biology, they assume that root biomass is 1.4 times leaf biomass... More importantly, the tree-level carbon budget should be closed and balanced. Needs to be specified that this correponds to the total FINE root biomass. Note that similarly, in ED2 (Xu et al. 2016), a constant ratio of 1 is assumed between leaf biomass and fine root biomass (whereas a constant ratio of 0.25 is assumed between coarse root and aboveground stem biomass). But leaf and fine root phenology follow slightly different dynamics (see Xu et al. 2016 SI) so unclear how this is consistent...?

    // Tree root depth
    // t_root_depth=fminf(ctx.soil.layer_depth[(ctx.soil.nblayers_soil-1)],pow((0.001*total_root_biomass),0.8)*0.2173913);       //this is the tree root depth, computed according to Arora & Boer 2003 Earth interactions (see equ. 11 therein). Instead of assuming a fixed root exponential profile, that typically results in a fixed root depth from recruitment to death, hence overestimating root vertical growth and depth for small saplings, they derived a root distribution profile that depends on root biomass. The simulated root depth and distribution have to be checked. 0.2173913=3/b with b=13.8
    // t_root_depth=fminf(ctx.soil.layer_depth[(ctx.soil.nblayers_soil-1)],pow((0.001*total_root_biomass),0.1)*0.7);
    // t_root_depth=fmaxf(2.0, fminf(ctx.soil.layer_depth[(ctx.soil.nblayers_soil-1)],pow((0.001*total_root_biomass),0.1)*0.7));
    // t_root_depth=fminf(ctx.soil.layer_depth[(ctx.soil.nblayers_soil-1)], 0.06410256*t_height+1.435897);  // this dependency of tree root depth with tree height is highly questionable -- cf. eg. Stahl et al. 2013 Oecologia, but see also Brum et al. 2019 Journal of Ecology-- but parcimonious. Note that neither Stahl et al. 2013 nor Brum et al. 2019 reported actual root depth, but the "depth of water uptake". Should probably depends more on the total soil depth, and average soil water availabilty -- cf. distribution of maximal root depth per biomes in Canadell/Jackson et al.
    // t_root_depth=0.22*pow(t_dbh*100, 0.54); // this is the allometry used in ED2, Xu et al. 2016, for evergreen trees. Parameter values were fixed based on Kenzo et al. 2009 Journal of Tropical Biology, which reported on data from excavated trees in wet secondary forests in Malaysia, although Xu et al. 2016 applied it in a seasonally dry tropical forests in Costa Rica. This allometry was then compared against the one obtained from excavated trees in a seasonally dry tropical forests in Costa Rica by Smith-Martin et al. 2019 New Phytologist: although the shape and the differences among the allometries for deciduous and evergreen trees were in overall agreement with observations, the ED2 allometry generally underestimated the observed root depth. As it seems that dry forest species have deeper root than wet forest species (see e.g. Holbrook et al. 1995, cited in Smith-Martin et al. 2019, or Markesteijn & Poorter 2009), this should be better suited to our application to wet forest although this remains to be confronted with data. This allometry gives a root depth of 22cm and 2.65 ctx.params.m for trees with dbh=1cm and 1m respectively. Note that in Xu et al., the first parameter (here 0.22) varied among phenological types, so that root depth of evergreen trees is about twice the ones of deciduous. This link between plant phenology and rooting depth was found empirically in Smith-Martin et al. 2019 New Phytologist, Hasselquist et al. 2010 Oecologia, while Markesteijn & Poorter found a correlation between root depth and stem density on first year seedling in Panama. This link between phenology and rooting depth with higher rooting depth for evergreen species was needed to sustain leaf cover in simulations with ED2, despite their higher leaf and stem drought tolerance (P50 and TLP; Smith-Martin et al. 2019). In absence of data, we here assumed this allometry suits to all species but this has to be discussed !
    t_root_depth = fminf(0.35 * pow(t_dbh * 100, 0.54), ctx.soil.layer_depth[(ctx.soil.nblayers_soil - 1)]); // this is the allometry used in ED2-hydro, Xu et al. 2016, for evergreen trees, with one parameter, b1=root_depth(dbh=1cm), changed to correct the overall underestimation of tree rootdepth it leads to (cf. Smith-Martin et al. 2019 New phytol.). In absence of data on root depth (in general, and even more at the species level, see email with V. Freycon & B. Ferry january 20th 2020), this should probably be calibrated/fine-tuned.
    float i_root_depth = 1.0 / t_root_depth;

    // Tree water potential in the root zone (or predawn water potential)
    float shallow_bound = 1.0;
    float deep_bound = 1.0;
    float sumG = 0.0;
    t_phi_root = 0.0;

    float layer_depth_previous = 0.0;
    float layer_thickness = 0.0;
    float root_area = t_CR * t_CR * PI;

    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {

        layer_thickness = ctx.soil.layer_depth[l] - layer_depth_previous;
        layer_depth_previous = ctx.soil.layer_depth[l];

        deep_bound = exp(-3.0 * i_root_depth * ctx.soil.layer_depth[l]);
        // deep_bound=exp(-4.02*i_root_depth*ctx.soil.layer_depth[l]); // by consistency with the root depth allometry, we here used the parameter value used by Xu et al. 2016 for evergreen trees. This value (log(beta)) was defined based on Jackson et al. 1996 Oecologia.
        // t_soil_layer_weight[l]=shallow_bound-deep_bound;                        //the weight for each soil ctx.diag.layer are the relative tree root biomass in this ctx.diag.layer. This is the most classic approach, however probably not fully relevant, inasmuch as trees may equilibrate "preferentially" with the wettest part of the soil, where water is more available and more easily exctractable. Hence another possible weighting integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2011). See below.
        t_root_biomass[l] = total_root_biomass * (shallow_bound - deep_bound); // this is the root biomass in ctx.diag.layer l, computed from the integration of the root biomass exponential profile in this ctx.diag.layer, following Arora & Boer 2003. This also corresponds to the shape of cumulative root fraction used by Jackson et al. 1996 Oecologia, who gathered a large global database of root distributions with this equation. It was thus used in many models, e.g. Duursma & Medlyn 2012; Xu et al. 2016.
        // t_soil_layer_weight[l]=t_root_biomass[l]*10.0/(-log(sqrt(PI*t_root_biomass[l]*10.0)*0.001)); // this soil ctx.diag.layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)

        if (t_LA > 0.0)
        {

            if (ctx.opt._WATER_RETENTION_CURVE == 1)
            {
                if (ctx.opt._SOIL_LAYER_WEIGHT == 0)
                { // soil ctx.diag.layer weights as a function of root biomass in each ctx.diag.layer only (cf. M1 in de Kauwe et al 2015)
                    t_soil_layer_weight[l] = t_root_biomass[l];
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];
                }
                else if (ctx.opt._SOIL_LAYER_WEIGHT == 1)
                {
                    t_soil_layer_weight[l] = t_root_biomass[l] * 10.0 * ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]] / (abs(log(sqrt(PI * t_root_biomass[l] * 10.0) * 0.001))); // this soil ctx.diag.layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];                                                                // the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil ctx.diag.layer in the DCELL where the tree stands. Note that ctx.soil.KsPhi was here not computed as ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]*ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                }
                else if (ctx.opt._SOIL_LAYER_WEIGHT == 2)
                { // soil ctx.diag.layer weight as in Duursma & Medlyn 2012 (cf M3 in de Kauwe et al. 2015)

                    // if (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=t_root_biomass[l]*10.0*ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0)*0.001)))*(ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]]+3); // this soil ctx.diag.layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012). Note that in their implementation of MAESPA, Christina et al. used minimum root water potential = -1.6 MPa and not -3 MPa as here, and they also added a gravimetric component, since they explore the effect of very deep root. Sensibility to this value and addition to be tested.
                    // if (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=t_root_biomass[l]*10.0*ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0/(ctx.grid.length_dcell*ctx.grid.length_dcell*layer_thickness))*0.001)))*(ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]]+3);
                    if (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] > (-3.0))
                        t_soil_layer_weight[l] = (t_root_biomass[l] * 10.0 / root_area) * ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]] / (abs(log(sqrt(PI * t_root_biomass[l] * 10.0 / (root_area * layer_thickness)) * 0.001))) * (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] + 3);
                    else
                        t_soil_layer_weight[l] = 0.0;
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]]; // the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil ctx.diag.layer in the DCELL where the tree stands. Note that ctx.soil.KsPhi was here not computed as ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]*ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                }
            }
            else if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {

                if (ctx.opt._SOIL_LAYER_WEIGHT == 0)
                { // soil ctx.diag.layer weights as a function of root biomass in each ctx.diag.layer only (cf. M1 in de Kauwe et al 2015)
                    t_soil_layer_weight[l] = t_root_biomass[l];
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];
                }
                else if (ctx.opt._SOIL_LAYER_WEIGHT == 1)
                {
                    t_soil_layer_weight[l] = t_root_biomass[l] * 10.0 / (abs(log(sqrt(PI * t_root_biomass[l] * 10.0) * 0.001))); // this soil ctx.diag.layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]];                                         // the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil ctx.diag.layer in the DCELL where the tree stands. Note that ctx.soil.KsPhi was here not computed as ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]*ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                    t_soil_layer_weight[l] *= ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]];
                }
                else if (ctx.opt._SOIL_LAYER_WEIGHT == 2)
                { // soil ctx.diag.layer weight as in Duursma & Medlyn 2012 (cf M3 in de Kauwe et al. 2015)

                    if (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] > (-3.0))
                        t_soil_layer_weight[l] = (t_root_biomass[l] * 10.0 / root_area) / (abs(log(sqrt(PI * t_root_biomass[l] * 10.0 / (root_area * layer_thickness)) * 0.001))) * (ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] + 3); // this soil ctx.diag.layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012). Note that in their implementation of MAESPA, Christina et al. used minimum root water potential = -1.6 MPa and not -3 MPa as here, and they also added a gravimetric component, since they explore the effect of very deep root. Sensibility to this value and addition to be tested.
                    else
                        t_soil_layer_weight[l] = 0.0;
                    t_phi_root += t_soil_layer_weight[l] * ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]]; // the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil ctx.diag.layer in the DCELL where the tree stands. Note that ctx.soil.KsPhi was here not computed as ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]*ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                    t_soil_layer_weight[l] *= ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]];
                }
            }

            if (isnan(t_phi_root))
            {
                cout << "nan t_phi_root: " << endl;
                cout << l << "\t" << t_soil_layer_weight[l] << "\t" << ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]] << ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] << "\t" << endl;
            }
        }
        else
        {
            if (l == 0)
            {
                t_soil_layer_weight[l] = 1;
                t_phi_root += ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];
            }
            else
                t_soil_layer_weight[l] = 0;
        }

        // if(t_soil_layer_weight[l]<=0) cout << "t_soil_layer_weight[l]=" << t_soil_layer_weight[l] << " t_root_biomass[l]=" << t_root_biomass[l] << " ctx.soil.Ks=" << ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]] << " ctx.soil.Ks*Phi soil =" << ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]] << " ctx.params.phi soil=" <<   ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] << endl;
        // t_soil_layer_weight[l]=t_root_biomass[l]*10/(-log(sqrt(PI*t_root_biomass[l]*10)*0.001));
        // t_phi_root+=t_soil_layer_weight[l]*fmaxf(0.0,(KsPhi2[l][ctx.grid.site_DCELL[t_site]]-ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]]*t_s->s_tlp));  //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil ctx.diag.layer
        // t_soil_layer_weight[l]*=fmaxf(0.0,(ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]]-ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]]*t_s->s_tlp));
        sumG += t_soil_layer_weight[l];
        shallow_bound = deep_bound;

        if (t_soil_layer_weight[l] < 0.0 || isnan(t_soil_layer_weight[l]))
        {
            cout << "incorrect soil_layer_weight: " << endl;
            cout << l << "\t" << ctx.grid.site_DCELL[t_site] << "\t" << t_soil_layer_weight[l] << "\t" << t_root_biomass[l] << "\t" << ctx.soil.Ks[l][ctx.grid.site_DCELL[t_site]] << "\t" << ctx.soil.KsPhi[l][ctx.grid.site_DCELL[t_site]] << "\t" << -log(sqrt(PI * t_root_biomass[l] * 10) * 0.001) << "\t" << deep_bound << "\t" << shallow_bound << "\t" << ctx.soil.SWC3D[l][ctx.grid.site_DCELL[t_site]] << "\t" << (ctx.soil.SWC3D[l][ctx.grid.site_DCELL[t_site]] - ctx.soil.Min_SWC[l]) / (ctx.soil.Max_SWC[l] - ctx.soil.Min_SWC[l]) << "\t" << ctx.soil.Ksat[l] << "\t" << ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];
            cout << endl;
        }
    }

    if (sumG > 0.0)
    {
        float isumG = 1.0 / sumG;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            t_soil_layer_weight[l] *= isumG;
        }
        t_phi_root *= isumG;
        if (isnan(t_phi_root))
        {
            cout << "nan t_phi_root" << endl;
            cout << isumG << endl;
        }
    }
    else
    {
        float iRootB = 1.0 / total_root_biomass;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            t_phi_root += t_root_biomass[l] * ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]];
            if (isnan(t_phi_root))
            {
                cout << "nan t_phi_root" << endl;
                cout << t_root_biomass[l] << "\t" << ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] << endl;
            }
            t_soil_layer_weight[l] = t_root_biomass[l] * iRootB;
        }
        t_phi_root *= iRootB;
        if (isnan(t_phi_root))
        {
            cout << "nan t_phi_root" << endl;
            cout << iRootB << endl;
        }
    }

    float sum_weight = 0.0;
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
    {
        sum_weight += t_soil_layer_weight[l];
    }
    if (sum_weight > 1.5)
    {
        cout << "Weird weights in Water_availability; sum_weight=" << sum_weight << " sumG=" << sumG << " total_root_biomass=" << total_root_biomass << " t_phi_root=" << t_phi_root;
        if (sumG > 0.0)
        {
            float isumG = 1.0 / sumG;
            cout << " isumG=" << isumG << endl;
        }
        cout << "t_soil_layer_weight[l]: ";
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            cout << t_soil_layer_weight[l] << "\t";
        }
        cout << endl;
        cout << "soil_phi: ";
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            cout << ctx.soil.soil_phi3D[l][ctx.grid.site_DCELL[t_site]] << "\t";
        }
        cout << endl;
    }

    t_phi_root -= 0.01 * t_height; // gravitational effect on leaf water potential and stress

    if (isnan(t_phi_root))
    {
        cout << "nan t_phi_root" << endl;
        cout << t_height << endl;
    }

    //! Tree water stress factors

    // t_WSF=fminf(1.0, fmaxf(0.0, ((ctx.soil.SWC3D[0][ctx.grid.site_DCELL[t_site]]-ctx.soil.Min_SWC[0])/(ctx.soil.Max_SWC[0]-ctx.soil.Min_SWC[0])))); // this is the simple linear WSF, with SWC as independent varible, used in lots of model (see Powell et al. 2013 New Phytol, De Kauwe et al. 2015 Biogeosciences, Laio et al. 2001 Advances in Water resources, Egea et al. 2011 AFM....)
    // t_WSF=fminf(1.0, fmaxf(0.0, (t_phi_root-(t_s->s_tlp))*(t_s->s_dWSF))); // this is the linear WSF, with wtare potential as independent variable, used in CLM model (see Powell et al. 2013 New Phytologist, Verhoef & Egea 2011 AFM)
    // t_WSF=exp((log(0.05)/t_s->s_tlp)*t_phi_root); // this is the WSF, to simulate stomatal limitation (ie. hinder ctx.params.g1), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences. If this version of WSF is finally adopted, the parameter b=log(0.05)/t_s->s_tlp, should be declare as species variable (instead of s_dWSF).
    t_WSF = exp(t_b * t_phi_root); // this is the WSF shape used to simulate stomatal limitation (ie. hinder ctx.params.g1), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences, but with a different parameterization for the parameter b (using the relationship between phi_gs90 and tlp from Martin-StPaul et al. 2017 Ecology letters, and assuming the WSF=0.9 at phi_gs90).
    // float par=t_s->s_tlp+1;
    // t_WSF_A=(1.0+exp(6.0*par))/(1.0+exp(6.0*(par-t_phi_root))); // this is the WSF, to simulate non-stomatal limitation (ie. hinder Vcmax and Jmax), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences
    // t_WSF_A=(1.0+exp(3.0*par))/(1.0+exp(3.0*(par-t_phi_root))); // this is the WSF, to simulate non-stomatal limitation (ie. hinder Vcmax and Jmax), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences
#ifdef WATER
    // t_g1= (-3.97 * t_wsg + 6.53)*t_WSF; this is the relationship provided by Lin et al. 2015
    t_g1 = t_g1_0 * t_WSF; // with the water stress factor added
#endif

    t_WSF_A = 1 / (1 + pow(t_phi_root * t_itlp, 6)); // this is the WSF used to hinder Vcmax and Jmax in Xu et al. 2016 (equ. S5) (ED2-hydro).

    if (t_WSF < 0.0 || t_WSF_A < 0.0 || t_WSF > 1.0 || t_WSF_A > 1.0 || t_phi_root > 0.0 || isnan(t_WSF) || isnan(t_WSF_A) || isnan(t_phi_root))
    {
        cout << "incorrect value in one of WSF, WSF_A, t_phi_root " << endl;
        cout << t_WSF << "\t" << t_phi_root << "\t" << t_tlp << "\t" << t_dbh << "\t" << t_height << "\t" << t_age << "\t" << t_phi_lethal << "\t" << sumG;
        cout << endl;
    }
    //}
}

// ###################################################################
//  Contribution of trees to the stand ctx.soil.Transpiration field. Called by UpdateField
// ####################################################################
//! - Adds up each tree contribution to the stand ctx.soil.Transpiration field, that is the water removed from the soil through all tree transpiration.
//! - Similar to CalcLAI, that adds up each tree contribution to ctx.field.LAI3D field.
void Tree::Water_uptake(Context &ctx)
{
    if (t_age > 0.0)
    {
        int l = 0;
        float depth = 0.0;
        float sum_weights = 0.0;
        while (sum_weights < 1.0 && l < ctx.soil.nblayers_soil)
        {
            ctx.soil.Transpiration[l][ctx.grid.site_DCELL[t_site]] += t_soil_layer_weight[l] * t_transpiration;
            if (t_soil_layer_weight[l] < 0.0 || t_transpiration < 0.0)
            {
                cout << setprecision(10);
                cout << "Problem with soil_layer_weight and transpiration at site: " << t_site << " ctx.diag.layer: " << l << " depth: " << depth << endl;
                cout << t_soil_layer_weight[l] << "\t" << t_transpiration << "\n";
            }
            sum_weights += t_soil_layer_weight[l];
            depth = ctx.soil.layer_depth[l];
            l++;
        }
    }
}

#endif

// ##################################################
//  Computation of PPFD right above the tree -- called by Tree::Birth and Tree::Growth
// ####################################################
//! Mean light flux received by the tree crown ctx.diag.layer at height h new version (PPFD symmetrical with T and VPD); PPFD, VPD, T and leafarea of a ctx.diag.layer are now local variables (not tree variables); v2.3.0, updated in v.2.5
//! - v.2.3.: Tree::Fluxh() computes the average light flux received by a tree crown ctx.diag.layer at height h , and also the average VPD and T (modified 1/02/2016)
//! - modified in v.2.4: weighting of each ctx.diag.layer's GPP by the leafarea in the ctx.diag.layer, this is also computed here
//! - modified in v.2.4 and v.2.5: introducing an alternative crown shape, "umbrella"-like, inspired by previous shell models and similar to the crown shapes in the PPA. If activated, crowns contain three layers of vegetation that, once the crown goes beyond 3m in depth, will bend downwards on the edges with a linear slope. Since v.2.5 all loops (CalcLAI, Fluxh, leafarea_max) are executed through the same template. This allows to implement other crown shapes in the future and ensures that modifications are carried through across the code
#ifdef CROWN_UMBRELLA
#ifdef WATER
void Tree::Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer, float &PPFD_incident, float &ExtinctLW)
{
#else
void Tree::Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer)
{
#endif

    int site_crowncenter = t_site + t_CrownDisplacement;
    int row_crowncenter = site_crowncenter / ctx.grid.cols;
    int col_crowncenter = site_crowncenter % ctx.grid.cols;
    int shell_fromtop = int(t_height) + 1 - h;

#ifdef WATER
    float canopy_environment_cumulated[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#else
    float canopy_environment_cumulated[4] = {0.0, 0.0, 0.0, 0.0};
#endif
    float fraction_filled_target = t_fraction_filled;

    LoopLayerUpdateCrownStatistic_template(ctx, row_crowncenter, col_crowncenter, t_height, t_CR, t_CD, fraction_filled_target, shell_fromtop, [&ctx](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); }, t_LAI, canopy_environment_cumulated, LAI2dens, [&ctx](int h, int s, float d, auto &ce){ GetCanopyEnvironment(ctx, h, s, d, ce); });

    float LA_layer = canopy_environment_cumulated[0];
    float iLA_layer;
    if (LA_layer > 0.0)
        iLA_layer = 1.0 / LA_layer;
    else
        iLA_layer = 0.0;

    leafarea_layer = LA_layer;
    PPFD = canopy_environment_cumulated[1] * iLA_layer;
    VPD = canopy_environment_cumulated[2] * iLA_layer;
    Tmp = canopy_environment_cumulated[3] * iLA_layer;
#ifdef WATER
    PPFD_incident = canopy_environment_cumulated[4] * iLA_layer;
    ExtinctLW = canopy_environment_cumulated[5] * iLA_layer;
#endif

    // if (PPFD <=0) {
    //  cout << "Warning PPFD <=0 in Fluxh !! PPFD=" << PPFD << "; iLA_layer=" << iLA_layer << "; canopy_environment_cumulated[1]=" << canopy_environment_cumulated[1] << "; VPD=" << VPD << "; Tmp=" << Tmp << "; h=" << h << "; Tree height=" << t_height << "; Tree_crowndepth=" << t_CD << "; Tree_leafarea=" << t_LA << endl;
    //}
}

#else
void Tree::Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer)
{

    int site_crowncenter = t_site + t_CrownDisplacement;
    int row_crowncenter = site_crowncenter / ctx.grid.cols;
    int col_crowncenter = site_crowncenter % ctx.grid.cols;

    float crown_area = fmaxf(PI * t_CR * t_CR, 0.0);
    int crown_intarea = int(crown_area);      // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea, 1);    // minimum area of crown (1)
    crown_intarea = min(crown_intarea, 1963); // maximum area of crown (radius 25), int(3.14*25*25)

    float fraction_filled_target = t_fraction_filled;
    float fraction_filled_actual = 0.0;

    int crown_intarea_gaps = 0;
    int crown_intarea_allocated_nogaps = 0;

    // loop over LookUp table, until crown_intarea is reached
    for (int i = 0; i < crown_intarea; i++)
    {
        if (fraction_filled_actual > fraction_filled_target)
        {
            fraction_filled_actual = (fraction_filled_actual * float(i)) / (float(i) + 1.0);
            crown_intarea_gaps++;
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
                float absorb_prev = ctx.field.LAI3D[h][site];
                float absorb_current = ctx.field.LAI3D[h - 1][site];
                float absorb_delta = absorb_current - absorb_prev;
                if (absorb_delta < 0.0)
                    absorb_delta = 0.0; //! eliminate rounding errors
                int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);
                PPFD += ctx.climate.WDailyMean * ctx.lookup.LookUp_flux_absorption[intabsorb];
                VPD += ctx.climate.VPDDailyMean * ctx.lookup.LookUp_VPD[intabsorb];
                Tmp += ctx.climate.tDailyMean - ctx.lookup.LookUp_T[intabsorb];
                crown_intarea_allocated_nogaps++;
            }
        }
    }

    // weighting has to be done via the actually allocated crown area (integer), not the floating point value
    float icrown_intarea_allocated = 1.0 / float(crown_intarea_allocated_nogaps);

    PPFD *= icrown_intarea_allocated;
    VPD *= icrown_intarea_allocated;
    Tmp *= icrown_intarea_allocated;

    // weighting of each ctx.diag.layer
    float crown_area_nogaps;
    if (fraction_filled_actual > fraction_filled_target)
        crown_area_nogaps = float(crown_intarea - crown_intarea_gaps);
    else
        crown_area_nogaps = crown_area - float(crown_intarea_gaps);

    int h_index = h - 1;
    int crown_top = int(t_height);
    int crown_base = int(t_height - t_CD);
    float fraction_abovetop = t_height - float(crown_top);
    float fraction_belowbase = float(crown_base + 1) - (t_height - t_CD);
    float dens_avg = t_LAI / t_CD;

    float dens_layer = dens_avg;
    if (crown_top == crown_base)
        dens_layer *= t_CD;
    else if (h_index == crown_top)
        dens_layer *= fraction_abovetop;
    else if (h_index == crown_base)
        dens_layer *= fraction_belowbase;

    // given constant density within ctx.diag.layer, we just need to multiply by area
    leafarea_layer = dens_layer * crown_area_nogaps;
}
#endif

// #############################
//  Empirical functions (traits)
// #############################

float Tree::CalcJmaxm(Context &ctx)
{
    float SLA = 10000.0 / t_LMA;                                                                                                                                // in cm2 g-1
    float Jmaxm = pow(10.0, fminf((-1.50 + 0.41 * log10(t_Nmass * 1000.0) + 0.45 * log10(SLA)), (-0.74 + 0.44 * log10(t_Pmass * 1000.0) + 0.32 * log10(SLA)))); // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7). s_Nmass and s_Pmass are given in g g-1, but should be in mg g-1 in equ 2 in Domingues et al. 2010, hence the mutiplication by 1000.
    return (Jmaxm);                                                                                                                                             //  in micromol C ctx.params.m-2 s-1
}

// #######################
//  Death rate calculation
// #######################
float Tree::DeathRateNDD(Context &ctx, float dbh, float carbon_starv, float ndd)
{
    float dr = 0;
    float basal = ctx.params.m * (1 - t_wsg);
    float dd = ctx.crown.deltaD * ndd * (1 - 2 * dbh / t_dbhmax);

    dr = basal;
    if (ctx.opt._LA_regulation == 0)
    {
        if (carbon_starv > t_leaflifespan)
            dr += 1.0 / ctx.time.timestep;
    }
    else
    {
        if (carbon_starv <= 0.0 && t_NPP <= 0.0)
            dr += 1.0 / ctx.time.timestep; // newIM 2021: carbon starvation occurs when the carbon stocks has been completly depleted, and carbon_starv represents t_carbon_storage (while it represents t_NPPneg when ctx.opt._LA_regulation==0)
    }
    if (dd > 0)
        dr += dd;
    return dr * ctx.time.timestep;
}

#ifdef WATER
float Tree::DeathRate(Context &ctx, float dbh, float carbon_starv, float phi_root)
{
    float dr = 0;
    float basal = fmaxf(ctx.params.m - ctx.params.m1 * t_wsg, 0.0);

    dr = basal;
    // if (carbon_starv > t_leaflifespan) dr+=1.0/ctx.time.timestep;
    if (ctx.opt._LA_regulation == 0)
    {
        if (carbon_starv > t_leaflifespan)
            dr += 1.0 / ctx.time.timestep;
    }
    else
    {
        if (carbon_starv <= 0.0 && t_NPP <= 0.0)
            dr += 1.0 / ctx.time.timestep; // newIM 2021: carbon starvation occurs when the carbon stocks has been completly depleted, and carbon_starv represents t_carbon_storage (while it represents t_NPPneg when ctx.opt._LA_regulation==0)
    }
    // if the water availabiliy in the root zone is below the lethal level, the tree dies, !!!: not that deterministic, right?
    if (phi_root < (t_phi_lethal))
        dr += 1.0 / ctx.time.timestep;
    if (ctx.time.iter == int(ctx.time.nbiter - 1))
        ctx.out.output[26] << t_wsg << "\t" << basal << "\t" << dbh << "\t" << dr << "\n";

    /*if (ctx.time.iter>=622 && dr*ctx.time.timestep>=0.8) {
        cout<< "high deathrate: wsg=" << t_wsg << "; basal=" << basal << "; dbh="  << dbh << "; dr="  << dr*ctx.time.timestep   << "; carbon_starv="  << carbon_starv   << "; NPP="  << t_NPP   << "; phi_root="  << phi_root   << "; ctx.S[t_sp_lab].s_phi_lethal=" << ctx.S[t_sp_lab].s_phi_lethal << "; t_WSF=" << t_WSF << "; t_WSF_A=" << t_WSF_A << "; t_LA=" << t_LA << endl;
    }*/

    return dr * ctx.time.timestep;
}
#else
float Tree::DeathRate(Context &ctx, float dbh, float carbon_starv)
{
    float dr = 0.0;
    float basal = fmaxf(ctx.params.m - ctx.params.m1 * t_wsg, 0.0);

    dr = basal;
    if (ctx.opt._LA_regulation == 0)
    {
        if (carbon_starv > t_leaflifespan)
            dr += 1.0 / ctx.time.timestep;
    }
    else
    {
        if (carbon_starv <= 0.0 && t_NPP <= 0.0)
            dr += 1.0 / ctx.time.timestep; // newIM 2021: carbon starvation occurs when the carbon stocks has been completly depleted, and carbon_starv represents t_carbon_storage (while it represents t_NPPneg when ctx.opt._LA_regulation==0)
    }

    return dr * ctx.time.timestep;
}
#endif

// ###############################################################
//  Farquhar von Caemmerer Berry model -- called by  Tree::GPPleaf
// ###############################################################

#ifdef WATER
//! - Function Tree::GPPleaf when the WATER option is off, Tree::FluxesLeaf when WATER option is on
//! - Includes the Farquhar model of photosynthesis and Medlyn et al. 2011 model of stomatal conductance -- see also Prentice et al 2014 Ecology Letters, Lin et al 2015 Nature Climate Change, Dewar et al. 2018 New phytologist; min added in order to prevent ci:ca bigger than 1 (even though Ehleringer et al 1986 reported some values above 1 (Fig3)
//! - Tree::FluxesLeaf iterates on the leaf conditions (Tl, Cs and Dleaf) and calls another function Tree::Photosyn, which itself implements the models of photosynthsis and stomatal conductance.
//! - When WATER option is on, the stomatal conductance parameter (ctx.params.g1) and photosynthetic capacities are reduced under water stressed conditions: stomatal and non-stomatal responses to water stress (see Egea et al. 2011 AFM, Zhou et al. 2013 AFM, Zhou et al. 2014 PCE, etc.). This is done by multipliying the values of ctx.params.g1 and of Vcmax and Jmax in absence of water stress by water stress factors (WSF and WSF_A respectively).
//! - Returns the primary productivity (assimilation) per unit leaf area, i.e. in micromoles C ctx.params.m-2 s-1, and water fluxes (transpiration) in micromol H20 ctx.params.m-2 s-1.
//! - It is converted into gC per ctx.params.m^2 of leaf per ctx.time.timestep by "ctx.time.nbhours_covered*15.7788*ctx.time.timestep" where 15.7788 = 3600*365.25*12/1000000 and ctx.time.nbhours_covered is the duration read in the daily variation file
//! - NB1: 12 is the molar mass of carbon NB2: ctx.time.timestep is given as fraction of a year, so what is computed is actually the full assimilation per year which, in turn, is multiplied by the fraction per year that is under consideration.
//! - BEWARE: slight inconsistency through use of 365.25 when daily ctx.time.timestep is likely to be given as 365, but not incorrect.
//! - It uses lookup tables for acceleration of computation of T dependencies. cf. Bernacchi et al 2003 PCE; von Caemmerer 2000

//!  the function Tree::FluxesLeaf implements an iterative algorithm to determine the leaf temperature, and CO2 concentration and VPD at the leaf surface, ie the ones that should be taken as arguments of the Farquhar and stomatal conductance models. The latters are themeselves implemented within the function Tree::Photosyn (analogy with the MAESTRA/MAESPA code from which this iterative scheme has been inspired, see Medlyn et al. 2007 Tree physiology). At each iteration, Tree::FluxesLeaf computes carbon and water fluxes at the leaf level using the Farqhuar model and the Penman-Monteith equation, and updates T, CO2 and VPD accordingly. Starting by assuming that leaf T, CO2 and VPD equal the ones of the surrounding air, the iterative scheme stops when the difference in leaf T between two consecutive iterations is negligble. The corresponding carbon assimilation and water transpiration are then retained and provided as outputs of Tree::FluxesLeaf.
leafFluxes Tree::FluxesLeaf(Context &ctx, float PPFD, float VPDa, float Ta, float WIND, float ExtinctLW, float PPFDtop, float PPFDinc, float Tatop, float VPDatop)
{

    float TLEAF = 0.0, CS = 0.0, DS = 0.0,                              // leaf temperature and CO2 concentration and vapour pressure deficit at the leaf surface.
        ALEAF = 0.0, ET = 0.0,                                          // net assimilation rate and leaf-level evapotranspiration
        GH = 0.0, GBV = 0.0, GSV = 0.0, GV = 0.0, GSC = 0.0, GBH = 0.0, // conductance terms for Penman-Monteith calculation
        Rnetiso = 0.0, lambdaET = 0.0, HDIVT = 0.0;                     // terms used for Penman-Monteth calculation

    int ITERMAX = 30;
    // int itmax=ITERMAX;
    int convTA = int(ctx.lookup.iTaccuracy * Ta); // discrete temperature to avoid repeated computation
    // int convVPDA=int(ctx.lookup.iVPDaccuracy*VPDa); // discrete VPD to avoid repeated computation
    int convTAtop = int(ctx.lookup.iTaccuracy * Tatop);       // discrete temperature to avoid repeated computation
    int convVPDAtop = int(ctx.lookup.iVPDaccuracy * VPDatop); // discrete VPD to avoid repeated computation
    // float ESAT=LookUp_ESAT[convTA]; //ESAT: Saturation vapour pressure (in kPa)
    // VPD = ESAT*(1.0 - RH); // Calculation of air VPD from relative humidity and air temperature, in kPa (!!), needed only if RH, and not VPD, is provided as an argument.

    // Computation of Penman-Monteith terms that do not depend on TLEAF, CS and DS :

    // float RSOL = PPFD*ctx.params.PPFDtoSW;   // Incoming solar energy in J ctx.params.m-2 s-1 -- ctx.params.PPFDtoSW is to account for the fact that PAR actually account for ca. 50% of the total solar energy at the Earth's surface (ie after its way through the atmosphere), and Penman-Monteith should include all the shortwave energy (ie PAR + NIR), and convert from micromol quanta to J (ctx.params.PPFDtoSW~ 2.0/ 4.57).
    float RSOL = 0.2188 * PPFD + 0.0015 * PPFDtop + 0.0174 * PPFDinc; // this is based on the relationship between NIR/PAR ratio in the understory and PARtop/PAR or LAI, reported in Kume et al. 2011 Journal of Plant Research (Fig. 4 therein), to account for the absorbed NIR. 0.2188=1/4.47; 0.0015=0.1/(4.57*14.48) with 0.1~leaf absorptance in NIR and 14.48 fitted parameter in Kumme et al; 0.0174=0.1*11.52/(4.57*14.48), with 11.52 the other fitted parameter in Kumme et al., and 1/4.57 is to convert from microml ctx.params.m-2 s-1 to J.
    // if (PPFDinc/PPFDtop <0.8) {
    // cout << "PPFDtop=" << PPFDtop << " PPFD=" << PPFD <<" PPFDinc=" << PPFDinc << " RSOL=" << RSOL << " RSOL_NIR=" << RSOL_NIR << endl;
    // }
    float LHV = (H2OLV0 - 2.365e3 * Ta) * H2OMW;               // Latent heat of water vapour at air temperature (J mol-1) (this is the "landa" in PM equation)
    float SLOPE = ctx.lookup.LookUp_SLOPE[convTA];                        /* Slope of the dependence of saturated vapor pressure with temperature (Jones 2013, Equation (5.15), p.102)
                                                                This is constant s in Penman-Monteith equation  (Pa K-1) */
    float CMOLAR = 1000.0 * ctx.params.PRESS / (RCONST * (Ta - ABSZERO)); // 1000 because ctx.params.PRESS is in kPa
    float GAMMA = 1000.0 * ctx.params.PRESS * CPAIR * AIRMA / LHV;        // Psychrometric constant; 1000.0 to convert kPa into Pa. !! IM to be checked !! factor AIRMA not consistent with Jérôme's note and Appendix 3 of Jones, but similar to Duursma 's package
    float GRADN = ctx.lookup.LookUp_GRADN[convTA];                        // Radiation conductance (mol ctx.params.m-2 s-1) at air temperature (Jones 2013, Eq 5.10 p. 101). // IM to be double-checked with Jones (not fully consistent with Medlyn et al. 2007)
    float GBHU = 0.003 * sqrt(WIND / t_wleaf) * CMOLAR;        // Boundary ctx.diag.layer conductance to heat transfer by forced convection (single sided) in mol ctx.params.m-2 s-1; Leuning et al (1995) PC&E 18:1183-1200 Eqn E1; Equation A2 in Medlyn et al. 2007
    // Rnetiso = ctx.params.absorptance_leaves*RSOL - ctx.lookup.LookUp_INLR[convTA][convVPDA]; // Calculation of isothermal net radiation (J ctx.params.m-2 s-1; Jones (2013) equation (5.4) p.100)
    Rnetiso = RSOL - ctx.lookup.LookUp_INLR[convTAtop][convVPDAtop] * ExtinctLW; // Calculation of isothermal net radiation (J ctx.params.m-2 s-1; Jones (2013) equation (5.4) p.100). No need of ctx.params.absorptance_leaves, since PPFD provided in argument is already the absorbed flux.

    // if (ExtinctLW<0.2) {
    //    cout << "Rnetiso=" << Rnetiso << " RSOL=" << RSOL << " INLR=" << ctx.lookup.LookUp_INLR[convTA][convVPDA]  << " INLR_top=" << ctx.lookup.LookUp_INLR[convTAtop][convVPDAtop] << " //ExtinctLW=" << ExtinctLW <<  endl;
    // }

    // Iterative scheme to determine TLEAF and CO2 concentration and VPD at leaf surface, and the corresponding Anet and ET:

    TLEAF = Ta; // Initialization of leaf T, CO2 & VPD at leaf surface equal to the ones of surrounding air
    DS = VPDa;
    CS = ctx.params.Cair; // IM: check whether we should not make ctx.params.Cair vary across the canopy.
    leafFluxes ps;

    for (int ITER = 0; ITER < ITERMAX; ITER++)
    { // Iterative loop:

        if (CS < 0 || DS < 0 || CS > 100000)
        {
            cout << "Warning in FluxesLeaf ! DS=" << DS << "; CS=" << CS << "; PPFD=" << PPFD << "; TLEAF=" << TLEAF << endl;
        }

        ps = Photosyn(ctx, PPFD, TLEAF, CS, DS); // Compute CO2 assimilation rate and stomatal conductance at TLEAF, CS and DS
        ALEAF = ps.carbon_flux;             // in micromolCO2 ctx.params.m-2 s-1
        GSC = ps.water_flux;                // in mol CO2 ctx.params.m-2 s-1 by Photosyn function

        if (isnan(GSC) || isnan(ALEAF))
        {
            cout << "Warning in FluxesLeaf: " << " ALEAF= " << ALEAF << ", GSC= " << GSC << ", CS= " << CS << ", PPFD=" << PPFD << ", DS= " << DS << ", TLEAF= " << TLEAF << ", ITER=" << ITER << endl;
        }

        // Computation of the different conductance terms
        float GBHF; // Boundary ctx.diag.layer conductance to heat transfer by free convection (single sided) in mol ctx.params.m-2 s-1. Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4. Equation A1 in Medlyn et al. 2007
        if (abs(TLEAF - Ta) > 0.01)
        { // this is to speed up: calculating GBHF involves a pow, so do only if difftemp is not zero.
            // float GRASHOF = 1.6e8 * abs(TLEAF-TAIR) * pow(WLEAF,3.0); // Grashof number
            // GBHF = 0.5 * DHEAT * pow(GRASHOF,0.25) / WLEAF * CMOLAR;
            GBHF = 56.2341 * DHEAT * pow(abs(TLEAF - Ta) / t_wleaf, 0.25) * CMOLAR;
        }
        else
            GBHF = 0.0;
        GBH = GBHU + GBHF;      // Total boundary ctx.diag.layer conductance to heat (single sided), in mol ctx.params.m-2 s-1
        GH = 2 * (GBH + GRADN); // GH heat and radiative conductance (free & forced & radiative components), mol ctx.params.m-2 s-1, to be used in the isothermal version of the Penman-Monteith equation (always two sided, hence the factor 2)
        // if (t_sp_lab==6) GBH*=2; // test to see if assuming that Cecropia obtusa is amphistomatous (and not hypostomatous) would solve the issue found when ALEAF is big, wleaf is big (and hence GBH small) that leads to negative CS.... IM 09/09/2021
        HDIVT = GH * CPAIR * AIRMA; // in J ctx.params.m-2 s-1 K-1, useful in several calculations below
        GBV = GBVGBH * GBH;         // Boundary ctx.diag.layer conductance for water vapour, for hypostomatous leaves (single-sided value). Note: we thus here assumed that all leaves are hypostomatous, this is a reasonable assumption as few species have been reported to have amphistomatous leaves in tropical forests (e.g. Drake et al. 2019 New phytologist; Muir 2015 Proc Roy Soc; email with Lawren Sack Oct. 29th 2020). Cf. comment in Leuning et al. 1995 after equ. E5.
        GSV = GSVGSC * GSC;         // Stomatal conductance for water vapour, in mol H2O ctx.params.m-2 s-1
        if ((GSV / GBV) > 100000)
            GV = GBV;
        else
            GV = (GBV * GSV) / (GBV + GSV); // Total conductance to water vapour (single-sided value: hypostomatous; stomatal & bdry ctx.diag.layer components in series), mol H2O ctx.params.m-2 s-1 //

        if (isnan(GSV) || isnan(GV) || (GSV < 0))
        {
            cout << "Warning in FluxesLeaf: " << " GV= " << GV << ", GSC= " << GSC << ", GBV= " << GBV << ", GSV=" << GSV << ", GSC=" << GSC << ", ALEAF=" << ALEAF << ", GSVGSC= " << GSVGSC << ", ALEAF=" << ALEAF << ", PPFD=" << PPFD << ", DS=" << DS << ", ITER=" << ITER << ", WIND=" << WIND << endl;
        }

        // Computation of leaf-level evapotranspiration (ET) using the Penman-Monteith equation (mol H2O ctx.params.m-2 s-1).
        if (GV > 0.0)
            lambdaET = (SLOPE * Rnetiso + 1000.0 * VPDa * HDIVT) / (SLOPE + GAMMA * GH / GV);
        else
            lambdaET = 0.0;

        // if (lambdaET<0) lambdaET = 0.0; // check whether this is ok. remy says keep negative value but force ET TO ZERO
        if (lambdaET < 0)
            ET = 0.0;
        else
            ET = lambdaET / LHV; // in mol H2O ctx.params.m-2 s-1

        // Update of TLEAF, CS, and DS
        CS = ctx.params.Cair - ALEAF * GBHGBC / GBH; // CO2 concentration at the leaf surface, in micromol mol-1
        if (CS < 0.0)
            CS = 1.0;
        DS = ET * ctx.params.PRESS / GSV; // in kPa; IM important change in regards to Duursma R package or MAESTRA/MAESPA code: GV replaced by GSV (cf for eg Medvigy et al. 2009 equ B16, cf. doc comments on Physiol V6)
        // DS = ET * ctx.params.PRESS / GV; // test cluster 14/05/22

        if (CS <= 0)
        {
            cout << "Warning in FluxesLeaf: " << " CS= " << CS << ", Cair=" << ctx.params.Cair << ", ALEAF=" << ALEAF << ", GBH=" << GBH << ", GBHGBC=" << GBHGBC << endl;
        }

        /*if(DS < 0 ) {

            cout << "Warning in FluxesLeaf: " << " DS= " << DS << ", VPDa= "<< VPDa << ", ITER=" << ITER << ", GSV=" << GSV << ", GSC=" << GSC << ", ALEAF=" << ALEAF << ", ET= " << ET << ", GV=" << GV << ", initial lambdaET=" << (SLOPE*Rnetiso+1000.0*VPDa*HDIVT)/(SLOPE+GAMMA*GH/GV) << ", Rnetiso=" << Rnetiso << "; PPFD=" << PPFD << "; Vcmax=" << t_Vcmax << "; Jmax=" << t_Jmax << "; t_Pmass=" << t_Pmass << "; ctx.S[spp].s_Pmass=" << ctx.S[t_sp_lab].s_Pmass << endl ;
        }*/

        // float TLEAF1 = Ta + (Rnetiso - lambdaET)/HDIVT; // Leaf temperature inferred from Jones (2013) equation (9.7), p.225. Note : in Vezy et al. 2018 they accounted for the rate of evaporation of water at the leaf sufarce when the leaf is wet in addition to the one of tranpiration. To do so they use Penman-Monteith with an "infinitly" large stomatal conductance to water vapor and use a weighted mean of both PM values...

        float TDIFF = (Rnetiso - lambdaET) / HDIVT;
        float TLEAF1 = Ta + 0.25 * TDIFF; // divide by 4 to slow down convergence and avoid big changes -- taken from code from Vezy et al. 2018 -- cf. discussion with Guerric Le Maire
        if (abs(TLEAF - TLEAF1) < 0.01)
        {
            // itmax=ITER;
            ITER = ITERMAX;

        } // If Tleaf diff < 0.01°C, stop search // check whether this threshold value is optimal
        TLEAF = TLEAF1; // update TLEAF for next iteration // cf. Vezy et al. and discussion with Guerric Lemaire: suggest to make smaller steps at each iteration to ease convergence (TDIFF -> TDIFF/4)

        if (isnan(DS) || isnan(ET) || isnan(ALEAF) || ET < 0.0 || TLEAF <= 0.0 || TLEAF > 50.0 || CS < 0 || DS < 0 || CS > 100000)
        {
            cout << "Warning in FluxesLeaf: " << " ALEAF= " << ALEAF << ", ET= " << ET << ", CS= " << CS << ", PPFD=" << PPFD << ", DS= " << DS << ", TLEAF= " << TLEAF << ", GSV=" << GSV << ", GV= " << GV << ", Rnetsiso=" << Rnetiso << "; lambdaET=" << lambdaET << "; HDIVT=" << HDIVT << "; TDIFF=" << TDIFF << "; Ta=" << Ta << "; Cair=" << ctx.params.Cair << "; GBH=" << GBH << "; Cair=" << ctx.params.Cair << "; GBHGBC=" << GBHGBC << "; t_site=" << t_site << "; t_sp_lab=" << t_sp_lab << "; ITER=" << ITER << endl;
        }
    }

    // cout << "itmax=" << itmax << endl;

    // At the end of loop, return water (ET, in umol H2O ctx.params.m-2 s-1) and carbon (assimilation) fluxes
    ET = ET * 1e6;

    // if (ctx.time.iter > 4000 && t_height> 10.0) {
    //     cout << " species:" << ctx.S[t_sp_lab].s_name << " Tree height=" << t_height << " ET=" << ET << " SLOPE=" << SLOPE << " Rnetiso=" << Rnetiso << " SLOPE*Rnetiso=" << SLOPE*Rnetiso <<  " 1000*VPDa*HDIVT=" << 1000.0*VPDa*HDIVT << " VPDa=" << VPDa << " HDIVT=" << HDIVT << " GAMMA*GH/GV=" << GAMMA*GH/GV << endl;
    // }

    leafFluxes outFluxes;
    outFluxes.carbon_flux = ALEAF; // in micromol CO2 ctx.params.m-2 s-1
    outFluxes.water_flux = ET;     // in micromol H20 ctx.params.m-2 s-1
    return outFluxes;
}

#ifdef G0

leafFluxes Tree::Photosyn(Context &ctx, float PPFD, float TLEAF, float CS, float DS)
{

    float Anet;

    // Parameters for the Farquhar model
    int convT = int(ctx.lookup.iTaccuracy * TLEAF);
    float KmT = ctx.lookup.LookUp_KmT[convT];                        // with temperature dependencies
    float GammaT = ctx.lookup.LookUp_GammaT[convT];                  // with temperature dependencies
    float Rday = t_Rdark * ctx.lookup.LookUp_Rleaf[convT] * DAYRESP; // leaf respiration with temperature dependencies and inhibition due to light (DAYRESP).

    ////////////// Model of stomatal conductance for CO2 /////////////
    float GSDIVA;
    if (DS < 0.0005)
        GSDIVA = (1.0 + t_g1 / sqrt(0.0005)) / CS; // VPDMIN in MAESPA/Plantecophys R package
    else
        GSDIVA = (1.0 + t_g1 / sqrt(DS)) / CS; // Medlyn et al. 2011 Global Change Biology. Note that ctx.params.g1 is now a tree class variable that is initiated at tree birth (in Tree::Birth) and updated in Tree::Water_Availability depending on the tree water stress.
    // float GSDIVA = (1.0 + ctx.params.g1/sqrt(VPD))/(CS-GammaT); //Dewar et al. 2018 New Phytologist
    if (GSDIVA <= 0 || isnan(GSDIVA))
    {
        cout << "Warning in Photosyn ! GSDIVA=" << GSDIVA << "; DS=" << DS << "; CS=" << CS << endl;
        cout << "In Photosyn: PPFD=" << PPFD << "; TLEAF=" << TLEAF << "; CS=" << CS << "; DS=" << DS << endl;
    }

    //////////// Farquhar model of photosynthesis //////////

    float A, B, C, AC, AJ;

    // Solution when Rubisco activity is limiting
    float VcmaxTW = t_Vcmax * ctx.lookup.LookUp_VcmaxT[convT] * t_WSF_A;
    A = ctx.params.g0 + GSDIVA * (VcmaxTW - Rday);
    B = (1.0 - CS * GSDIVA) * (VcmaxTW - Rday) + ctx.params.g0 * (KmT - CS) - GSDIVA * (VcmaxTW * GammaT + KmT * Rday);
    C = -(1.0 - CS * GSDIVA) * (VcmaxTW * GammaT + KmT * Rday) - ctx.params.g0 * KmT * CS;
    float CIC = QUAD(A, B, C, 1); // Finds larger root of quadratic function
    if ((CIC <= 0.0) || (CIC - CS > 0.01))
    { // newIM: have replaced CIC>CS by (CIC-CS>0.01) to avoid error due to comparing floats (maybe there would be a best way to fix this...)
        AC = 0.0;
    }
    else
        AC = VcmaxTW * (CIC - GammaT) / (CIC + KmT);

    // Solution when electron transport rate is limiting
    float JmaxTW = t_Jmax * ctx.lookup.LookUp_JmaxT[convT] * t_WSF_A;
    float I = ctx.params.alpha * PPFD;
    float J = QUAD(ctx.params.theta, -(I + JmaxTW), JmaxTW * I, -1); // ctx.params.theta is the convexity term for electron transport rates (dimensionless, 0–1), here provided in input (default value =0.7, as in von Caemmerer 2000, but a value of 0.85 is used in Duursma' R package
    float VJ = J * 0.25;
    A = ctx.params.g0 + GSDIVA * (VJ - Rday);
    B = (1.0 - CS * GSDIVA) * (VJ - Rday) + ctx.params.g0 * (2.0 * GammaT - CS) - GSDIVA * (VJ * GammaT + 2.0 * GammaT * Rday);
    C = -(1.0 - CS * GSDIVA) * GammaT * (VJ + 2.0 * Rday) - ctx.params.g0 * 2.0 * GammaT * CS;
    float CIJ = QUAD(A, B, C, 1); // Finds larger root of quadratic function
    AJ = VJ * (CIJ - GammaT) / (CIC + 2.0 * GammaT);
    if ((AJ - Rday) < 1.0e-6)
    {
        CIJ = CS; // Below light compensation point
        AJ = VJ * (CIJ - GammaT) / (CIJ + 2.0 * GammaT);
    }

    // Co-limitation
    Anet = fminf(AC, AJ) - Rday; // compute the *net* assimilation rate, in micromol C02 ctx.params.m-2 s-1

    if (isnan(Anet) || isnan(GSDIVA))
    {
        cout << "Warning in Photosyn: " << ", Anet= " << Anet << ", GSDIVA=" << GSDIVA << "; Anet*GSDIVA=" << Anet * GSDIVA << ", AC= " << AC << ", AJ= " << AJ << ", Rday=" << Rday << ", CS=" << CS << "; CIC=" << CIC << "; CIJ=" << CIJ << ", DS= " << DS << ", TLEAF= " << TLEAF << ", t_WSF_A=" << t_WSF_A << ", t_WSF= " << t_WSF << "; GammaT=" << GammaT << "; KmT=" << KmT << "; J=" << J << "; PPFD=" << PPFD << "; T_Rdark=" << t_Rdark << "; convT=" << convT << endl;
    }

    if (isnan(Anet * GSDIVA) && Anet == 0.0)
    { // to solve nan issue when Anet=0 and GSDIVA=inf
        Anet = 1.0e-5;
    }

    leafFluxes outPhotosyn;
    outPhotosyn.carbon_flux = Anet;              // in micromol C02 ctx.params.m-2 s-1
    outPhotosyn.water_flux = ctx.params.g0 + GSDIVA * Anet; // in mol CO2 ctx.params.m-2 s-1, the conversion to mol H20 is made in function FluxesLeaf.
    if (outPhotosyn.water_flux < ctx.params.g0)
        outPhotosyn.water_flux = ctx.params.g0;

    if (isnan(outPhotosyn.water_flux) || isnan(outPhotosyn.carbon_flux))
    {
        cout << "Warning in Photosyn: " << ", Anet= " << Anet << ", AC= " << AC << ", AJ= " << AJ << ", Rday=" << Rday << ", GSDIVA=" << GSDIVA << ", CS=" << CS << "; CIC=" << CIC << "; CIJ=" << CIJ << ", DS= " << DS << ", TLEAF= " << TLEAF << ", t_WSF_A=" << t_WSF_A << ", t_WSF= " << t_WSF << "; GammaT=" << GammaT << "; KmT=" << KmT << "; J=" << J << "; PPFD=" << PPFD << "; t_Rdark=" << t_Rdark << "; convT=" << convT << endl;
    }

    return outPhotosyn;
}

#else

leafFluxes Tree::Photosyn(Context &ctx, float PPFD, float TLEAF, float CS, float DS)
{

    float CI, Anet;

    // Parameters for the Farquhar model
    int convT = int(ctx.lookup.iTaccuracy * TLEAF);
    float KmT = ctx.lookup.LookUp_KmT[convT];                        // with temperature dependencies
    float GammaT = ctx.lookup.LookUp_GammaT[convT];                  // with temperature dependencies
    float Rday = t_Rdark * ctx.lookup.LookUp_Rleaf[convT] * DAYRESP; // leaf respiration with temperature dependencies and inhibition due to light (DAYRESP).

    // Model of stomatal conductance for CO2
    float GSDIVA = (1.0 + t_g1 / sqrt(DS)) / CS; // Medlyn et al. 2011 Global Change Biology. Note that ctx.params.g1 is now a tree class variable that is initiated at tree birth (in Tree::Birth) and updated in Tree::Water_Availability depending on the tree water stress.
    // float GSDIVA = (1.0 + ctx.params.g1/sqrt(VPD))/(CS-GammaT); //Dewar et al. 2018 New Phytologist
    if (GSDIVA <= 0 || isnan(GSDIVA))
    {
        cout << "Warning in Photosyn ! GSDIVA=" << GSDIVA << "; DS=" << DS << "; CS=" << CS << endl;
        cout << "In Photosyn: PPFD=" << PPFD << "; TLEAF=" << TLEAF << "; CS=" << CS << "; DS=" << DS << endl;
    }

    // Internal CO2 concentration, assuming that minimal conductance (cf. Duursma et al. 2019 New Phytologist, ~cuticular conductance) is negligeable (ctx.params.g0~0) [Note that, if ctx.params.g0>0, then a quadratic solution must be computed, which complexifies and should lengthen the speed of the code ==> to try later]. Expression derived using Medlyn et al. 2011 GCB model (see equ. 7 in Marechaux & Chave 2017, and equ. 13 in Medlyn et al. 2011 corrigendum)
    CI = CS * t_g1 / (t_g1 + sqrt(DS));

    // Farquhar model of photosynthesis
    // Solution when Rubisco activity is limiting
    float AC, AJ;
    float VcmaxTW = t_Vcmax * ctx.lookup.LookUp_VcmaxT[convT] * t_WSF_A;
    if ((CI <= 0.0) || (CI - CS > 0.01))
    { // newIM: have replaced CI>CS by (CI-CS>0.01) to avoid error due to comparing floats (maybe there would be a best way to fix this...)
        AC = 0.0;
    }
    else
        AC = VcmaxTW * (CI - GammaT) / (CI + KmT);
    // Solution when electron transport rate is limiting
    float JmaxTW = t_Jmax * ctx.lookup.LookUp_JmaxT[convT] * t_WSF_A;
    float I = ctx.params.alpha * PPFD;
    float J = (I + JmaxTW - sqrt((JmaxTW + I) * (JmaxTW + I) - 4.0 * ctx.params.theta * JmaxTW * I)) * 0.5 / ctx.params.theta; // THETA is the convexity term for electron transport rates (dimensionless, 0–1), here provided in input (default value =0.7, as in von Caemmerer 2000, but a value of 0.85 is used in Duursma' R package
    AJ = 0.25 * J * (CI - GammaT) / (CI + 2.0 * GammaT);
    if ((AJ - Rday) < 1.0e-6)
    {
        CI = CS; // Below light compensation point
        AJ = 0.25 * J * (CI - GammaT) / (CI + 2.0 * GammaT);
    }
    Anet = fminf(AC, AJ) - Rday; // compute the *net* assimilation rate, in micromol C02 ctx.params.m-2 s-1

    if (isnan(Anet) || isnan(GSDIVA))
    {
        cout << "Warning in Photosyn: " << ", Anet= " << Anet << ", GSDIVA=" << GSDIVA << "; Anet*GSDIVA=" << Anet * GSDIVA << ", AC= " << AC << ", AJ= " << AJ << ", Rday=" << Rday << ", CS=" << CS << "; CI=" << CI << ", DS= " << DS << ", TLEAF= " << TLEAF << ", t_WSF_A=" << t_WSF_A << ", t_WSF= " << t_WSF << "; GammaT=" << GammaT << "; KmT=" << KmT << "; J=" << J << "; PPFD=" << PPFD << "; T_Rdark=" << t_Rdark << "; convT=" << convT << endl;
    }

    if (isnan(Anet * GSDIVA) && Anet == 0.0)
    { // to solve nan issue when ANet=0 and GSDIVA=inf
        Anet = 1.0e-5;
    }

    leafFluxes outPhotosyn;
    outPhotosyn.carbon_flux = Anet; // in micromol C02 ctx.params.m-2 s-1
    if (fabs(GSDIVA * Anet) < 1.0e-6)
        outPhotosyn.water_flux = 1.0e-5; // GSC should not have null values
    else
        outPhotosyn.water_flux = fabs(GSDIVA * Anet); // in mol CO2 ctx.params.m-2 s-1, the conversion to mol H20 is made in function FluxesLeaf. Note that it is here assumed that the minimal stomatal conductance (ctx.params.g0) is negligeable.
    // if (outPhotosyn.water_flux<ctx.params.g0) outPhotosyn.water_flux=ctx.params.g0; // just as a reminder: in both Duursma's pacakage and Vezy's MAESPA code, such line is added (in Duursma, it is written: "# Extra step here; GS can be negative", in Vezy, it is written: "! Set nearly zero conductance (for numerical reasons)". this results in no possible negative values or zero values for GSC.

    if (isnan(outPhotosyn.water_flux) || isnan(outPhotosyn.carbon_flux))
    {
        cout << "Warning in Photosyn: " << ", Anet= " << Anet << ", AC= " << AC << ", AJ= " << AJ << ", Rday=" << Rday << ", GSDIVA=" << GSDIVA << ", CS=" << CS << "; CI=" << CI << ", DS= " << DS << ", TLEAF= " << TLEAF << ", t_WSF_A=" << t_WSF_A << ", t_WSF= " << t_WSF << "; GammaT=" << GammaT << "; KmT=" << KmT << "; J=" << J << "; PPFD=" << PPFD << "; T_Rdark=" << t_Rdark << "; convT=" << convT << endl;
    }

    return outPhotosyn;
}

#endif

//! - Function Tree::dailyGPPleaf when the WATER option is off, Tree::dailyFluxesLeaf when WATER option is on
//! - Returns the *daily* primary productivity (assimilation; computed from Tree::GPPleaf, or Tree::FluxesLeaf) averaged across the daily fluctuations in climatic conditions, per unit leaf area, in micromoles C/ctx.params.m^2/s.
leafFluxes Tree::dailyFluxesLeaf(Context &ctx, float PPFD, float VPD, float T, float W, float ExtinctLW, float PPFDinc, float &wf_1016)
{
    float dailyA = 0.0;
    float dailylT = 0.0;

    // cout << "In dailyFluxesLeaf PPFD=" << PPFD << "; T=" << T << "; VPD=" << VPD << "; W=" << W << endl;

    for (int i = 0; i < ctx.time.nbsteps_varday; i++)
    {

#ifdef FULL_CLIMATE

        float ppfd_vardaytimestep = 0.0;
        float ppfd_top_vardaytimestep = 0.0;
        float t_top_vardaytimestep = 0.0;
        float vpd_top_vardaytimestep = 0.0;
        float ppfd_inc_vardaytimestep = 0.0;
        float vpd_vardaytimestep = 0.0;
        float t_vardaytimestep = 0.0;
        float ws_vardaytimestep = 0.0;

        if (ctx.time.iter == -1)
        {
            ppfd_vardaytimestep = PPFD * ctx.climate.WDailyMean_all[i] * ctx.params.SWtoPPFD;
            ppfd_top_vardaytimestep = ctx.climate.WDailyMean_all[i] * ctx.params.SWtoPPFD;
            t_top_vardaytimestep = ctx.climate.tDailyMean_all[i];
            vpd_top_vardaytimestep = ctx.climate.VPDDailyMean_all[i];
            ppfd_inc_vardaytimestep = ctx.climate.WDailyMean_all[i] * ctx.params.SWtoPPFD;
            vpd_vardaytimestep = VPD * ctx.climate.VPDDailyMean_all[i];
            t_vardaytimestep = ctx.climate.tDailyMean_all[i] - T;
            ws_vardaytimestep = W * ctx.climate.windDailyMean_all[i];
        }
        else
        {
            int c = (ctx.time.iter % ctx.time.nbdays) * ctx.time.nbsteps_varday + i;
            ppfd_vardaytimestep = PPFD * ctx.climate.varday_light[c] * ctx.params.SWtoPPFD;
            ppfd_top_vardaytimestep = ctx.climate.varday_light[c] * ctx.params.SWtoPPFD;
            t_top_vardaytimestep = ctx.climate.varday_T[c];
            vpd_top_vardaytimestep = ctx.climate.varday_vpd[c];
            ppfd_inc_vardaytimestep = PPFDinc * ctx.climate.varday_light[c] * ctx.params.SWtoPPFD;
            vpd_vardaytimestep = VPD * ctx.climate.varday_vpd[c];
            t_vardaytimestep = ctx.climate.varday_T[c] - T;
            ws_vardaytimestep = W * ctx.climate.varday_WS[c];
        }

        if (vpd_vardaytimestep <= 0)
        {
            cout << "Warning in dailyFluxesLeaf !!! vpd_wardaytimestep <=0" << endl;
        }
        if (ppfd_vardaytimestep <= 0)
        {
            cout << "Warning ppfd_vardaytimestep  <=0 in dailyFluxesLeaf !!!" << endl;
        }

#else

        if (ctx.time.iter == -1)
        {
            float ppfd_vardaytimestep = PPFD * ctx.climate.varday_light[i];
            float ppfd_top_vardaytimestep = ctx.climate.WDailyMean_year * ctx.climate.varday_light[i];
            float t_top_vardaytimestep = ctx.climate.tDailyMean_year * ctx.climate.varday_T[i];
            float vpd_top_vardaytimestep = ctx.climate.VPDDailyMean_year * ctx.climate.varday_vpd[i];
            float ppfd_inc_vardaytimestep = PPFDinc * ctx.climate.varday_light[i];
            float vpd_vardaytimestep = VPD * ctx.climate.varday_vpd[i];
            float t_vardaytimestep = T * ctx.climate.varday_T[i];
            float ws_vardaytimestep = W * ctx.climate.varday_WS[i];
        }
        else
        {
            float ppfd_vardaytimestep = PPFD * ctx.climate.varday_light[i];
            float ppfd_top_vardaytimestep = ctx.climate.WDailyMean * ctx.climate.varday_light[i];
            float t_top_vardaytimestep = ctx.climate.tDailyMean * ctx.climate.varday_T[i];
            float vpd_top_vardaytimestep = ctx.climate.VPDDailyMean * ctx.climate.varday_vpd[i];
            float ppfd_inc_vardaytimestep = PPFDinc * ctx.climate.varday_light[i];
            float vpd_vardaytimestep = VPD * ctx.climate.varday_vpd[i];
            float t_vardaytimestep = T * ctx.climate.varday_T[i];
            float ws_vardaytimestep = W * ctx.climate.varday_WS[i];
        }

        if (vpd_vardaytimestep <= 0)
        {
            cout << "Warning in dailyFluxesLeaf !!! i=" << i << "; vpd_vardaytimestep=" << vpd_vardaytimestep << "; VPD=" << VPD << "; ctx.climate.varday_vpd[i]=" << ctx.climate.varday_vpd[i] << "; ppfd_vardaytimestep=" << ppfd_vardaytimestep << "; PPFD=" << PPFD << "; ctx.climate.varday_light[i]=" << ctx.climate.varday_light[i] << endl;
        }
        if (ppfd_vardaytimestep <= 0)
        {
            cout << "Warning ppfd_vardaytimestep  <=0 in dailyFluxesLeaf !!! ppfd_vardaytimestep=" << ppfd_vardaytimestep << "; PPFD=" << PPFD << "; ctx.climate.varday_light[i]=" << ctx.climate.varday_light[i] << endl;
        }

#endif

        leafFluxes fluxes_vardaytimestep = FluxesLeaf(ctx, ppfd_vardaytimestep, vpd_vardaytimestep, t_vardaytimestep, ws_vardaytimestep, ExtinctLW, ppfd_top_vardaytimestep, ppfd_inc_vardaytimestep, t_top_vardaytimestep, vpd_top_vardaytimestep);

        dailyA += fluxes_vardaytimestep.carbon_flux;
        dailylT += fluxes_vardaytimestep.water_flux;
        if (i > 5 && i < 18.5)
        {
            wf_1016 += fluxes_vardaytimestep.water_flux;
        }

        if (isnan(fluxes_vardaytimestep.carbon_flux) || isnan(fluxes_vardaytimestep.water_flux) || isnan(dailyA) || isnan(dailylT))
        {
            cout << "Warning in dailyFluxesLeaf: " << " fluxes_vardaytimestep.carbon_flux= " << fluxes_vardaytimestep.carbon_flux << ", fluxes_vardaytimestep.water_flux= " << fluxes_vardaytimestep.water_flux << ", dailyA=" << dailyA << ", dailylT=" << dailylT << ", i=" << i << endl;
        }

        // cout << dailyA << endl;
        //  deprecated in v.2.4.1: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
        //  Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
    }
    dailyA *= ctx.time.inv_nbsteps_varday;
    dailylT *= ctx.time.inv_nbsteps_varday;

    // cout << "In dailyFluxesLeaf, final daily average: dailyA=" << dailyA << "; dailylT=" << dailylT << endl;

    leafFluxes dailyF;
    dailyF.carbon_flux = dailyA; // in micromol C02 ctx.params.m-2 s-1
    dailyF.water_flux = dailylT; // in micromol H20 ctx.params.m-2 s-1

    if (isnan(dailyA) || isnan(dailylT))
    {
        cout << "Warning in dailyFluxesLeaf:" << " dailyA=" << dailyA << ", dailylT=" << dailylT << ", PPFD=" << PPFD << endl;
    }

    return dailyF;
}

#else

//! - Function Tree::GPPleaf when the WATER option is off, Tree::FluxesLeaf when WATER option is on
//! - Includes the Farquhar model of photosynthesis and Medlyn et al. 2011 model of stomatal conductance: formula for s_fci (ci/ca) -- see also Prentice et al 2014 Ecology Letters and Lin et al 2015 Nature Climate Change; min added in order to prevent ci:ca bigger than 1 (even though Ehleringer et al 1986 reported some values above 1 (Fig3)
//! - Returns the primary productivity (assimilation) per unit leaf area, i.e. in micromoles C/ctx.params.m^2/s.
//! - It is converted into gC per ctx.params.m^2 of leaf per ctx.time.timestep by "ctx.time.nbhours_covered*15.7788*ctx.time.timestep" where 15.7788 = 3600*365.25*12/1000000 and ctx.time.nbhours_covered is the duration read in the daily variation file
//! - NB1: 12 is the molar mass of carbon NB2: ctx.time.timestep is given as fraction of a year, so what is computed is actually the full assimilation per year which, in turn, is multiplied by the fraction per year that is under consideration.
//! - BEWARE: slight inconsistency through use of 365.25 when daily ctx.time.timestep is likely to be given as 365, but not incorrect. Commented version below was in use prior to version 2.3.0
//! lookup tables for acceleration of T dependence. cf. Bernacchi et al 2003 PCE; von Caemmerer 2000
float Tree::GPPleaf(Context &ctx, float PPFD, float VPD, float T)
{
    // v.2.3.0: ctx.params.theta defined as a global variable
    // ctx.params.theta=0.7;   // this is the fixed value of ctx.params.theta used by von Caemmerer 2000

    // float ctx.params.theta=0.76+0.018*T-0.00037*T*T;         // ctx.params.theta, but temperature dependent cf. Bernacchi et al 2003 PCE

    // Parameters for Farquhar model, with temperature dependencies
    int convT = int(ctx.lookup.iTaccuracy * T); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---

    // if(convT>500 || isnan(convT) || convT <0) cout << t_site << " | convT: " << convT << " | T: " << T << " | PPFD: " << PPFD << " | VPD: " << VPD << endl;
    float KmT = ctx.lookup.LookUp_KmT[convT];
    float GammaT = ctx.lookup.LookUp_GammaT[convT];

    // float ctx.params.g1 = -3.97 * t_wsg + 6.53 (Lin et al. 2015)

    float t_fci = ctx.params.g1 / (ctx.params.g1 + sqrt(VPD));
    float VcmaxT = t_Vcmax * ctx.lookup.LookUp_VcmaxT[convT];
    float JmaxT = t_Jmax * ctx.lookup.LookUp_JmaxT[convT];

    // Farquhar - -von Caemmerer - Berry model of carbon assimilation rate
    float I = ctx.params.alpha * PPFD;
    float J = (I + JmaxT - sqrt((JmaxT + I) * (JmaxT + I) - 4.0 * ctx.params.theta * JmaxT * I)) * 0.5 / ctx.params.theta;
    float A = fminf(VcmaxT / (t_fci + KmT), 0.25 * J / (t_fci + 2.0 * GammaT)) * (t_fci - GammaT);

    return A;
}

//! - Function Tree::dailyGPPleaf when the WATER option is off, Tree::dailyFluxesLeaf when WATER option is on
//! - Includes Medlyn et al. 2011 model of stomatal conductance: formula for s_fci (ci/ca) -- see also Prentice et al 2014 Ecology Letters and Lin et al 2015 Nature Climate Change; min added in order to prevent ci:ca bigger than 1 (even though Ehleringer et al 1986 reported some values above 1 (Fig3)
//! - Returns the *daily* primary productivity (assimilation; computed from Tree::GPPleaf)averaged across the daily fluctuations in climatic conditions, per unit leaf area, in micromoles C/ctx.params.m^2/s.
float Tree::dailyGPPleaf(Context &ctx, float PPFD, float VPD, float T)
{
    float dailyA = 0.0;

    for (int i = 0; i < ctx.time.nbsteps_varday; i++)
    {
        // cout << t_site << " i: " << i << " tempRday: " << tempRday << endl;
        float ppfd_vardaytimestep = PPFD * ctx.climate.varday_light[i];
        float vpd_vardaytimestep = VPD * ctx.climate.varday_vpd[i];
        float t_vardaytimestep = T * ctx.climate.varday_T[i];
        if (ppfd_vardaytimestep > 0.1)
            dailyA += Tree::GPPleaf(ctx, ppfd_vardaytimestep, vpd_vardaytimestep, t_vardaytimestep);
        // else { cout << endl << t_site << " species: " << t_s->s_name << " t_age: " << t_age << " PPFD: " << ppfd_vardaytimestep << " vpd_vardaytimestep " << vpd_vardaytimestep << " t_vardaytimestep: " << t_vardaytimestep << " GPPleaf: " << Tree::GPPleaf(ctx, ppfd_vardaytimestep,vpd_vardaytimestep,t_vardaytimestep) << endl;}
        //  deprecated in v.2.4.1: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
        //  Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)

        // the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        // float hhA=0;
        // hhA=GPPleaf(PPFD*vardaytime_light[i], VPD*vardaytime_vpd[i], T*vardaytime_T[i]);
        // float ctx.params.alpha=ctx.params.phi*PPFD*vardaytime_light[i]/hhA;
        // float D=ctx.params.klight*ctx.params.dens*CD;
        // hhA*=ctx.params.alpha/(D*(ctx.params.alpha-1))*log(ctx.params.alpha/(1+(ctx.params.alpha-1)*exp(-D)));
        // dailyA+=hhA;
    }
    dailyA *= ctx.time.inv_nbsteps_varday;
    return dailyA;
}

//! - Function Tree::dailyGPPcrown similar to Tree::dailyGPPleaf
//! Faster, whole crown GPP calculation
float Tree::dailyGPPcrown(Context &ctx, float PPFD, float VPD, float T, float LAI)
{
    float ppfde, dailyA = 0.0;

    for (int i = 0; i < ctx.time.nbsteps_varday; i++)
    {
        ppfde = PPFD * ctx.climate.varday_light[i];
        if (ppfde > 0.1)
            // new v.2.3.0: compute GPP only if enough light is available threshold is arbitrary, but set to be low: in full sunlight ppfd is aroung 700 W/m2, and even at dawn, it is ca 3% of the max value, or 20 W/m2. The minimum threshold is set to 0.1 W/m2
            // Future update: compute slightly more efficiently, using 3-hourly values? This will have to be aligned with climate forcing layers (e.g. NCAR)
            dailyA += Tree::GPPleaf(ctx, ppfde, VPD * ctx.climate.varday_vpd[i], T * ctx.climate.varday_T[i]);
        // the 6 lines in comment below corresponds to a finer version in which the multiplier is computed and used every 48 half hour, ie. with the corresponding environment instead of assuming a constant multiplier correponding the one at maximum incoming irradiance
        // float hhA=0;
        // hhA=GPPleaf(PPFD*vardaytime_light[i], VPD*vardaytime_vpd[i], T*vardaytime_T[i]);
        // float ctx.params.alpha=ctx.params.phi*PPFD*vardaytime_light[i]/hhA;
        // float D=ctx.params.klight*LAI;
        // hhA*=ctx.params.alpha/(D*(ctx.params.alpha-1))*log(ctx.params.alpha/(1+(ctx.params.alpha-1)*exp(-D)));
        // dailyA+=hhA;
    }
    // vardaytime_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for vardaytime_vpd and vardaytime_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model

    float alpha_crown = ctx.params.phi * PPFD / GPPleaf(ctx, PPFD, VPD, T); // RENAMED: ctx.params.alpha → alpha_crown (FvCB canopy integration ratio; distinct from global ctx.params.alpha = apparent quantum yield)
    float D = ctx.params.klight * LAI;                                                         // D is a non-dimensional figure used to compute the multiplier below, update in v.2.5: replaced ctx.params.dens * CD by LAI
    dailyA *= alpha_crown / (D * (alpha_crown - 1)) * log(alpha_crown / (1 + (alpha_crown - 1) * exp(-D))); // the FvCB assimilation rate computed at the top of the tree crown is multiplied by a multiplier<1, to account for the lower rate at lower light level within the crown depth. This multiplier is computed assuming that change in photosynthetic assimilation rate within a tree crown is mainly due to light decrease due to self-shading following a Michealis-menten relationship (ie. we assume that 1/ the change is not due to changes in VPD or temperature, which are supposed homogeneous at the intra-crown scale, and 2/ that other tree contributions to light decrease is neglected).

    dailyA *= ctx.time.inv_nbsteps_varday;
    return dailyA;
}

#endif

//! NEW in v. 2.4.0: Separate function for calculation of Rday
//! separation necessitates extra function calls, but allows for decoupling of photosynthesis and leaf respiration in special cases (i.e. no photosynthesis at very low ppfd, but still respiration
//! at low light levels, should leaf respiration not be closer to dark respiration than to respiration in full sunlight, i.e. should we not lose the factor 0.4 then?

#ifdef WATER
// Function Rdayleaf has not been modified when WATER is defined: there is no clear consensus on the effect of water shortage on respiration and how to implement such potential effect in models, but I have to further check this and it could be modified. IM 4 dec 2019.
// Rdayleaf is actually not used anymore when WATER is activated
#endif

float Tree::Rdayleaf(Context &ctx, float T)
{
    int convT = int(ctx.lookup.iTaccuracy * T);
    // if(T < 0 || isnan(T) || T > 50) cout << t_site << " species: " << t_s->s_name << " convT: " << T << endl;
    float Rday_leaf = t_Rdark * ctx.lookup.LookUp_Rleaf[convT]; // new IM: no redundancy anymore between LookUp_Rday and LookUp_Rnight
    return Rday_leaf;
}

float Tree::dailyRdayleaf(Context &ctx, float T)
{
    float Rdayleaf_daily = 0.0;
    for (int i = 0; i < ctx.time.nbsteps_varday; i++)
        Rdayleaf_daily += Tree::Rdayleaf(ctx, T * ctx.climate.varday_T[i]);
    Rdayleaf_daily *= 0.0417;
    return Rdayleaf_daily;
}

// Calculation of above ground biomass (in kg)
//  !!!: if updated, also update CalcIncrementDBH, cf. below
float Tree::CalcAGB(Context &ctx)
{
    // allometric equations from Chave et al. 2014 Global Change Biology to compute above ground biomass (conversion from dbh^2 in cm2 to m2, ie. factor e4)
    // float agb = 0.0673*pow(t_wsg*t_height*ctx.grid.LV*t_dbh*t_dbh*ctx.grid.LH*ctx.grid.LH*10000.0, 0.976);
    float agb = 0.0559 * t_wsg * t_height * ctx.grid.LV * t_dbh * t_dbh * ctx.grid.LH * ctx.grid.LH * 10000.0; // simplified assumption of cylinder, 0.0559 accounts for stem taper cf. Chave et al. 2014)
    return (agb);
}

// Calculation of the increment of dbh from assimilated carbon/biomass (in ctx.params.m)
float Tree::CalcIncrementDBH(Context &ctx, float delta_agb)
{
    float ddbh = fmaxf((delta_agb / (0.559 * t_wsg * 1.0e6 * t_dbh * ctx.grid.LH * t_height * ctx.grid.LV * (3.0 - t_dbh / (t_dbh + t_ah)))), 0.0) * ctx.grid.NH;
    return (ddbh);
}

// Calculation of the maximum amount of carbon stored in a tree
float Tree::CalcCarbonStorageMax(Context &ctx)
{
    // we assume non-structural carbohydrates (NSC) is around 10% of the whole tree's carbon. We assume that only half of NSC (i.e. 5% of whole tree carbon) can be remobilized, since NSC has important metabolic functions or can be stored in tissues that are no longer accessible (heartwood). So NSC for remobilization should never exceed 0.05 * above ground biomass (in units of carbon). Cf. Martínez-Vilalta 2016, Ecological Monographs
    // to calculate the 5% value, we use CalcAGB(), and reconvert kg into gram (1000.0), then convert biomass into carbon by multiplying with 0.5, and take 5%
    // newIM 2021: trees actually store a lot of NSC in roots ! and the 10% estimate account for belowground stock. Add a factor 1.25 to account for belowground biomass, with 0.25 being an estimate of root:shott ration from Ledo et al. 2018 New Phytologist (although this ratio  varies with dbh, species and environment...).
    float carbon_storage_max = 1000.0 * CalcAGB(ctx) * 0.5 * 0.05 * 1.25;
    return (carbon_storage_max);
}

// Calculation of the treefall threshold, if ctx.opt._BASICTREEFALL is activated
// Slightly updated in v.3.1
float Tree::CalcCt(Context &ctx)
{
    float dbhrealmax = t_dbhmax * 1.5;
    float hrealmax = t_mult_height * CalcHeightBaseline(t_ah, t_hmax, dbhrealmax); // realized maximum height
    float vC_intraspecific = ctx.params.vC / 1.5 - 1.0 / (2.3 * t_mult_height) + 1.0 / 2.3;   //! since v.2.5: adjusting ctx.params.vC to intraspecific height variation. If ctx.params.vC was not modified, tall trees would start falling at much larger heights than smaller trees of the same species and with the same dbh, despite a much worse height/dbh ratio. The default assumption is now that the minimum onset of treefalls should be around the same height threshold irrespective of the height multiplier, but stronger assumptions would be justified too (i.e. tall trees falling more easily). The formula is derived as follows: assuming that the onset of treefall can be described by the 99.5 percentile of the sqrt(-log(uniform)) distribution, which is 2.3, we calculate the corresponding Ct_min and impose the condition that it stays equal irrespective of t_mult_height. In this case, we can solve for vC_intraspecific. This is a conservative assumption, likely a stronger dependence on dbh/height ratio would be found, but probably superseded by E-Ping's module anyways
    float Ct = fminf(float(ctx.grid.HEIGHT - 1), hrealmax * fmaxf(1.0 - vC_intraspecific * sqrt(-log(gsl_rng_uniform_pos(ctx.rng.gslrand))), 0.0));
    return (Ct);
}

// Determines leaf life span, either from empirical function or from Kikuzawa model
void Tree::CalcLeafLifespan(Context &ctx)
{
    if (ctx.opt._LL_parameterization == 0)
    {                                                                   // prescribed relationship for LL
        t_leaflifespan = 12.755 * exp(0.007 * t_LMA - 0.565 * t_Nmass); // expression developed by Sylvain Schmitt, avoids problem of very low leaflifespans at low LMA
        t_leaflifespan = fmaxf(t_leaflifespan, 3.0);
        // t_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA))));    //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).
        // t_leaflifespan=1.5+pow(10,(7.18+3.03*log10(t_LMA*0.0001)));   //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
        // t_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(t_LMA)));        //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).
    }
    else
    { // relationship based on optimal theory for LL
        t_leaflifespan = predLeafLifespanKikuzawa(ctx);
    }

    t_leaflifespan *= 0.08333333 * ctx.time.iterperyear; // Converts leaflifespan from month unit to ctx.time.timestep units, this is needed for UpdateLeafDynamics and nppneg (0.08333333=1/12)

    // float time_young=fminf(t_leaflifespan/3.0,1.0);
    float time_young = fminf(t_leaflifespan / 3.0, 1.0 * 0.08333333 * ctx.time.iterperyear); // modified IM jan23
    float time_mature = t_leaflifespan / 3.0;
    float time_old = t_leaflifespan - time_mature - time_young;

    t_lambda_young = 1.0 / time_young;
    t_lambda_mature = 1.0 / time_mature;
    t_lambda_old = 1.0 / time_old;
}

void Tree::InitialiseLeafPools(Context &ctx)
{
    t_youngLA = t_LA / (t_lambda_young * t_leaflifespan);
    t_matureLA = t_LA / (t_lambda_mature * t_leaflifespan);
    t_oldLA = t_LA / (t_lambda_old * t_leaflifespan);
    t_litter = 0.0;
}

// Determine sapwood area, limited by increase in dbh (ddbh) (in m2)
//! - Options: fixed percentage (option ctx.opt._sapwood == 0) or based on tree's leaf area (option ctx.opt._sapwood > 0). In the case ctx.opt._sapwood>0, a tree cannot retroactively convert heartwood into sapwood, just because it could allocate more leaves
//! - TODO: should there be another limit on leaf area? (i.e. if sapwood cannot grow more than ddbh increment, then leaf area should probably not be allowed to grow beyond what is reasonable through the current sapwood area)
void Tree::UpdateSapwoodArea(Context &ctx, float ddbh)
{
    if (ctx.opt._sapwood > 0)
    {
        float sapwood_area_new = PI * 0.5 * ddbh * (t_dbh - 0.5 * ddbh); // correction from previous equation in v.3.0, derived from (0.5 * t_dbh) * (0.5 * t_dbh) * PI - (0.5 * (t_dbh - ddbh)) * (0.5 * (t_dbh - ddbh)) * PI, this presupposes that t_dbh has already been updated (i.e. t_dbh += ddbh). If dbh has grown entirely from zero, then the equation reduces to PI * 0.5 * dbh * 0.5 * dbh, i.e. the whole stem area
        t_sapwood_area += sapwood_area_new;
        t_sapwood_area = fminf(t_sapwood_area, 0.0001 * 2.0 * t_LA / (0.066 + 0.017 * t_height - 0.18 + 1.6 * t_wsg)); // upper bound on sapwood area, either through previously existing sapwood area (sapwood cannot grow quicker than the rest of the tree) or through the amount of sapwood needed from from Fyllas et al. 2014, based on inversion of pipe model, multiplication with 0.0001 to convert cm2 to m2
        float sapwood_minimum;
        if (t_dbh < 0.01)
            sapwood_minimum = t_dbh * t_dbh * 0.25 * PI;
        else
            sapwood_minimum = 0.005 * (t_dbh - 0.005) * PI;
        t_sapwood_area = fmaxf(t_sapwood_area, sapwood_minimum);
    }
    else
    {
        float sapthick;
        if (t_dbh < 0.08)
            sapthick = 0.5 * t_dbh;
        else
            sapthick = 0.04;
        t_sapwood_area = PI * sapthick * (t_dbh - sapthick);
    }
}

// Updates t_height, based on t_dbh
void Tree::UpdateHeight(Context &ctx)
{
    float height_baseline = CalcHeightBaseline(t_ah, t_hmax, t_dbh);
    t_height = fminf(t_mult_height * height_baseline, ctx.grid.HEIGHT - 1);
}

// Updates t_CR, based on t_dbh
void Tree::UpdateCR(Context &ctx)
{
    t_CR = CalcCRBaseline(ctx, t_dbh) * t_mult_CR;
    t_CR = fmaxf(ctx.params.CR_min, t_CR);
}

// Updates t_CD based on t_height
void Tree::UpdateCD(Context &ctx)
{
    // Since v.2.5, simplification of the computation of the crown depth, in accordance with the Canopy Constructor algorithm
    t_CD = CalcCDBaseline(ctx, t_height) * t_mult_CD;
    t_CD = fminf(t_CD, 0.5 * t_height);
}

// Determines the maximum LAI that the tree should reach, given the Farquhar model, at a theoretical average day, up to which point leaves can be allocated and until which the costs (self-shading) are lower than the benefits (additional assimilation)
//! - LAImax (maximum LAI) is assumed to lie between 0 and 10, and is such that adding more leaves results in net loss of carbon
//! - Range of LAImax is narrowed down by bisection (sequentially halving the possible range)
//! - Uses actually absorbed PPFD instead of incident PPFD, since leaves are not perfectly illuminated, but often in the lower canopy layers
void Tree::CalcLAImax(Context &ctx)
{
    float LAI_lowerbound = 0.0;
    float LAI_upperbound = 10.0;
    float LAImax_temp = 0.5 * (LAI_lowerbound + LAI_upperbound);

    // tree LAImax to not be greater 10. 10 is also the number of iterations to get the best value. Since we narrow down in half steps the precision should be ca. 0.5^10 ~ 0.001, which should be high enough for LAImax
    for (int i = 0; i < 10; i++)
    {
        float absorb_prev = LAImax_temp;
        float absorb_delta = 0.5; // calculate everything, assuming a medium leaf density in the lower canopy, TODO: if the crown has a non-cylindric shape, this should probably be taken as a third of LAImax, or if LAI_gradient is activated, as 25% of LAImax
        int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);

#ifdef FULL_CLIMATE // if FULL_CLIMATE is defined, then the terms just below should be only the relative attenuation due to the vertical structure of the canopy. Currently, CalcLAImax is defined depending on the climate conditions on the day of birth (see ctx.time.iter, in dailyFluxesLeaf), and not to the climate yearly mean as before. This needs to be rethought and changed.

        // get PPFD, VPD, and temperature at each discretisation step
        float PPFD_LAI = ctx.lookup.LookUp_flux_absorption[intabsorb];
        float VPD_LAI = ctx.lookup.LookUp_VPD[intabsorb];
        float Tmp_LAI = ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
        int intincident = CalcIntabsorb(absorb_prev);
        float Wind_LAI = exp(-0.5 * absorb_prev); // to be thoroughly checked and computed using a look up table as well.
        float ExtinctLW_LAI = ctx.lookup.LookUp_ExtinctLW[intincident];
        float PPFD_LAI_inc = ctx.lookup.LookUp_flux[intincident];
#endif // WATER

#else // FULL_CLIMATE

        // get PPFD, VPD, and temperature at each discretisation step
        float PPFD_LAI = ctx.climate.WDailyMean_year * ctx.lookup.LookUp_flux_absorption[intabsorb];
        float VPD_LAI = ctx.climate.VPDDailyMean_year * ctx.lookup.LookUp_VPD[intabsorb];
        float Tmp_LAI = ctx.climate.tDailyMean_year - ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
        int intincident = CalcIntabsorb(absorb_prev);
        float Wind_LAI = ctx.climate.windDailyMean_year * exp(-0.5 * absorb_prev); // to be thoroughly checked and computed using a look up table as well.
        float ExtinctLW_LAI = ctx.lookup.LookUp_ExtinctLW[intincident];
        float PPFD_LAI_inc = ctx.climate.WDailyMean_year * ctx.lookup.LookUp_flux[intincident];
#endif // WATER

#endif // FULL_CLIMATE

        // calculate the GPP
#ifdef WATER
        float wf = 0.0;
        float GPP_LAI = Tree::dailyFluxesLeaf(ctx, PPFD_LAI, VPD_LAI, Tmp_LAI, Wind_LAI, ExtinctLW_LAI, PPFD_LAI_inc, wf).carbon_flux;
#else
        float GPP_LAI = Tree::dailyGPPleaf(ctx, PPFD_LAI, VPD_LAI, Tmp_LAI);
        float Rday_LAI = Tree::dailyRdayleaf(ctx, Tmp_LAI);
#endif
        float effLA = 0.66 * ctx.time.nbhours_covered * 15.7788 * ctx.time.timestep;                // convert  from micromoles C/ctx.params.m^2/s into gC per ctx.params.m^2 of leaf per ctx.time.timestep by "ctx.time.nbhours_covered*15.7788*ctx.time.timestep" where 15.7788 = 3600*365.25*12/1000000 (seconds, days, and mass of carbon) and ctx.time.nbhours_covered is the amount of time that is covered by the daily variation file. We also assume that one third of the leaves are mature and that the rest of the leaves have half the assimilation rates, so we derive a factor 0.66
        float effLA_night = 0.83 * (24.0 - ctx.time.nbhours_covered) * 15.7788 * ctx.time.timestep; // same as during the day, but inverse of hours covered (we assume that non-covered hours are night values), assuming that respiration rate of yound and old leaves are 75% that of mature leaves.

        // float effLA = 0.5 * ctx.time.nbhours_covered * 15.7788 * ctx.time.timestep; //convert  from micromoles C/ctx.params.m^2/s into gC per ctx.params.m^2 of leaf per ctx.time.timestep by "ctx.time.nbhours_covered*15.7788*ctx.time.timestep" where 15.7788 = 3600*365.25*12/1000000 (seconds, days, and mass of carbon) and ctx.time.nbhours_covered is the amount of time that is covered by the daily variation file. here we assume that the leaves that will determine the LAImax are the one at the bottom of the crown and taht these are all old leaves, with hald the assimilation of mature leaves. Note that this is not the case for phenological strategies that exchange all their leaves. So when moving to a variable t_LAImax, we should account for the leaf area of the last ctx.diag.layer and fill it with old leaves at maximum.
        // float effLA_night = 0.75 * (24.0 - ctx.time.nbhours_covered) * 15.7788 * ctx.time.timestep;  //same as during the day, but inverse of hours covered (we assume that non-covered hours are night values), assuming that respiration rate of yound and old leaves are 75% that of mature leaves.

        GPP_LAI *= effLA;

        // get the night respiration
#ifdef FULL_CLIMATE
        int convTnight = int(ctx.lookup.iTaccuracy * ctx.climate.tnight);
#else
        int convTnight = int(ctx.lookup.iTaccuracy * ctx.climate.Tnight_year);
#endif
        float Rnight_LAI = t_Rdark * effLA_night * ctx.lookup.LookUp_Rleaf[convTnight];

        // add up the two components of leaf respiration, and multiply the result by 1.5 (fine root respiration), since in TROLL, this cannot be separated from leaf respiration
#ifdef WATER
        float Rleaf_LAI = Rtot_by_Rabove * Rtotleaf_by_Rdark * Rnight_LAI + (Rtotleaf_by_Rdark - 1) * Rnight_LAI;
#else
        Rday_LAI *= effLA * 0.4;                                    // inhibition of respiration by ca. 40%, cf. Atkin et al. 2000
        float Rleaf_LAI = Rtot_by_Rabove * (Rday_LAI + Rnight_LAI); // Rleaf=Rday+Rnight is multiplied by 1.5 to also account for fine root respiration (cf as in Fyllas et al 2014 and Malhi 2012).
#endif

        // calculate effective npp
        float NPP_leaf = 0.7 * (GPP_LAI - Rleaf_LAI); // growth respiration, v. 2.4.1: 0.75 replaced by 0.7. According to Cannell and Thornley 2000, higher values should be used only if other sources of respiration (phloem loading etc.) are explicitly accounted for. For leaf construction, typically even a lower factor is a assumed (cf. Villar & Merino 2001, avg construction costs of 1.52 gC/gC)

        // update boundaries
        if (NPP_leaf > 0.0)
            LAI_lowerbound = LAImax_temp;
        else
            LAI_upperbound = LAImax_temp;

        // now update LAImax_temp
        LAImax_temp = 0.5 * (LAI_lowerbound + LAI_upperbound);
    }

    t_LAImax = LAImax_temp;
}

// Determines the maximum LAI that the tree should reach, similar to CalcLAImax, but it directly calculates the maximum leafarea for the current light environment the tree experiences rather than for theoretical day
//! - the idea is to compute the LAI that the tree experiences at the crown top and compare it to its maximum LAI, then the tree adjusts its LAI or leaf area accordingly
//! - importantly, since PPFD is an exponential function of leaf density (nonlinear), averaging LAI across the whole crown would not correspond to the LAI experienced by the tree. A heterogeneous environment with a mixture of light specks and dense spots above the crown will contribute considerably more PPFD (and have a lower effective LAI) than a homogeneously filled crown. We calculate the effective LAI instead of the real LAI above the crown top by first calculating average PPFD and then converting back to LAI
void Tree::CalcLAmax(Context &ctx, float &LAIexperienced_eff, float &LAmax)
{
    float crown_area = PI * t_CR * t_CR;
    float crown_area_nogaps = GetCrownAreaFilled(crown_area);

    if (ctx.opt._LA_regulation == 1)
    {
        // this is the case where we limit leaf allocation to the tree's maximum LAI, but have no dynamic reactions to the environment
        LAIexperienced_eff = 0.0;
    }
    else
    {
        int site_crowncenter = t_site + t_CrownDisplacement;
        int row_crowncenter = site_crowncenter / ctx.grid.cols;
        int col_crowncenter = site_crowncenter % ctx.grid.cols;

#ifdef CROWN_UMBRELLA
        // get the PPFD that the tree experiences and the area that is looped
        float ppfd_CA[2] = {0.0, 0.0};
        float noinput = 0.0;

        int shell_fromtop = 0;
        float fraction_filled_target = t_fraction_filled;

        LoopLayerUpdateCrownStatistic_template(ctx, row_crowncenter, col_crowncenter, t_height, t_CR, t_CD, fraction_filled_target, shell_fromtop, [&ctx](float CR, float e, float p){ return GetRadiusSlope(ctx, CR, e, p); }, noinput, ppfd_CA, KeepFloatAsIs, [&ctx](int h, int s, float ni, float (&ca)[2]){ GetPPFDabove(ctx, h, s, ni, ca); });

        float ppfd_experienced = ppfd_CA[0];
        float crown_area_looped = ppfd_CA[1];
#else
        int height_abovetop = int(t_height) + 1;
        int crown_intarea = GetCrownIntarea(t_CR);

        float fraction_filled_actual = 0.0;
        float fraction_filled_target = t_fraction_filled;

        float ppfd_experienced = 0.0;
        int crown_area_looped = 0;

        // loop over LookUp table, until crown_intarea is reached
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
                    // first get voxel field densities
                    float absorb_prev = ctx.field.LAI3D[height_abovetop][site + ctx.grid.SBORD];
                    int intabsorb = CalcIntabsorb(absorb_prev);

                    // obtain PPFD for the voxel, and also record the circled area
                    ppfd_experienced += ctx.climate.WDailyMean * ctx.lookup.LookUp_flux[intabsorb];
                    crown_area_looped++;
                }
            }
        }
#endif
        if (crown_area_looped > 0)
        {
            float icrown_area_looped = 1.0 / crown_area_looped;
            ppfd_experienced *= icrown_area_looped;
        }
        else
        {
            // this should not happen, potentially introduce a break or assert() here
            ppfd_experienced = ctx.climate.WDailyMean;
        }

        // now calculate the effective LAI above the tree
        // importantly, we use the yearly average here, as the maximum tree LAI is also defined with respect to yearly averages
        LAIexperienced_eff = -log(ppfd_experienced / ctx.climate.WDailyMean_year) / ctx.params.kpar;
        if (LAIexperienced_eff < 0.0001)
            LAIexperienced_eff = 0.0;
    }

    // now take the difference between the maximum value and the actual LAI experienced. Will be positive, if there is still room for leaves, and negative, if it's already too dense
    // then calculate the relative size of the difference with respect to the tree's own LAI and calculate by how much it can still increase or has to reduce its leafarea
    float LAIdiff = t_LAImax - LAIexperienced_eff; // this is the amount of LAI the tree has to loose or can add relative to its currently experienced LAI. It translates directly into the relative leafarea the tree can add or should lose, since all the calculations before are linear functions of the density
    float LAdiff = LAIdiff * crown_area_nogaps;

    // in the worst case (if the difference in LAI cannot be compensated for by the tree anymore, i.e. it would have to have negative leafarea), we will get zero LA_max, and the tree will allocate everything to storage
    float LAmaxphysiology = t_LAImax * crown_area_nogaps;

    LAmax = fminf(LAmaxphysiology, LAdiff);
    LAmax = fmaxf(0.0, LAmax);

#ifdef TRACK_INDIVIDUALS
    if (t_month_born >= 0)
    {
        t_LAIabove_effavgyear += LAIexperienced_eff * ctx.time.timestep;
        t_LAIeffcum += LAIexperienced_eff;
        t_LAIeffsquared_cum += LAIexperienced_eff * LAIexperienced_eff;
    }
#endif
}

// Kikuzawa model for leaf lifespan
//! - based on Kikuzawa, K. (1991). A cost-benefit analysis of leaf habit and leaf longevity of trees and their geographical pattern. The American Naturalist, 138(5), 1250-1263.
//! - the model assumes that the daily photosynthesis of leaves declines (linearly) with leaf age, then calculates the optimal leaf lifespan, i.e. the one that maximizes lifetime carbon gain for a leaf
//! - this model is here implemented by calculating the LAI that the tree could have under optimal conditions (i.e. no leaf area above, and a LAI equivalent to LAImax), then calculating the effective daily assimilation (GPP) per unit leaf area (i.e. excluding stem respiration), and solving the model equation for leaf lifespan
//! - there are a few caveats/things to consider in further developing/exploring this model a) A crucial parameter is the leaf age when carbon assimilation would reach 0, often denoted as b. This parameter has been shown to correlate with LMA in Xu et al. 2017, but, like LL-SLA relations, the relation comes with a lot of unexplained variation (on log scales, so multiplicative), and when I (FF) looked at their data, I obtained a better fit with Vcmax_25_mass, which we calculate anyways in TROLL, so it seemed to be internally consistent to use this relation. However, variation is still enormous (exp(0.6)), and this is clearly far from mechanistic. So while the Kikuzawa model partially solves the problem of deriving LL from LMA and the large uncertainties as well as potentially misleading inferences (i.e. very short leaf lifespans), part of the problem with large variation and weak correlations still remains and is now just shifted onto the derivation of b! Any more research on b and how to derive it for trees would thus be important  b) A second, crucial question is, what exactly the costs are that should be considered in the cost-benefit analysis. While leaf construction cost is uncontroversial (here taken to be LMA, i.e. mass per unit leaf area * 1.5 gC/gC to account for additional respiration in construction), there is a host of other costs that could be considered. We here also consider costs for leaf maintenance and for fine roots transpiration, but do not account for other cost such as building of support structures (Kikuzawa and Ackerly 2002, https://esj-journals.onlinelibrary.wiley.com/doi/abs/10.1046/j.1442-1984.1999.00005.x). Conversely, some papers, discount any respiration (e.g. Wang et al. 2021, https://www.biorxiv.org/content/10.1101/2021.02.07.430028v1.full). This is fundamentally about how modular/interconnected a tree is, and whether leaves could be modelled, at least as a first order approximation, as isolated. c) There is a choice to be made about what kind of environment should be supposed for leaf life span calculations. Are tree leaves optimized for an optimal light environment or the light environment when they are produced?  d) Finally, a minor point: I simply transferred the linear model of decline in leaf photosynthesis from Kikuzawa, etc. into TROLL. But this should probably be checked again, since we have an alternative way of simulating decline in leaf photosynthesis through young, mature and old leaves
float Tree::predLeafLifespanKikuzawa(Context &ctx)
{

    float absorb_prev = 0.0;
    float absorb_delta;
    if (ctx.opt._LA_regulation > 0)
        absorb_delta = t_LAImax;
    else
        absorb_delta = 3.0;

    int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);

#ifdef FULL_CLIMATE

    float PPFD = ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD = ctx.lookup.LookUp_VPD[intabsorb];
    float T = ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float W = exp(-0.5 * absorb_prev); // to be thoroughly re-thought, and computed using a look-up table as well -- IM June 2021
    int intincident = CalcIntabsorb(absorb_prev);
    float PPFDinc = ctx.lookup.LookUp_flux_absorption[intincident];
    float ExtinctLW = ctx.lookup.LookUp_ExtinctLW[intincident];
    float wf = 0.0;
    float GPP = Tree::dailyFluxesLeaf(ctx, PPFD, VPD, T, W, ExtinctLW, PPFDinc, wf).carbon_flux;
#else
    float GPP = Tree::dailyGPPleaf(ctx, PPFD, VPD, T);
    float Rday = Tree::dailyRdayleaf(ctx, T) * 0.4; // inhibition of respiration by ca. 40%, cf. Atkin et al. 2000
#endif // WATER

    int convTnight = int(ctx.lookup.iTaccuracy * ctx.climate.tnight);

#else // FULL_CLIMATE

    float PPFD = ctx.climate.WDailyMean_year * ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD = ctx.climate.VPDDailyMean_year * ctx.lookup.LookUp_VPD[intabsorb];
    float T = ctx.climate.tDailyMean_year - ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float W = ctx.climate.windDailyMean_year * exp(-0.5 * absorb_prev); // to be thoroughly re-thought, and computed using a look-up table as well -- IM June 2021
    int intincident = CalcIntabsorb(absorb_prev);
    float PPFDinc = ctx.climate.WDailyMean_year * ctx.lookup.LookUp_flux_absorption[intincident];
    float ExtinctLW = ctx.lookup.LookUp_ExtinctLW[intincident];
    float wf_fc = 0.0;
    float GPP = Tree::dailyFluxesLeaf(ctx, PPFD, VPD, T, W, ExtinctLW, PPFDinc, wf_fc).carbon_flux;
#else
    float GPP = Tree::dailyGPPleaf(ctx, PPFD, VPD, T);
    float Rday = Tree::dailyRdayleaf(ctx, T) * 0.4; // inhibition of respiration by ca. 40%, cf. Atkin et al. 2000
#endif // WATER

    int convTnight = int(ctx.lookup.iTaccuracy * ctx.climate.Tnight_year);

#endif // FULL_CLIMATE
    float Rnight = t_Rdark * ctx.lookup.LookUp_Rleaf[convTnight];

#ifdef WATER
    float GPP_effective = (GPP - Rtot_by_Rabove * Rtotleaf_by_Rdark * Rnight + (Rtotleaf_by_Rdark - 1) * Rnight);
#else
    float GPP_effective = (GPP - Rtot_by_Rabove * Rday - Rtot_by_Rabove * Rnight); // in gC ctx.params.m^-2 day-1  we assume that Rday/Rnight decreases linearly as well, so we will use a GPP_eff = GPP - Rday;
#endif

    GPP_effective *= ctx.time.nbhours_covered * 3600.0 * 12.0 / 1000000.0; // we convert micromoles C/ctx.params.m^2/s into gC/ctx.params.m^2/day (factor 12.0 for conversion into gC, 3600 for second-hour conversion, and 10^6 to convert micromoles to moles; finally, we divide by the total leaf area that intercepted the light, since

    int convT = int(ctx.lookup.iTaccuracy * 25.0); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    float Vcmax_25_mass = t_Vcmax * ctx.lookup.LookUp_VcmaxT[convT] / t_LMA;

    // Core model
    float b = exp(5.467025 - 1.138354 * log(Vcmax_25_mass) + gsl_ran_gaussian(ctx.rng.gslrand, 0.6112195)); // parameter as in Xu et al. 2017, fit to their data
    float LL = 1.0 + 0.0333333 * fminf(b, sqrt(3.0 * 0.5 * t_LMA * b / GPP_effective));             // cf. Kikuzawa 1991 for formula, and Xu et al. 2017. As opposed to Xu et al. 2017 we assume decline in respiration rates. Also, LL can never be bigger than b
    LL = fmaxf(LL, 3.0);

    return (LL);
}

// #############################################
//  Tree growth
// #############################################
void Tree::Growth(Context &ctx)
{
    // update age
    t_age += ctx.time.timestep; // new v.2.2: increments are not 1 yr, but the duration of the ctx.time.timestep (usually 1 or <1, i.e. 1/12 if monthly, 1/365 if daily

    // Set growth carbon to zero
    t_carbon_biometry = 0.0;
#ifdef WATER
    t_transpiration = 0.0;
#endif

    // calculate GPP and respiration
    CalcRespGPP(ctx);

    //! calculate NPP
    CalcNPP(ctx);

#ifdef CHECK_CARBON
    carbon_net_total += t_NPP;
    carbon_assimilated_total += t_GPP;
#endif

    if (ctx.opt._LA_regulation == 0)
    {
        // Classic case of no regulation of leaf area, all NPP > 0.0 is converted into leaves
        if (t_NPP < 0.0)
        {
            t_NPPneg++;
#ifdef TRACK_INDIVIDUALS
            if (t_month_born >= 0)
            {
                t_time_carbonstarvation_year++;
                t_time_carbonstarvation++;
            }
#endif
            t_NPP = 0.0;
        }
        else
        {
            t_NPPneg = 0;
            // NPP allocation to wood and tree size increment
            UpdateTreeBiometry(ctx);
        }

        // NPP allocation to leaves
        UpdateLeafDynamics(ctx);
    }
    else
    {
        // leaf area is dynamically adjusted
        if (t_NPP < 0.0)
        {
            // as long as the tree has enough storage to compensate for negative NPP, we will not count it as negative NPP
            // the storage will, however, not be used for growth, only for upkeep
            float carbon_storage_leftover = t_NPP + t_carbon_storage;
            if (carbon_storage_leftover > 0.0)
            {
                t_NPP = 0.0;
                t_carbon_storage = carbon_storage_leftover;
            }
            else
            {
                t_carbon_storage = 0.0;
            }
        }

        // if NPP is still negative, even after compensating for it through stored carbon, then carbon starvation sets in.
        if (t_NPP < 0.0)
        {
            // t_NPPneg++; //newIM 2021: if NPP still negative, this means that the tree has completly depleted its minimal carbon stock pool, ie. is not viable anymore (Martinez-Vilalta et al. 2016) ==> death at the next ctx.time.timestep (the previous version made a mix between the use of NPPneg (that was introduced in absence of a variable for carbon storage and a dynamic monitoring and allocation to this carbon storage pool: NPPneg is not needed anymore when ctx.opt._LA_regulation >0)
            t_NPP = 0.0;
            // v.2.3.0 -- Line of code below was odd. If NPP <0.0, then to ensure C balance it should be simply reset to NPP=0 at this stage
            // t_NPP=t_GPP - 1.5*(t_Rday+t_Rnight+t_Rstem); REMOVED AS OF v.2.3.0.a4
            // NPP allocation to leaves
            UpdateLeafDynamics(ctx);
        }
        else
        {
            // t_NPPneg=0; //newIM 2021: no use of NPPneg anymore when ctx.opt._LA_regulation >0
            // in v.2.4.1: allocation to leaves is done before updating biometry
            // idea is: at least leaves need to be sustained, otherwise the tree cannot come back from carbon stress
            // NPP allocation to leaf, storage, wood and tree size increment
            UpdateLeafDynamics(ctx);
            UpdateTreeBiometry(ctx);
        }
    }

    UpdateVolumeDensity();

#ifdef WATER
    if (ctx.time.iter == (ctx.time.nbiter - 90))
        OutputTreeStandard(ctx, ctx.out.output[28]);
    if (ctx.time.iter == (ctx.time.nbiter - 45))
        OutputTreeStandard(ctx, ctx.out.output[29]);
    if (ctx.time.iter == (ctx.time.nbiter - 1))
        OutputTreeStandard(ctx, ctx.out.output[30]);

    if (t_site == 4)
        OutputTreeStandard(ctx, ctx.out.output[12]);
    if (t_site == 10380)
        OutputTreeStandard(ctx, ctx.out.output[13]);
    if (t_site == 100950)
        OutputTreeStandard(ctx, ctx.out.output[14]);
    if (t_site == 12090)
        OutputTreeStandard(ctx, ctx.out.output[15]);
    if (t_site == 120090)
        OutputTreeStandard(ctx, ctx.out.output[16]);
    if (t_site == 150667)
        OutputTreeStandard(ctx, ctx.out.output[17]);
#endif
}

// ###################################################
//  GPP and respiration calculation, called by Tree::Growth
// #####################################################
//! - new v.2.4: t_GPP and t_Rday are updated in separate function (thus adaptable for different modules)
//! - but this is not the case when WATER is activated: indeed Rday is now computed within Photosyn and FluxesLeaf functions and thus included in GPP computation.
void Tree::CalcRespGPP(Context &ctx)
{

#ifdef WATER
    float tree_transpiration_1016 = 0.0;
#endif

    if (t_LA > 0.0)
    {
#ifdef WATER
#else
        if (ctx.opt._GPPcrown)
        {
            // v.2.3.1 -- fast GPP calculation option.
            float PPFD = 0.0, VPD = 0.0, Tmp = 0.0, leafarea_layer;
            Fluxh(ctx, int(t_height) + 1, PPFD, VPD, Tmp, leafarea_layer);
            t_GPP = Tree::dailyGPPcrown(ctx, PPFD, VPD, Tmp, t_LAI);
            t_Rday = Tree::dailyRdayleaf(ctx, Tmp);
        }
        else
        {
#endif
        // Photosynthesis and respiration, set to zero
        t_GPP = 0.0;
#ifdef WATER // the computation of leaf respiration during the day is accounted for in the Photosyn function
#else
            t_Rday = 0.0;
#endif

        int crown_above_base = int(t_height - t_CD) + 1; // for flux above crown base
        int crown_above_top = int(t_height) + 1;         // for flux above crown top

#ifdef CROWN_UMBRELLA
        int h_stop = max(crown_above_base, crown_above_top - 3);
#else
            int h_stop = crown_above_base;
#endif
        float leafarea_cumulated = 0.0;

#ifdef WATER
        float leafarea_toprounding = 0.0; // newIM
        float MeanDCELLHeight = ctx.soil.Canopy_height_DCELL[ctx.grid.site_DCELL[t_site]];
#endif

        for (int h = crown_above_top; h >= h_stop; h--)
        {
            float PPFD = 0.0, VPD = 0.0, Tmp = 0.0, leafarea_layer = 0;
#ifdef WATER
            float PPFD_incident = 0.0, ExtinctLW = 0.0;
            // float W=WS*exp(-0.5*ctx.soil.LAI_DCELL[h][ctx.grid.site_DCELL[t_site]]); // computation of wind speed, which is here assumed to declined exponentially with the cumulative LAI within a given neighborhood (here taken as the tree's dcell) (see equation B11 in Medvigy et al. 2009; see Leuning et al. 1995 PCE equ. E2). We could alternatively used the aerodynamic momentum transfer model (Monteith and Unsworth 2008) using Lorey'sheight (see E-Ping's second manuscript) -- to be discussed. // to be computed with a lookup table maybe, and within Fluxh probably.
            int convHratio = int(ctx.lookup.iHaccuracy * h / MeanDCELLHeight);
            float W = ctx.soil.TopWindSpeed_DCELL[ctx.grid.site_DCELL[t_site]] * ctx.lookup.LookUp_Wind[convHratio];
            if (W <= 0)
            {
                cout << " Wind=" << W << " h=" << h << " MeanDCELLHeight=" << MeanDCELLHeight << " convHratio=" << convHratio << " ctx.lookup.LookUp_Wind[convHratio]=" << ctx.lookup.LookUp_Wind[convHratio] << endl;
            }

            // cout << " Wind=" << W << " h=" << h << " MeanDCELLHeight=" << MeanDCELLHeight <<" convHratio=" << convHratio << " ctx.lookup.LookUp_Wind[convHratio]=" << ctx.lookup.LookUp_Wind[convHratio] << " ctx.soil.TopWindSpeed_DCELL[ctx.grid.site_DCELL[t_site]]=" << ctx.soil.TopWindSpeed_DCELL[ctx.grid.site_DCELL[t_site]] << endl;

            Fluxh(ctx, h, PPFD, VPD, Tmp, leafarea_layer, PPFD_incident, ExtinctLW);
#else

                Fluxh(ctx, h, PPFD, VPD, Tmp, leafarea_layer);
#endif

#ifdef WATER
            if (leafarea_toprounding > 0)
            {
                // cout << "Rounding when h=" <<h << " with crown_above_top=" << crown_above_top << " and leafarea_toprounding=" << leafarea_toprounding << endl;
                leafarea_layer += leafarea_toprounding;
                leafarea_toprounding = 0.0;
            }

            if (leafarea_layer < 0.0)
            {
                cout << "Warning leafarea <=0 in CalcRespGPP !! PPFD=" << PPFD << "; leafarea_layer=" << leafarea_layer << " when h=" << h << " with crown_above_top=" << crown_above_top << " and leafarea_toprounding=" << leafarea_toprounding << endl;
            }

            if (leafarea_layer > 0.0 && PPFD > 0)
            {
                if (PPFD < 0)
                {
                    cout << "Warning PPFD <0 in CalcRespGPP !! PPFD=" << PPFD << endl;
                }
                float water_flux_1016 = 0.0;
                leafFluxes dailyF = Tree::dailyFluxesLeaf(ctx, PPFD, VPD, Tmp, W, ExtinctLW, PPFD_incident, water_flux_1016);
                t_GPP += leafarea_layer * dailyF.carbon_flux;
                t_transpiration += leafarea_layer * dailyF.water_flux;
                tree_transpiration_1016 += leafarea_layer * water_flux_1016;

                if (isnan(t_GPP) || isnan(t_transpiration) || t_transpiration < 0.0)
                    cout << "Problem at site: " << t_site << ", t_GPP= " << t_GPP << ", transpiration= " << t_transpiration << ", leafarea_layer= " << leafarea_layer << ", PPFD=" << PPFD << ", VPD= " << VPD << ", T= " << Tmp << ", t_WSF_A=" << t_WSF_A << ", t_WSF= " << t_WSF << endl;

#else
                t_GPP += leafarea_layer * Tree::dailyGPPleaf(ctx, PPFD, VPD, Tmp);
                // if(isnan(t_GPP) || t_GPP < 0.0 || PPFD == 0.0) cout << "Problem at site: " <<  t_site << " leafarea_layer: " << leafarea_layer << " GPP: " << t_GPP << " PPFD: " << PPFD << " VPD: " << VPD << " T: " << Tmp << endl << endl;
#endif

#ifdef WATER
#else
                t_Rday += leafarea_layer * Tree::dailyRdayleaf(ctx, Tmp);
#endif
                leafarea_cumulated += leafarea_layer;
                // if(t_dbh > 0.2) cout << crown_above_top - h << "PPFD: " << PPFD << " VPD: " << VPD << " Tmp: " << Tmp << " leafarea_layer: " << endl;

#ifdef WATER
            }
            else
            {
                if (leafarea_layer > 0.0)
                    leafarea_toprounding = leafarea_layer;
            }
#endif
        }

        // Averaging across layers. To check consistency, leafarea_layer can be added up and the sum compared to t_LA
        float ileafarea_cumulated;
        if (leafarea_cumulated > 0.0)
            ileafarea_cumulated = 1.0 / leafarea_cumulated;
        else
            ileafarea_cumulated = 0.0;

        t_GPP *= ileafarea_cumulated;
#ifdef WATER
        t_transpiration *= ileafarea_cumulated;
        tree_transpiration_1016 *= ileafarea_cumulated;
#else
            t_Rday *= ileafarea_cumulated;
        }
#endif
    }

    float effLA = 0.5 * (t_LA + t_matureLA) * ctx.time.nbhours_covered * 15.7788 * ctx.time.timestep;                        // we convert  from micromoles C/ctx.params.m^2/s into gC per ctx.params.m^2 of leaf per ctx.time.timestep by "ctx.time.nbhours_covered*15.7788*ctx.time.timestep" where 15.7788 = 3600*365.25*12/1000000 (seconds, days, and mass of carbon) and ctx.time.nbhours_covered is the amount of time that is covered by the daily variation file.
    float effLA_night = (0.75 * t_LA + 0.25 * t_matureLA) * (24.0 - ctx.time.nbhours_covered) * 15.7788 * ctx.time.timestep; // same as during the day, but inverse of hours covered (we assume that non-covered hours are night values), and assuming that respiration of young and old leaves are 75% the one of mature leaves:  Kitajima et al. 2002 found no or small reductions in respiration with leaf age, Reich et al. 2009 did (Rdark/Amax constant), but often respiration declines are less steep than photosynthetic capacity, cf. Villar et al. 1995.

    t_GPP *= effLA;

#ifdef WATER
    float effLAT = (0.75 * t_LA + 0.25 * t_matureLA) * ctx.time.nbhours_covered * 23.652 * 0.000001 * ctx.time.timestep; // assuming that the decline in stomatal conductance of young and old leaves is less than for photosynthetic capacities (cf. ALbert et al. 2018)(75% the one of mature leaves), ie the water use efficiencies is lower of young and old leaves. Note that idally the Farquhar model should be applied to each leaf cohorts separately... 23.652=18.10^-6 * 3600*365, where 18.10^-6 is to convert micromol H20 into g, 3600*365 to convert s into year (as ctx.time.timestep is given in year); and 0.000001 is to convert g in m3.
    t_transpiration *= effLAT;                                                                         // this is the amount of water transpired by the tree during the ctx.time.timestep in m3.  CHECK: from our WILT data of sapflow of several canopy mature trees, and assuming a constant sapwood thickness of 4cm, this can be 5-110 liter/day, ie. 0.005-0.11 m3/day. See also values in Granier et al. 1996 (0.250-0.3 m3/day for a big dominant canopy tree), or Andrade et al. 1998 Oecologia (46-379 kg/day, ie. 0.046-0.379 m3/day). Note, in this version, we never account for nighttime transpiration... to be discussed.
    tree_transpiration_1016 *= (0.75 * t_LA + 0.25 * t_matureLA) * 0.0324 * 0.000001;                  // where 0.0324= 18.10^-6*3600*0.5
    ctx.soil.transpiration_1016 += tree_transpiration_1016;
#endif

    int convT = int(ctx.lookup.iTaccuracy * ctx.climate.tDailyMean);  // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    int convTnight = int(ctx.lookup.iTaccuracy * ctx.climate.tnight); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---

    t_Rstem = t_sapwood_area * (t_height - t_CD) * 0.5 * (ctx.lookup.LookUp_Rstem[convT] + ctx.lookup.LookUp_Rstem[convTnight]);

    // the following two lines should maybe be revised, because the decline in respiration might not be the same as in photosynthetic activity. Kitajima et al. 2002 found no or small reductions in respiration with leaf age, Reich et al. 2009 did (Rdark/Amax constant), but often respiration declines are less steep than photosynthetic capacity, cf. Villar et al. 1995
    t_Rnight = t_Rdark * effLA_night * ctx.lookup.LookUp_Rleaf[convTnight];
#ifdef WATER
#else
    t_Rday *= effLA * 0.4; // inhibition of respiration by ca. 40%, cf. Atkin et al. 2000
#endif
}

#ifdef WATER

void Tree::CalcNPP(Context &ctx)
{
    t_NPP = 0.7 * (t_GPP - Rtot_by_Rabove * (Rtotleaf_by_Rdark * t_Rnight + t_Rstem) + (Rtotleaf_by_Rdark - 1) * t_Rnight); // t_GPP already accounts for daytime leaf respiration. would be good to go back to total respiration computation whith this new scheme (here is a succession of approximation..)
}

#else

void Tree::CalcNPP(Context &ctx)
{
    t_NPP = 0.7 * (t_GPP - Rtot_by_Rabove * (t_Rday + t_Rnight + t_Rstem)); // growth respiration, v. 2.4.1: 0.75 replaced by 0.7. According to Cannell and Thornley 2000, higher values should be used only if other sources of respiration (phloem loading etc.) are explicitly accounted for. For leaf construction, typically even a lower factor is a assumed (cf. Villar & Merino 2001, avg construction costs of 1.52 gC/gC)
    // Rleaf=Rday+Rnight is multiplied by 1.5 to also account for fine root respiration (cf as in Fyllas et al 2014 and Malhi 2012); Rstem is multiplied by 1.5 to account for coarse root respiration (according to the shoot root biomass ratio of 0.2 - Jérôme's paper in prep- and also to branch respiration (Meir & Grace 2002, Cavaleri 2006, Asao 2015).
}

#endif

// ####################################################
//  Leaf dynamics and C allocation, called by Tree::Growth
// #####################################################
//! - NPP allocation to leaves
//! - In this current scheme of leaf demography and phenology in three leaf age classes: only the old leaves generate litterfall, and the dynamic of leaves cycle is generated by the dynamic of NPP, with a total leaf biomass varying - as opposed to De Weirdt et al 2012 in ORCHIDEE, but as in Wu et al 2016 but importantly without prescribing litterfall-
void Tree::UpdateLeafDynamics(Context &ctx)
{
    float SLA = 1.0 / t_LMA; // newIM, but should be added as a tree variable to avoid multiple computation

    float ileafdem_resolution = 1.0 / float(ctx.grid.leafdem_resolution);
    float flush = 2.0 * fmaxf(t_NPP, 0.0) * ctx.params.falloccanopy * 0.68 * SLA; // 0.68 is the fraction of NPP allocated to leaves (the other 32% being allocated to twigs and fruits), retrieved from Chave et al. 2008, 2010.
    // fine resolution of flush
    float flush_fine = ileafdem_resolution * flush;
    float lambda_young = ileafdem_resolution * t_lambda_young;
    float lambda_mature = ileafdem_resolution * t_lambda_mature;
    float lambda_old = ileafdem_resolution * t_lambda_old;

    if (ctx.opt._LA_regulation == 0)
    {
        // leaf demography, without any dynamic leaf area regulation
        t_litter = 0.0;

        for (int i = 0; i < ctx.grid.leafdem_resolution; i++)
        {
            float new_litter = lambda_old * t_oldLA;
            float new_young = flush_fine;
            float new_mature = t_youngLA * lambda_young;
            float new_old = t_matureLA * lambda_mature;

            t_youngLA += new_young - new_mature;
            t_matureLA += new_mature - new_old;
            t_oldLA += new_old - new_litter;

            t_litter += new_litter * t_LMA;
        }
    }
    else
    {
        // dynamic leaf area regulation (!!!: TODO, literature references)
        // the main idea is that the fixed fractions with which allocation to leaves and stem are modelled ("ctx.params.falloccanopy", "ctx.params.fallocwood", as fractions of newly assimilated carbon) are only regulative, i.e. they apply if the tree grows in sunlight without much constraint and still has space for new leaves. If a deviation from these conditions occurs, then the allocation patterns will change. 1) Reaching maximum leaf area: If a tree has filled its whole crown with leaves and any additional tree leaf would only result in more self-shading and a net loss of carbon, then the tree will only allocate leaves for keeping up the maximum leaf area and a) allocate excess carbon to a non-structural carbon (NSC) storage b) or allocate it to growth in diameter, if the storage is full already 2) Small or no potential for photosynthesis (maximum leaf area is smaller than the tree's actual leaf area, down to 0). Then, the tree will not allocate any carbon to leaf construction and reserve carbon mostly for respiration, storage or stem growth, until a gap appears. The module as a whole prioritizes leaf upkeep over growth, which is a conversative strategy. In reality trees or species may have more aggressive strategies (i.e. growing faster when shaded), but this would be worth a research project on its own
        // first, get LAmax, i.e. the maximum of leafarea that the tree could allocate before creating too much self-shading
        float LAIexperienced_eff;
        CalcLAmax(ctx, LAIexperienced_eff, t_LAmax);

#ifdef PHENO_DROUGHT

        /* version where the multiplier modifies leafarea_max, similarly to Xu et al. 2016: drawbacks: no possibilities for simultaneous increase of leaf production and old leaves litterfal.
         if (t_phi_root < ctx.S[t_sp_lab].s_tlp) t_Ndays_dry++;
         else t_Ndays_dry=0;
         if (t_phi_root > 0.5*ctx.S[t_sp_lab].s_tlp) t_Ndays_wet++;
         else t_Ndays_wet=0;
         if (t_Ndays_dry >=5) t_pheno_factor -=0.05;
         if (t_Ndays_wet >=5) t_pheno_factor +=0.05;
         t_pheno_factor=fmaxf(0.0, fminf(1.0, t_pheno_factor));
         t_LAmax = fmaxf(0.0,t_pheno_factor*t_LAmax);
         */

        /* version where the multiplier modifies lambda_old, ie residence time in old leaf pool, more similar to Chen et al. 2020, and could be supported by a less negative TLP of old leaves than younger leaves. */
        // but still using an approach similar to Xu et al. 2016, ie. counting number of "dry days" to trigger an increase of old leaf fall
        float Thres = fminf(ctx.params.pheno_a0 * t_tlp, -0.01 * t_height - ctx.params.pheno_b0);
        if (t_phi_root < Thres)
            t_Ndays_dry++;
        else
            t_Ndays_dry = 0;
        if (t_phi_root > Thres)
            t_Ndays_wet++;
        else
            t_Ndays_wet = 0;
        if (t_Ndays_dry >= 1)
            t_pheno_factor -= ctx.params.pheno_delta;
        if (t_Ndays_wet >= 1)
            t_pheno_factor += ctx.params.pheno_delta;
        t_pheno_factor = fmaxf(0.001, fminf(1.0, t_pheno_factor));
        lambda_old = fmin(1.0, fmaxf(0.0, lambda_old / t_pheno_factor));
        t_LAmax = fmaxf(0.0, t_LAmax);

#else
        t_LAmax = fmaxf(0.0, t_LAmax);
#endif

        // carbon intended for tree growth, mobilized as if it was used for flushing leaves (in case it is needed, if not it will be backconverted at the end of the procedure), in leaf area units
        float flush_biometry = 2.0 * fmaxf(t_NPP, 0.0) * ctx.params.fallocwood * 0.6 * SLA; // new in v. 2.4.0: only 60% of woody npp is actually used for construction, the rest is for branch fall repair (cf. Malhi et al. 2011)
        float flush_biometry_fine = ileafdem_resolution * flush_biometry;

        // storage carbon, mobilized as if it was used for flushing leaves (in case it is needed, if not, backconversion at end of procedure), in leaf area units
        float flush_storage = 2.0 * t_carbon_storage * SLA;
        float flush_storage_fine = ileafdem_resolution * flush_storage;

        // We have temporarily converted all carbon into carbon for leaf flushing, so we set the pools to zero, they will be filled up again later with all the carbon that was not required. Also overall leaf litter is set to zero
        t_carbon_biometry = 0.0;
        t_carbon_storage = 0.0;
        t_litter = 0.0;

        // get the maximum leafarea the tree can currently support
        // given that we assume that the leaf density environment only changes once every iteration, this can also be calculated only once every iteration

        // float new_young = 0.0;
        for (int i = 0; i < ctx.grid.leafdem_resolution; i++)
        {
            float new_litter = lambda_old * t_oldLA;

            if (new_litter < 0)
            {
                cout << "Warning in UpdateLeafDynamics, newlitter <0, lambda_old=" << lambda_old << " t_oldLA=" << t_oldLA << " t_pheno_factor=" << t_pheno_factor << endl;
            }

            // first check whether tree can still support more leaves (does not have excess leafarea)
            // then calculate the flush needed, which is at least the amount of leaves to keep up the current leafarea

            ///* newIM 2021: modification

            float la_excess = t_LA - t_LAmax;
            float flush_needed = new_litter; // flush_needed is set equal to litter, but will be reduced, if there is excess leaves on the tree

            float new_young;
            // if leafarea does not exceed the maximum leafarea, just allocate all the flush, otherwise adjust flush_needed
            if (la_excess < 0.0)
            {
                la_excess = 0.0;
                new_young = flush_fine; // incoming flux
            }
            else
            {
                flush_needed = fmaxf(flush_needed - la_excess, 0.0); // newIM 2021
                new_young = fminf(flush_fine, flush_needed);         // incoming flux
            }

            // cout << "flush_needed=" << flush_needed << " t_LAmax=" << t_LAmax << " t_LA=" << t_LA << " new_litter=" << new_litter << " new_young=" << new_young << endl;

            // End modification newIM 2021*/

            /*new version newIM 2021

             float flush_needed=t_LAmax - t_LA + new_litter;
             float new_young= fminf( fmaxf(flush_needed, 0), flush_fine);

             //cout << "flush_needed=" << flush_needed << " t_LAmax=" << t_LAmax << " t_LA=" << t_LA << " new_litter=" << new_litter << " new_young=" << new_young << endl;

             end new version newIM 2021*/

            // subtract from flush whatever has been allocated
            flush -= new_young;

            // if flush needed is higher than incoming flush, then get more flush from carbohydrates usually reserved for diameter growth or storage
            float flush_diff = flush_needed - new_young;

            if (flush_diff > 0.0)
            {
                float flush_from_biometry = fminf(flush_biometry_fine, flush_diff); // flush_diff is the extra flush needed, so limit the amount of flush taken from biometry to what is needed
                new_young += flush_from_biometry;                                   // add extra flush to new leaves
                flush_biometry -= flush_from_biometry;                              // and deduce it from the source

                flush_diff = flush_needed - new_young;

                if (flush_diff > 0.0)
                {
                    float flush_from_storage = fminf(flush_storage_fine, flush_diff); // flush_diff is the extra flush needed, so limit the amount of flush taken from storage to what is needed
                    new_young += flush_from_storage;                                  // add extra flush to new leaves
                    flush_storage -= flush_from_storage;                              // and deduce it from the source
                }
            }

            // transferring leaves across pools
            float new_mature = t_youngLA * lambda_young;
            float new_old = t_matureLA * lambda_mature;

            t_youngLA += new_young - new_mature;
            t_matureLA += new_mature - new_old;
            t_oldLA += new_old - new_litter;

            t_litter += new_litter * t_LMA;

            if (t_litter < 0 || t_oldLA < 0)
            {
                cout << "Warning in UpdateLeafDynamics, new_litter=" << new_litter << ", lambda_old=" << lambda_old << " t_oldLA=" << t_oldLA << " t_pheno_factor=" << t_pheno_factor << " t_LMA=" << t_LMA << " new_old=" << new_old << " t_Ndays_dry=" << t_Ndays_dry << " t_Ndays_wet=" << t_Ndays_wet << endl;
            }

            t_LA = t_youngLA + t_matureLA + t_oldLA;
        }

        // This is to ensure that trees do not get infinitely small leaf areas (let's assume minimum size of 2.5cm * 2cm = 0.0005 m2)
        if (t_LA < 0.0005)
            t_LA = 0.0;

        // refill storage, if storage has been tapped
        // in light of floating point calculations not being 100% exact, we have to make sure that quantities are above zero (and not -...e-9)
        float carbon_storage_max = CalcCarbonStorageMax(ctx);
        float carbon_from_flush = 0.5 * (flush_storage + flush) * t_LMA;
        float carbon_excess = carbon_from_flush - carbon_storage_max;

        if (carbon_excess > 0.0)
        {
            t_carbon_storage = carbon_storage_max;
            if (ctx.opt._seedsadditional == 0 || t_dbh < t_dbhmature)
            {
                t_carbon_biometry += carbon_excess; // by default, excess carbon that cannot be stored is allocated to growth
            }
            else
            {
                float seedcarbon = t_NPP * ctx.params.falloccanopy * 0.08 * 0.5;
                t_multiplier_seed = int((carbon_excess + seedcarbon) / seedcarbon);
            }
        }
        else if (carbon_from_flush > 0.0)
        {
            t_carbon_storage = carbon_from_flush;
        }
        float carbon_biometry_total = 0.5 * flush_biometry * t_LMA;
        if (carbon_biometry_total > 0.0)
            t_carbon_biometry += carbon_biometry_total;
    }
}

// Compute biometric relations, including allometry
//! - New standalone function in v.2.3.0
//! - volume in ctx.params.m^3: the first factor of 2 is to convert C into biomass. the 1/s_ wsg to convert biomass into volume (g/cm^3). the 1e-6 term converts cm^3 into ctx.params.m^3 (the sole metric unit in the model). ctx.params.fallocwood is the fraction of biomass allocated to aboveground wood (stem + branches) growth. For the time being, we shall assume that a fixed proportion of NPP is allocated into AGB production. Currently, 0.20=%biomasse allocated to stem increment could be a global variable, even though this % allocation could in fact vary with resouce variation/co-limitation
void Tree::UpdateTreeBiometry(Context &ctx)
{

    //! taking into account wood elements recycling (ex. fallen branches etc...)
    //! t_ddbh = flor( volume* 4.0/( 3.0*PI*t_dbh*ctx.grid.LH*t_height*ctx.grid.LV ) )* ctx.grid.NH;

    float delta_agb;
    if (ctx.opt._LA_regulation == 0)
    {
        delta_agb = 2.0 * t_NPP * ctx.params.fallocwood * 0.6; // new in v. 2.4.0: only 60% of woody npp is actually used for construction, the rest is for branch fall repair (cf. Malhi et al. 2011)
    }
    else
    {
        delta_agb = 2.0 * t_carbon_biometry;
    }

    if (t_dbh > t_dbhmax)
        delta_agb *= fmaxf(3.0 - 2.0 * t_dbh / t_dbhmax, 0.0);

    // Tree dbh increment
    float ddbh = CalcIncrementDBH(ctx, delta_agb); // moved to a separate empirical function CalcIncrementDBH
    // With V=pi*r^2*h, increment of volume = dV = 2*pi*r*h*dr + pi*r^2*dh
    // With isometric growth assumption (ddbh/dbh=dh/h)and dbh=2*r: dV=3/4*pi*dbh*h*ddbh, ddbh in ctx.params.m, it follows: ddbh = 4/3 * V = 4/3 * 1/(pi*dbh*h)
    if (t_dbh + ddbh > 0.1 && t_dbh < 0.1)
        ctx.S[t_sp_lab].s_nbind10++;
    if (t_dbh + ddbh > 0.3 && t_dbh < 0.3)
        ctx.S[t_sp_lab].s_nbind30++;

    t_dbh += ddbh;
    UpdateSapwoodArea(ctx, ddbh);

    UpdateHeight(ctx);
    UpdateCR(ctx);
    UpdateCD(ctx);

#ifdef Output_ABC
    ctx.S[t_sp_lab].s_dbhmax_realized = fmaxf(ctx.S[t_sp_lab].s_dbhmax_realized, t_dbh);
#endif
}

void Tree::UpdateVolumeDensity()
{
    float crown_area = PI * t_CR * t_CR;
    float crown_area_nogaps = GetCrownAreaFilled(crown_area);
    t_LAI = t_LA / crown_area_nogaps;
}

// ####################################################
//  Tree death, called by Tree::Update
// ####################################################
//! - This function records basic properties of a tree that has died and empties the variables at the tree site.
//! - It does not, however, apply the destructor to the tree object
void Tree::Death(Context &ctx)
{

#ifdef TRACK_INDIVIDUALS
    if (t_timeofyear_born >= 0)
    {
        float agb = 1000.0 * CalcAGB(ctx);
        if (t_dbh >= 0.1)
        {
            ctx.out.output_track[1] << t_site << "\t" << t_timeofyear_born << "\t" << ctx.time.iter << "\t" << t_age << "\t" << t_seedsproduced_sumyear << "\t" << t_seedsproduced << "\t" << t_time_carbonstarvation_year << "\t" << t_time_carbonstarvation << "\t" << t_dbh << "\t" << t_dbh - t_dbh_tracked << "\t" << t_height << "\t" << t_height - t_height_tracked << "\t" << t_CR << "\t" << t_CR - t_CR_tracked << "\t" << agb << "\t" << agb - t_agb_tracked << "\t" << t_GPP_sumyear << "\t" << t_GPPsquared_sumyear << "\t" << t_NPP_sumyear << "\t" << t_NPPsquared_sumyear << "\t" << t_Rday_sumyear << "\t" << t_Rnight_sumyear << "\t" << t_Rstem_sumyear << "\t" << t_LAIabove_effavgyear << "\t" << t_carbon_storage_avgyear << endl;
        }

        ctx.out.output_track[2] << t_site << "\t" << t_timeofyear_born << "\t" << ctx.time.iter << "\t" << t_age << "\t" << t_seedsproduced << "\t" << t_time_carbonstarvation << "\t" << t_dbh << "\t" << t_height << "\t" << t_CR << "\t" << agb << "\t" << t_GPPcum << "\t" << t_NPPcum << "\t" << t_LAIcum << "\t" << t_LAIeffcum << "\t" << t_GPPsquared_cum << "\t" << t_NPPsquared_cum << "\t" << t_LAIsquared_cum << "\t" << t_LAIeffsquared_cum << endl;
    }
#endif

#ifdef MIP_Lichstein
    if (ctx.opt._FromInventory || (!ctx.opt._FromInventory && ctx.time.iter >= (ctx.time.nbiter - 100 * ctx.time.iterperyear)))
    {
        if (t_dbh * ctx.grid.LH >= 0.01 && t_inInventory == 1)
        {
            float agb = 0.5 * CalcAGB(ctx); // in kg C
            ctx.out.output_MIP_ind << ctx.time.iter << "\t" << ctx.S[t_sp_lab].s_name << "\t" << -9999 << "\t" << 0.0 << "\t" << 1.0 << "\t" << t_dbh * 100 << "\t" << t_height << "\t" << -9999 << "\t" << agb << "\t" << 1000 * t_wsg << "\t" << 1000 / t_LMA << "\t" << t_Nmass << "\t" << t_Pmass << "\t" << t_dbhmax << "\t" << t_tlp << "\t" << t_leafarea << endl;
            t_inInventory = 0;
        }
    }
#endif

    // new v.2.4: statistics are now calculated inside the Death() function
    // tree death statistics
    ctx.diag.nbdead_n1++;
    ctx.diag.nblivetrees--;
    if ((ctx.S[t_sp_lab].s_nbind) > 0)
        (ctx.S[t_sp_lab].s_nbind)--;
    if (t_dbh * ctx.grid.LH > 0.1)
    {
        ctx.diag.nbdead_n10++;
        if ((ctx.S[t_sp_lab].s_nbind10) > 0)
            (ctx.S[t_sp_lab].s_nbind10)--;
#ifdef Output_ABC
        int row = t_site / ctx.grid.cols;
        int col = t_site % ctx.grid.cols;
        if (row >= row_start && row < row_end && col >= col_start && col < col_end)
            nbdead_n10_abc++;
#endif
    }
    if (t_dbh * ctx.grid.LH > 0.3)
    {
        ctx.diag.nbdead_n30++;
        if ((ctx.S[t_sp_lab].s_nbind30) > 0)
            (ctx.S[t_sp_lab].s_nbind30)--;
    }
    // New v.2.2. new outputs
    if (ctx.opt._OUTPUT_extended)
    {
        if (ctx.time.iter == 2)
            ctx.out.output[23] << "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_height << "\n";
        if (ctx.time.iter == int(ctx.time.nbiter / 2))
            ctx.out.output[24] << "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_height << "\n";
        if (ctx.time.iter == int(ctx.time.nbiter - 1))
            ctx.out.output[25] << "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_height << "\n";
    }

    t_sp_lab = 0;
    t_age = 0;
    t_hurt = 0;
    t_NPP = t_GPP = t_Rday = t_Rnight = t_Rstem = 0.0; //! new v.2.3
    t_dbh = t_height = t_CR = t_CD = 0.0;
    t_CrownDisplacement = 0;

    if (ctx.opt._BASICTREEFALL)
        t_Ct = 0.0;

#ifdef Output_ABC
    t_dbh_previous = 0.0;
#endif
}

// ##################################
//  Tree death and growth
// ##################################
//! - This routine calls the appropriate DeathRate modules; if the death condition is met, function Tree::Death() is called, otherwise function Tree::Growth is called
void Tree::Update(Context &ctx)
{
    int death;
    if (t_age)
    {
        if (t_dbh > 0.1)
            ctx.diag.nbtrees_n10++;
        if (t_dbh > 0.3)
            ctx.diag.nbtrees_n30++;

#ifdef WATER
        // !!!: changed by FF, t_PPFD has been removed in v.2.5 (along with t_VPD, t_T). It was never updated, because, although formally passed on to DeathRate(), it was never used in the DeathRate() function. Furthermore, these quantities were dangerous, because they were defined at the tree level, but could change throughout the crown. Now passed by reference only where they are needed. This makes checking whether they are actually needed easier as well. Whether this affects any procedures in WATER module, needs to be checked)
        // Start new in v3.0 IM
        // Fluxh(int(t_height)+1);            //strangly, t_PPFD was not updated at the beginning of Tree::Update (as I here suggest to do now),but only in Tree::Growth and CalcRespGPP, even though t_PPFD was already used in the first part of Tree::Update (death).
        // End new in v3.0 IM
        Water_availability(ctx); // here, t_phi_root and WSF are updated, which is needed to compute Deathrate, carbon assimilation and transpiration
#endif
        // v.2.4.0: outputs have been moved to Death() function
        if (ctx.opt._NDD)
            death = int(gsl_rng_uniform(ctx.rng.gslrand) + DeathRateNDD(ctx, t_dbh, t_NPPneg, t_NDDfield[t_sp_lab]));
        else
#ifdef WATER // note that I did not include a version with both drought-induced mortality/carbon starvation and NDD effect on mortality, to be done if needed.
             // death = int(gsl_rng_uniform(ctx.rng.gslrand)+DeathRate(ctx, t_dbh, t_NPPneg, t_phi_root));
            if (ctx.opt._LA_regulation == 0)
                death = int(gsl_rng_uniform(ctx.rng.gslrand) + DeathRate(ctx, t_dbh, t_NPPneg, t_phi_root));
            else
                death = int(gsl_rng_uniform(ctx.rng.gslrand) + DeathRate(ctx, t_dbh, t_carbon_storage, t_phi_root)); // newIM 2021: directly use the t_carbon_stoarge variable instead of NPPneg
#else
            if (ctx.opt._LA_regulation == 0)
            death = int(gsl_rng_uniform(ctx.rng.gslrand) + DeathRate(ctx, t_dbh, t_NPPneg));
        else
            death = int(gsl_rng_uniform(ctx.rng.gslrand) + DeathRate(ctx, t_dbh, t_carbon_storage)); // newIM 2021: directly use the t_carbon_storage variable instead of NPPneg
#endif
        if (death)
            Death(ctx);
        else
            Growth(ctx); // v.2.4: t_hurt is now updated in the TriggerTreefallSecondary() function
    }
}

// ####################################
//  Tree falling function, called by TriggerTreefall
// ####################################
//! - Tree falling routine, formerly FallTree(),  ctx.opt._BASICTREEFALL, changed in v.2.4.0
//! - Creates a treefall (but no longer treefall probability). Takes angle as argument and can now be used for primary, secondary treefalls, forestry, or other disturbances
//! - NEW in TROLL v.2.4: FallTree() function has become Treefall() function, calculation of angle and treefall outside of function, and damages are now added up from several treefalls
void Tree::Treefall(Context &ctx, float angle)
{
    // treefall statistics
    ctx.diag.nbTreefall1++;
#ifdef Output_ABC
    if (t_dbh * ctx.grid.LH > 0.1)
    {
        ctx.diag.nbTreefall10++;
        int row = t_site / ctx.grid.cols;
        int col = t_site % ctx.grid.cols;
        if (row >= row_start && row < row_end && col >= col_start && col < col_end)
            nbTreefall10_abc++;
    }
#else
    if (t_dbh * ctx.grid.LH > 0.1)
        ctx.diag.nbTreefall10++;
#endif
    if (t_dbh * ctx.grid.LH > 0.3)
        ctx.diag.nbTreefall30++;
    int xx, yy;
    int row0, col0, h_int, r_int;
    float h_true = t_height * ctx.grid.LV;
    h_int = int(h_true * ctx.grid.NH);
    row0 = t_site / ctx.grid.cols;
    col0 = t_site % ctx.grid.cols;

    // update of ctx.field.Thurt field at the site of the tree, for consistency
    // ctx.field.Thurt[0][t_site+ctx.grid.sites] = max(int(t_height),ctx.field.Thurt[0][t_site+ctx.grid.sites]);
    // fallen stem destructs other trees
    for (int h = 1; h < h_int; h++)
    {                                                // loop on the fallen stem (horizontally)
        xx = int(fmaxf(col0 + h * cos(angle), 0.0)); // get projection in col (= xx) direction, where xx is absolute location
        if (xx < ctx.grid.cols)
        {
            yy = int(row0 + h * sin(angle)); // get projection in row (= yy) direction, where yy is absolute location
            ctx.field.Thurt[0][xx + (yy + ctx.grid.rows) * ctx.grid.cols] = max(int(t_height), int(ctx.field.Thurt[0][xx + (yy + ctx.grid.rows) * ctx.grid.cols]));
            // ctx.field.Thurt[0] where the stem fell, calculation: xx+(yy+ctx.grid.rows)*ctx.grid.cols= xx + yy*ctx.grid.cols + ctx.grid.rows*ctx.grid.cols = xx + yy*ctx.grid.cols + ctx.grid.sites / NEW in v.2.4: addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }

    // fallen crown destructs other trees, less damaging than stem
    xx = col0 + int((h_true * ctx.grid.NH - t_CR) * cos(angle));
    yy = row0 + int((h_true * ctx.grid.NH - t_CR) * sin(angle));
    r_int = int(t_CR);
    for (int col = max(0, xx - r_int); col < min(ctx.grid.cols, xx + r_int + 1); col++)
    { // loop on the fallen crown (horizontally)
        for (int row = max(0, yy - r_int); row < min(ctx.grid.rows, yy + r_int + 1); row++)
        {
            if ((col - xx) * (col - xx) + (row - yy) * (row - yy) < r_int * r_int)
                ctx.field.Thurt[0][col + (row + ctx.grid.rows) * ctx.grid.cols] = max(int((t_height - t_CR * ctx.grid.NV * ctx.grid.LH) * 0.5), int(ctx.field.Thurt[0][col + (row + ctx.grid.rows) * ctx.grid.cols])); // less severe damage than stem / NEW in v.2.4: max() or addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }
    // v.2.4.0: outputs have been moved to Death() function
    Death(ctx);
}

// ####################################################
//  Computes Average and OutputField
// ####################################################
//  - Short routine that basically only updates the vector s_output_field
void Tree::Average(Context &ctx)
{
    if (t_age > 0)
    {
        if (t_dbh * ctx.grid.LH >= 0.1)
        {
            (ctx.S[t_sp_lab].s_sum10)++;
            ctx.S[t_sp_lab].s_ba10 += t_dbh * ctx.grid.LH * t_dbh * ctx.grid.LH * 3.1415 * 0.25;
        }
        if (t_dbh * ctx.grid.LH >= 0.3)
            (ctx.S[t_sp_lab].s_sum30)++;
        ctx.S[t_sp_lab].s_ba += t_dbh * ctx.grid.LH * t_dbh * ctx.grid.LH * 3.1415 * 0.25;
        ctx.S[t_sp_lab].s_npp += t_NPP * 1.0e-6;
        ctx.S[t_sp_lab].s_gpp += t_GPP * 1.0e-6;
        float agb = CalcAGB(ctx);
        ctx.S[t_sp_lab].s_agb += agb;
        ctx.S[t_sp_lab].s_rday += t_Rday * 1.0e-6;
        ctx.S[t_sp_lab].s_rnight += t_Rnight * 1.0e-6;
        ctx.S[t_sp_lab].s_rstem += t_Rstem * 1.0e-6;
        ctx.S[t_sp_lab].s_litterfall += t_litter * 1.0e-6;

#ifdef WATER
        int crown_top = int(t_height);
        int crown_base = int(t_height - t_CD);
        float grad = 1 / float(crown_top - crown_base + 1);
        for (int l = crown_base; l < (crown_top + 1); l++)
        {
            ctx.soil.LAI_young[l] += t_youngLA * grad;
            ctx.soil.LAI_mature[l] += t_matureLA * grad;
            ctx.soil.LAI_old[l] += t_oldLA * grad;
        }

        ctx.soil.abund_phi_root += t_phi_root;
        if (t_dbh * ctx.grid.LH >= 0.1)
            ctx.soil.abund10_phi_root += t_phi_root;
        ctx.soil.agb_phi_root += agb * t_phi_root;
#endif

#ifdef MIP_Lichstein
        if (ctx.time.iter % ctx.time.iterperyear == 364 && (ctx.opt._FromInventory || (!ctx.opt._FromInventory && ctx.time.iter >= (ctx.time.nbiter - 100 * ctx.time.iterperyear))))
        {
            if (t_dbh * ctx.grid.LH >= 0.01)
            {
                t_inInventory = 1;
                ctx.out.output_MIP_ind << ctx.time.iter << "\t" << ctx.S[t_sp_lab].s_name << "\t" << -9999 << "\t" << 1.0 << "\t" << 0.0 << "\t" << t_dbh * 100 << "\t" << t_height << "\t" << -9999 << "\t" << 0.5 * agb << "\t" << 1000 * t_wsg << "\t" << 1000 / t_LMA << "\t" << t_Nmass << "\t" << t_Pmass << "\t" << t_dbhmax << "\t" << t_tlp << "\t" << t_leafarea << endl;
            }
        }
#endif
    }
}

// Computation of dbh histograms
void Tree::histdbh(Context &ctx)
{
    if (t_age)
        ctx.diag.nbdbh[int(100. * t_dbh * ctx.grid.LH)]++;
    // where dbh is in cm (it is in number of horizontal cells throughout the code)
    // values are always rounded down (so ctx.diag.nbdbh[30] gives you trees with more than 30 cm dbh, and less than 31))
}

#ifdef WATER
// Standard outputs during the simulation -- written to file
void Tree::OutputTreeStandard(Context &ctx, fstream &output)
{
    output << ctx.time.iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_height << "\t" << t_dbh << "\t" << t_litter << "\t" << t_age << "\t" << t_LA << "\t" << t_youngLA << "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_CR << "\t" << t_CD << "\t" << t_GPP << "\t" << t_NPP << "\t" << t_Rstem << "\t" << t_Rnight << "\t" << ctx.field.LAI3D[int(t_height)][t_site + ctx.grid.SBORD] << "\t" << ctx.field.LAI3D[int(t_height - t_CD) + 1][t_site + ctx.grid.SBORD] << "\t" << t_root_depth << "\t" << t_phi_root << "\t" << t_WSF << "\t" << t_WSF_A << "\t" << t_transpiration << "\t" << t_LAImax << "\t" << t_LAmax;
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        output << "\t" << t_root_biomass[l];
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        output << "\t" << t_soil_layer_weight[l];
    output << endl;
}
// Standard outputs during the simulation -- written to screen in real time
void Tree::OutputTreeStandard(Context &ctx)
{
    cout << ctx.time.iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_height << "\t" << t_dbh << "\t" << t_litter << "\t" << t_age << "\t" << t_LA << "\t" << t_youngLA << "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_CR << "\t" << t_CD << "\t" << t_GPP << "\t" << t_NPP << "\t" << t_Rstem << "\t" << t_Rday << "\t" << t_Rnight << "\t" << ctx.field.LAI3D[int(t_height)][t_site + ctx.grid.SBORD] << "\t" << ctx.field.LAI3D[int(t_height - t_CD) + 1][t_site + ctx.grid.SBORD] << "\t" << t_root_depth << "\t" << t_phi_root << "\t" << t_WSF;
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        cout << "\t" << t_root_biomass[l];
    for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        cout << "\t" << t_soil_layer_weight[l];
    cout << endl;
}
#endif

// Calculate the crown area filled by leaves (only relevant for crown gap fractions > 0.0)
float Tree::GetCrownAreaFilled(float crown_area)
{
    // for now calculated explicitly. Ideally, we would replace the loop by an equation that expresses the underlying logic
    int crown_intarea = int(crown_area);      // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea, 1);    // minimum area of crown (1)
    crown_intarea = min(crown_intarea, 1963); // maximum area of crown (radius 25), int(3.14*25*25)

    int crown_intarea_gaps = 0;
    float fraction_filled_target = t_fraction_filled;
    float fraction_filled_actual = 0.0;

    for (int i = 0; i < crown_intarea; i++)
    {
        if (fraction_filled_actual > fraction_filled_target)
        {
            fraction_filled_actual = (fraction_filled_actual * float(i)) / (float(i) + 1.0);
            crown_intarea_gaps++;
        }
        else
            fraction_filled_actual = (fraction_filled_actual * float(i) + 1.0) / (float(i) + 1.0);
    }

    // now determine crown_area_filled, depending on whether the next voxel is filled or not filled
    float crown_area_filled;
    if (fraction_filled_actual > fraction_filled_target)
    {
        crown_area_filled = float(crown_intarea - crown_intarea_gaps);
    }
    else
        crown_area_filled = crown_area - float(crown_intarea_gaps);

    return (crown_area_filled);
}

#ifdef TRACK_INDIVIDUALS
// Diagnostic function to track trees born at a reference year
float Tree::StartTracking(Context &ctx)
{
    // Only tracks trees born in a mature forest at year 501
    // currently hardcoded
    if (ctx.time.iter >= 6000 && ctx.time.iter < 6012)
        t_timeofyear_born = ctx.time.iter % ctx.time.iterperyear;
    else
        t_timeofyear_born = -1;

    if (t_timeofyear_born >= 0)
    {
        // initialise variables
        t_seedsproduced = 0;
        t_seedsproduced_sumyear = 0;
        t_time_carbonstarvation = 0;
        t_time_carbonstarvation_year = 0;

        t_GPP_sumyear = 0.0;
        t_NPP_sumyear = 0.0;
        t_GPPsquared_sumyear = 0.0; // for standard deviation
        t_NPPsquared_sumyear = 0.0; // for standard deviation
        t_Rday_sumyear = 0.0;
        t_Rnight_sumyear = 0.0;
        t_Rstem_sumyear = 0.0;
        t_LAIabove_effavgyear = 0.0;
        t_carbon_storage_avgyear = 0.0;

        t_LAIcum = 0.0;
        t_LAIeffcum = 0.0;
        t_GPPcum = 0.0;
        t_NPPcum = 0.0;
        t_LAIsquared_cum = 0.0;
        t_LAIeffsquared_cum = 0.0;
        t_GPPsquared_cum = 0.0;
        t_NPPsquared_cum = 0.0;

        t_dbh_tracked = t_dbh;
        t_height_tracked = t_height;
        t_CR_tracked = t_CR;
        t_agb_tracked = 1000.0 * CalcAGB(ctx);

        ctx.out.output_track[0] << t_site << "\t" << t_timeofyear_born << "\t" << t_site % ctx.grid.cols << "\t" << t_site / ctx.grid.cols << "\t" << t_s->s_name << "\t" << t_dbh << "\t" << t_CR << "\t" << t_height << "\t" << t_agb_tracked << "\t" << t_mult_CR << "\t" << t_mult_height << "\t" << t_wsg << "\t" << t_Nmass << "\t" << t_Pmass << "\t" << t_LMA << "\t" << t_dev_wsg << "\t" << t_mult_N << "\t" << t_mult_P << "\t" << t_mult_LMA << "\t" << t_Vcmax << "\t" << t_Jmax << "\t" << t_Rdark << "\t" << t_LAImax << "\t" << t_leaflifespan << endl;
    }
}
#endif

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
#ifdef MPI
    // Treefall field passed to the n.n. procs
    MPI_ShareTreefall(ctx.field.Thurt, ctx.grid.sites);
#endif
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        // Update of Field hurt
        if (ctx.T[site].t_age)
        {
            ctx.T[site].t_hurt = max(ctx.field.Thurt[0][site + ctx.grid.sites], ctx.T[site].t_hurt); // NEW in v.2.4: addition of damages, alternative: max()
#ifdef MPI
            if (mpi_rank)
                ctx.T[site].t_hurt = max(ctx.T[site].t_hurt, ctx.field.Thurt[1][site]); // ? v.2.4: Update needed, ctx.field.Thurt[1], why max?
            if (mpi_rank < mpi_size - 1)
                ctx.T[site].t_hurt = max(ctx.T[site].t_hurt, ctx.field.Thurt[2][site]);
#endif
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

#ifdef MPI
    //! Treefall field passed to the n.n. procs
    MPI_ShareTreefall(ctx.field.Thurt, ctx.grid.sites);
#endif
}

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

//! Global functions to compute Vcmax (on a mass basis) and leaf dark respiration (on a leaf area basis)
// newIM : to avoid repeating same functions in Tree and Species class. Was needed for the Species class to correctly initialise Rdark and compute LCP

float CalcVcmaxm(float LMA, float Nmass, float Pmass)
{
    float SLA = 10000.0 / LMA;                                                                                                                               // in cm2 g-1
    float Vcmaxm = pow(10.0, fminf((-1.56 + 0.43 * log10(Nmass * 1000.0) + 0.37 * log10(SLA)), (-0.80 + 0.45 * log10(Pmass * 1000.0) + 0.25 * log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA). Nmass and Pmass are given in g g-1, but should be in mg g-1 in equ 2 in Domingues et al. 2010, hence the mutiplication by 1000.
    return (Vcmaxm);                                                                                                                                         //  in micromol C g-1 s-1
}
float CalcRdark(float LMA, float Nmass, float Pmass, float Vcmax)
{
    float Parea = Pmass * LMA;                                                                         // in g ctx.params.m-2
    float Narea = Nmass * LMA;                                                                         // in g ctx.params.m-2
    float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2)); // in micromolC ctx.params.m-2 s-1 //since v.2.5: correction of Atkin et al. 2015 New phytologist formula. The original formula was based on mean-centered values (cf. Atkin et al. 2015, and the clarification/correction published afterwards), unfortunately only one formula (PFT-specific formula) with absolute values is provided and this is the one used here, cf. corrigendum TableS4): https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.13253&file=nph13253-sup-0001-SupInfo.pdf Warning: Vcmax should be provided on an area basis.
    return (Rdark);
}

#ifdef CROWN_UMBRELLA

// Global function: linear decrease of crown radius
float GetRadiusSlope(Context &ctx, float CR, float crown_extent, float crown_position)
{
    float crown_slope = CR * (1.0 - ctx.crown.shape_crown) / crown_extent;
    float radius = CR - crown_slope * float(crown_position);
    return (radius);
}

// Global function: not currently used, but returns the input radius
float GetRadiusCylinder(float CR, float crown_extent, float crown_position)
{
    return (CR);
}

// Global function: converts floating point crown area into integer value, imposing lower and upper limits
int GetCrownIntarea(float crown_radius)
{
    // crown area
    float crown_area = PI * crown_radius * crown_radius;
    int crown_intarea = int(crown_area);      // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea, 1);    // minimum area of crown (1)
    crown_intarea = min(crown_intarea, 1963); // maximum area of crown (radius 25), int(3.14*25*25)
    return (crown_intarea);
}

// Global function: deduces within-crown densities from LAI with a gradient from 50% in top ctx.diag.layer to 25% in belowtop and 25% in all shells underneath (1 ctx.diag.layer for umbrella-like shape)
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens_layer) // RENAMED: ctx.params.dens → dens_layer (was shadowing global ctx.params.dens/ctx.params.ctx.params.dens)
{
    if (CD < 2.0)
    {
        dens_top = dens_belowtop = dens_layer = LAI / CD;
    }
    else if (CD < 3.0)
    {
        dens_top = 0.5 * LAI;
        dens_belowtop = dens_layer = 0.5 * LAI / (CD - 1.0);
    }
    else
    {
        dens_top = 0.5 * LAI;
        dens_belowtop = 0.25 * LAI;
#ifdef CROWN_UMBRELLA
        dens_layer = 0.25 * LAI;
#else
        dens_layer = 0.25 * LAI / (CD - 2.0);
#endif
    }
}

// Global function: deduces within-crown density from LAI, assuming uniform leaf distribution
void GetDensityUniform(float LAI, float CD, float &dens_layer) // RENAMED: ctx.params.dens → dens_layer
{
#ifdef CROWN_UMBRELLA
    float crownshells_limit = fminf(CD, 3.0);
    dens_layer = LAI / crownshells_limit;
#else
    dens_layer = LAI / CD;
#endif
}

// Global function: dummy function when no modification is needed
void KeepFloatAsIs(float input, float &output, float CD, float height, int layer_fromtop)
{
    output = input;
}

void KeepIntAsIs(int input, int &output, float CD, float height, int layer_fromtop)
{
    output = input;
}

#endif //CROWN_UMBRELLA

#ifdef G0

// Solves the quadratic equation. ROOT: finds smaller (-1), larger (1) root.
float QUAD(float A, float B, float C, int ROOT)
{
    float quad = 0.0;
    float discriminant = B * B - 4.0 * A * C;
    // ROOT=-1: smaller root; ROOT=1: larger root
    if (discriminant < 0)
    {
        cerr << "IMAGINARY ROOTS IN QUADRATIC" << endl;
        quad = 0.0;
    }
    if (A == 0.0)
    {
        if (B == 0.0)
            quad = 0.0;
        else
            quad = -C / B;
    }
    else
        quad = (-B + ROOT * sqrt(discriminant)) / (2.0 * A);
    return quad;
}

#endif




