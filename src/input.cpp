#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"
#include "internal/param_registry.hpp"
#include <functional>

//! Global function: This function provides limits and defaults for species-specific parameters
// ============================================================================
//   Species parameter registry (replaces AssignSpeciesParam if/else chain)
// ============================================================================
static std::unordered_map<std::string, std::function<void(Species &, std::string, std::string)>> species_param_registry;

//! Register the species-file columns (name, type, allowed range) used to read species parameters
static void BuildSpeciesRegistry()
{
    species_param_registry["s_name"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_name, std::string("indet_indet"), true);
    };
    species_param_registry["s_LMA"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_LMA, 0.0f, 1000.0f, 100.0f, true);
    };
    species_param_registry["s_Nmass"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_Nmass, 0.0f, 1.0f, 0.02f, true);
    };
    species_param_registry["s_Pmass"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_Pmass, 0.0f, 1.0f, 0.0005f, true);
    };
    species_param_registry["s_wsg"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_wsg, 0.0f, 1.5f, 0.6f, true);
    };
    species_param_registry["s_dbhmax"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_dbhmax, 0.0f, 2.5f, 0.5f, true);
    };
    species_param_registry["s_hmax"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_hmax, 0.0f, 100.0f, 50.0f, true);
    };
    species_param_registry["s_ah"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_ah, 0.0f, 10.0f, 0.3f, true);
    };
    species_param_registry["s_seedmass"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_seedmass, 0.0f, 10000.0f, 1.0f, true);
    };
    species_param_registry["s_regionalfreq"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_regionalfreq, 0.0f, 1.0f, 1.0f, true);
    };
    species_param_registry["s_tlp"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_tlp, -10.0f, 0.0f, -2.0f, true);
    };
    species_param_registry["s_leafarea"] = [](Species &S, std::string n, std::string v) {
        SetParameter(n, v, S.s_leafarea, 0.0f, 800.0f, 80.0f, true);
    };
}

//! @param S [in,out] species whose parameter is set
//! @param name column name in the species file
//! @param value value read for this species (text)
void AssignSpeciesParam(Species &S, const std::string &name, const std::string &value)
{
    static bool built = false;
    if (!built) { BuildSpeciesRegistry(); built = true; }

    auto it = species_param_registry.find(name);
    if (it == species_param_registry.end())
    {
        std::cerr << "Warning: Unknown species parameter '" << name << "'\n";
        return;
    }
    it->second(S, name, value);
}

// added v.3.1.6
//! @param ctx simulation context (all model state)
//! @param parameter_name name of the point-cloud parameter
//! @param parameter_value value read from the point-cloud parameter file (text)
void AssignValuePointcloud(Context &ctx, string parameter_name, string parameter_value)
{
    // we set parameters to values that have been read, or to their defaults, if outside of range or not the right type
    bool quiet = 1; //! only applies to successful initialization, warnings are always given

    if (parameter_name == "ctx.pc.mean_beam_pc")
    {
        SetParameter(parameter_name, parameter_value, ctx.pc.mean_beam_pc, 0.1f, 100.0f, 10.0f, quiet);
    }
    else if (parameter_name == "ctx.pc.sd_beam_pc")
    {
        SetParameter(parameter_name, parameter_value, ctx.pc.sd_beam_pc, 0.0f, 100.0f, 10.0f, quiet);
    }
    else if (parameter_name == "ctx.pc.klaser_pc")
    {
        SetParameter(parameter_name, parameter_value, ctx.pc.klaser_pc, 0.1f, 0.9f, 0.5f, quiet);
    }
    else if (parameter_name == "ctx.pc.transmittance_laser")
    {
        SetParameter(parameter_name, parameter_value, ctx.pc.transmittance_laser, 0.0f, 1.0f, 0.4f, quiet);
    }
    else if (parameter_name == "ctx.pc.iter_pointcloud_generation")
    {
        SetParameter(parameter_name, parameter_value, ctx.pc.iter_pointcloud_generation, 0, ctx.time.nbiter - 1, ctx.time.nbiter - 1, quiet);
    }
}

//! Global function: This function reads inputs from the general input file
//! - Also sets up initial values for some parameters computed from input parameters
//! - !!!: TO IMPLEMENT: automatic stop of program under error to avoid running empty simulations on cluster etc.
//! - !!!: TO IMPLEMENT: separate error messages (i.e. no input file provided, empty input file, etc.)
//! @param ctx simulation context (all model state)
void ReadInputGeneral(Context &ctx)
{
    fstream In(ctx.fileio.inputfile, ios::in);
    if (In)
    {
        cout << endl
             << "Reading in file: " << ctx.fileio.inputfile << endl;

        fstream In(ctx.fileio.inputfile, ios::in);
        if (In)
        {
            In.getline(ctx.fileio.buffer, 256, '\n'); // skip header
            std::string parameter_name, parameter_value;

            while (In >> parameter_name >> parameter_value)
            {
                In.getline(ctx.fileio.buffer, 256, '\n'); // consume rest of line
                AssignParamFromRegistry(parameter_name, parameter_value);
            }
        }
        else
        {
            cout << "ERROR. General input file could not be read." << endl;
        }
        // update derived parameters
        ctx.grid.sites = ctx.grid.rows * ctx.grid.cols;    // leave global for now
        ctx.grid.sites_per_dcell = ctx.grid.length_dcell * ctx.grid.length_dcell;
        ctx.grid.nbdcells = int(ctx.grid.sites / ctx.grid.sites_per_dcell);
        ctx.grid.linear_nb_dcells = int(ctx.grid.cols / ctx.grid.length_dcell);
        cout << "ctx.grid.rows: " << ctx.grid.rows << " ctx.grid.cols: " << ctx.grid.cols << " ctx.grid.HEIGHT: " << ctx.grid.HEIGHT << endl;
        cout << "Number of dcells: " << ctx.grid.nbdcells << endl;
        cout << "Lin number of dcells: " << ctx.grid.linear_nb_dcells << endl;
        cout << ctx.opt._WATER_RETENTION_CURVE << " " << ctx.opt._SOIL_LAYER_WEIGHT << " " << ctx.params.Cair << endl;

#ifdef WATER
        ctx.grid.i_sites_per_dcell = 1.0 / float(ctx.grid.sites_per_dcell);
        ctx.params.PPFDtoSW = 1 / ctx.params.SWtoPPFD;
#ifdef G0
        ctx.params.g0 *= 0.001 * GSVGSC; // ctx.params.g0 is provided in mmolH20 ctx.params.m-2 s-1, but will be used in Photosyn function for stomatal conductance to CO2 in molCO2 ctx.params.m-2 s-1 (and then converted back to flux of H20 in FluxesLeaf)
#endif
#endif
        ctx.grid.LV = 1.0 / ctx.grid.NV;
        ctx.grid.LH = 1.0 / ctx.grid.NH;
        if (ctx.time.nbout)
            ctx.time.freqout = max(1, ctx.time.nbiter / ctx.time.nbout); // at least 1: runs shorter than nbout iterations crashed on "iter % 0"

        ctx.params.kpar = ctx.params.klight * ctx.params.absorptance_leaves; // ctx.params.kpar is the ctx.params.klight factor times the absorptance of leaves
        // convert correlations to covariances
        ctx.intra.cov_N_P = ctx.intra.corr_N_P * ctx.intra.sigma_N * ctx.intra.sigma_P;
        ctx.intra.cov_N_LMA = ctx.intra.corr_N_LMA * ctx.intra.sigma_N * ctx.intra.sigma_LMA;
        ctx.intra.cov_P_LMA = ctx.intra.corr_P_LMA * ctx.intra.sigma_P * ctx.intra.sigma_LMA;

        if (ctx.intra.cov_N_P == 0.0 || ctx.intra.cov_N_LMA == 0.0 || ctx.intra.cov_P_LMA == 0.0)
        {
            cerr << "\nCovariance matrix N,P,LMA could not be decomposed. Using uncorrelated variation of trait values instead" << endl;
            ctx.rng.covariance_status = 0;
        }
        else
        {
            cout << "Correlation status. ctx.intra.corr_N_P: " << ctx.intra.corr_N_P << " ctx.intra.cov_N_LMA: " << ctx.intra.corr_N_LMA << " ctx.intra.corr_P_LMA: " << ctx.intra.corr_P_LMA << endl;
            ctx.rng.covariance_status = 1;
            // Initialise covariance matrix for N, P, LMA
            ctx.rng.mcov_N_P_LMA = gsl_matrix_alloc(3, 3);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 0, 0, ctx.intra.sigma_N * ctx.intra.sigma_N);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 0, 1, ctx.intra.cov_N_P);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 0, 2, ctx.intra.cov_N_LMA);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 1, 0, ctx.intra.cov_N_P);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 1, 1, ctx.intra.sigma_P * ctx.intra.sigma_P);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 1, 2, ctx.intra.cov_P_LMA);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 2, 0, ctx.intra.cov_N_LMA);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 2, 1, ctx.intra.cov_P_LMA);
            gsl_matrix_set(ctx.rng.mcov_N_P_LMA, 2, 2, ctx.intra.sigma_LMA * ctx.intra.sigma_LMA);

            // Cholesky decomposition for multivariate draw
            cout << "\nCovariance matrix N,P,LMA: " << endl;
            for (int mrow = 0; mrow < 3; mrow++)
                cout << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 2) << endl;
            gsl_linalg_cholesky_decomp1(ctx.rng.mcov_N_P_LMA);
            cout << "\nCovariance matrix N,P,LMA (after Cholesky decomposition) " << endl;
            for (int mrow = 0; mrow < 3; mrow++)
                cout << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(ctx.rng.mcov_N_P_LMA, mrow, 2) << endl;
            // allocating mean and result vectors for multivariate draw
            ctx.rng.mu_N_P_LMA = gsl_vector_alloc(3);
            for (int j = 0; j < 3; j++)
                gsl_vector_set(ctx.rng.mu_N_P_LMA, j, 0.0); // zero means
            ctx.rng.variation_N_P_LMA = gsl_vector_alloc(3);
        }
        ctx.crown.crown_gap_fraction = fmaxf(ctx.crown.crown_gap_fraction, 0.000001); // ctx.crown.crown_gap_fraction is prevented from becoming zero in order to avoid division by zero. Given that crown area is currently limited to 1963 (int(3.14 * 25.0 * 25.0), the lowest ctx.crown.crown_gap_fraction that could potentially have an effect would be 1/1963, which is ~ 0.0005
        ctx.params.iCair = 1.0 / ctx.params.Cair;
        ctx.params.DBH0 *= ctx.grid.NH;
        ctx.params.H0 *= ctx.grid.NV;
        ctx.params.CR_min *= ctx.grid.NH;
        ctx.params.CD0 *= ctx.grid.NV;
#ifdef WATER
        // ctx.params.alpha = ctx.params.absorptance_leaves*(1-LSQ)*0.5; // apparent quantum yield to electron transport. 0.5 because each photosystem absorbs half of the photons (von Caemmerer 2000, p 35). Note that, ctx.params.alpha was previously computed as 4*the apparent quantum yield for C fixation, ctx.params.phi in molC/mol photons, with ctx.params.phi=0.06 (in Marechaux & Chave 2017) or 0.075 (in FFischer's version), which provided ctx.params.alpha=0.3 and 0.24 respectively. This now gives ctx.params.alpha=0.37. Duursma et al. 2015 used ctx.params.alpha=0.24, while values reported or used by von Caemmerer 2000, Medlyn et al. 2002 are typically 0.3-0.37. Domingues et al. 2014 Plant Ecology & Diversity reported values for ctx.params.phi ranging from 0.041 to 0.098 (ie ctx.params.alpha ranging from 0.164 to 0.392, within a same tropical forest community (Tapajos, Brazil).
        ctx.params.alpha = (1 - LSQ) * 0.5; // no need of ctx.params.absorptance_leaves since PPFD provided in argument to Photosyn is already the absorbed flux
#else
        ctx.params.alpha = 4.0 * ctx.params.phi; // apparent quantum yield to electron transport in mol e-/mol photons see Mercado et al 2009 , the conversion of the apparent quantum yield in micromolCO2/micromol quantum into micromol e-/micxromol quantum is done by multipliyng by 4, since four electrons are needed to regenerate RuBP. ctx.params.alpha is fixed at 0.3 mol e-/mol photons in Medlyn et al 2002, but see equ8 and Appendix 1 in Farquahr et al 1980: it seems that ctx.params.alpha should vary with leaf thickness: there is a fraction of incident light which is lost by absorption by other leaf parts than the chloroplast lamellae, and this fraction f may increase with leaf thickness. With the values of the paper: ctx.params.alpha= 0.5*(1-f)=0.5*(1-0.23)=0.385, but this is a theoretical value and observations often report lower values (see ex discussion in medlyn et al 2005 Tree phsyiology, Lore Veeryckt values, Mercado et al 2009 Table 10, Domingues et al. 2014)
#endif

        // new in v.3.1.2: visual extent
        int maxextent_visual = min(ctx.grid.rows, ctx.grid.cols);
        if (ctx.crown.extent_visual > maxextent_visual)
            ctx.crown.extent_visual = maxextent_visual; // make sure visualization does not exceed boundaries of simulation
        if (ctx.crown.extent_visual > 0)
        {
            // define boundaries
            int rowextent_slice = 10;
            ctx.crown.mincol_visual = ctx.grid.cols / 2 - ctx.crown.extent_visual / 2;
            ctx.crown.maxcol_visual = ctx.grid.cols / 2 + ctx.crown.extent_visual / 2;
            ctx.crown.minrow_visual = ctx.grid.rows / 2 - ctx.crown.extent_visual / 2;
            ctx.crown.maxrow_visual = ctx.grid.rows / 2 + ctx.crown.extent_visual / 2;
            ctx.crown.minrow_visual_slice = max(ctx.grid.rows / 2 - rowextent_slice / 2, ctx.crown.minrow_visual);
            ctx.crown.maxrow_visual_slice = min(ctx.grid.rows / 2 + rowextent_slice / 2, ctx.crown.maxrow_visual);
        }
    }
    else
    {
        cout << "ERROR. General input file could not be read." << endl;
    }
}

//! Global function: This function reads inputs from the species input file
//! @param ctx simulation context (all model state)
void ReadInputSpecies(Context &ctx)
{
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_species << endl;
    fstream InSpecies(ctx.fileio.inputfile_species, ios::in);
    if (InSpecies)
    {
        // possible parameters to initialise vector<string> parameter_names{"s_name","s_LMA","s_Nmass","s_Pmass","s_wsg","s_dbhmax","s_hmax","s_ah","s_seedmass","s_regionalfreq","s_tlp","s_leafarea"};
        //        int nb_parameters = int(parameter_names.size()); only works from C++11 onwards
        string parameter_names[12] = {"s_name", "s_LMA", "s_Nmass", "s_Pmass", "s_wsg", "s_dbhmax", "s_hmax", "s_ah", "s_seedmass", "s_regionalfreq", "s_tlp", "s_leafarea"};
        int nb_parameters = 12;

        // first get parameter names
        string line;
        getline(InSpecies, line);
        istringstream firstlinestream(line);
        string parameter_name;
        vector<int> parameter_index;
        int nb_parameters_fromfile = 0;

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
        }

        int nb_parameterlines = 0;
        Species species_dummy; // !!!: currently, species at 0 is not initialized, for legacy reasons, maybe to change
        ctx.S.push_back(species_dummy);

        // we go through all lines in the input file
        while (getline(InSpecies, line))
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
            Species species_new;
            // now we assign values
            for (int i = 0; i < nb_parameters; i++)
            {
                AssignSpeciesParam(species_new, parameter_names[i], parameter_values[i]);
            }
            ctx.S.push_back(species_new);
            nb_parameterlines++;
        }

        if (nb_parameterlines > 0)
        {
            ctx.grid.nbspp = nb_parameterlines;
            cout << "Successfully initialised " << ctx.grid.nbspp << " species from file." << endl;
        }
        else
        {
            ctx.grid.nbspp = 1;
            cout << "WARNING! Species file was empty. A default species is initialized." << endl;
            vector<string> parameter_values(nb_parameters, "");
            Species species_default;
            for (int i = 0; i < nb_parameters; i++)
            {
                AssignSpeciesParam(species_default, parameter_names[i], parameter_values[i]);
            }
            ctx.S.push_back(species_default);
        }

        for (int sp = 1; sp <= ctx.grid.nbspp; sp++)
        {

            ctx.S[sp].Init(ctx);
            // cout << ctx.S[sp].s_name << " LCP: " << ctx.S[sp].s_LCP << endl;
        }
        cout << "Successfully read in species file." << endl;
    }
    else
    {
        cout << "ERROR with the species file" << endl;
    }

    InSpecies.close();
}

#ifdef FULL_CLIMATE

//! Global function: This function reads inputs from the environmental daily variation input file
//! @param ctx simulation context (all model state)
void ReadInputDailyvar(Context &ctx)
{
    // currently very simple reading in, only basic error checking
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_daytimevar << endl;

    fstream InDaily(ctx.fileio.inputfile_daytimevar, ios::in);

    if (InDaily)
    {
        InDaily.getline(ctx.fileio.buffer, 256, '\n'); // read in header

        // we go through all lines in the input file
        ctx.time.nbsteps_varday = 0;
        ctx.time.nbhours_covered = 12.0; // quick fixed, but I should provide both endtime and startime for each half-hour iteration.

        float day_current, starttime_current, varday_light_current, varday_vpd_current, varday_T_current, varday_WS_current;

        string line;
        while (getline(InDaily, line))
        {
            istringstream linestream(line);

            linestream >> day_current >> starttime_current >> varday_T_current >> varday_light_current >> varday_vpd_current >> varday_WS_current;
            // float nbhours_current;
            // if(endtime_current > starttime_current) nbhours_current = endtime_current - starttime_current;  // just calculate the time window
            // else nbhours_current = 24.0 - (starttime_current - endtime_current); // i.e. if starttime is 23.5 and endtime is 0.5, etc. Calculate the inverse time window and then subtract from 24
            //  else nbhours_current =0.5;

            if (varday_vpd_current <= 0.01)
            {
                varday_vpd_current = 0.01;
                cout << "Warning, VPD is <=0.01, on day " << day_current << " at " << starttime_current << "h." << endl;
            }

            if (varday_light_current <= 1)
            {
                varday_light_current = 1;
                cout << "Warning, light is <=1, on day " << day_current << " at " << starttime_current << "h." << endl;
            }

            if (varday_WS_current <= 0.1)
            {
                varday_WS_current = 0.1;
                cout << "Warning, wind is <=0.1, on day " << day_current << " at " << starttime_current << "h." << endl;
            }

            // ctx.time.nbhours_covered += nbhours_current;
            ctx.climate.varday_light.push_back(varday_light_current);
            ctx.climate.varday_vpd.push_back(varday_vpd_current);
            ctx.climate.varday_T.push_back(varday_T_current);
            ctx.climate.varday_WS.push_back(varday_WS_current);
            // cout << "at: " << starttime_current << " vardaytime_light: " << ctx.climate.varday_light[ctx.time.nbsteps_varday] << " vardaytime_vpd: " << ctx.climate.varday_vpd[ctx.time.nbsteps_varday] << " vardaytime_T: "<< ctx.climate.varday_T[ctx.time.nbsteps_varday] << " vardaytime_WS: "<< ctx.climate.varday_WS[ctx.time.nbsteps_varday] << endl;
            ctx.time.nbsteps_varday++;
        }

        ctx.time.nbsteps_varday /= ctx.time.nbdays;

        cout << "Read in: " << ctx.time.nbsteps_varday << " timesteps per day, covering " << ctx.time.nbhours_covered << " hours of the day." << endl;
        cout << "Successfully read in daytime variation file" << endl;
        ctx.time.inv_nbsteps_varday = 1.0 / float(ctx.time.nbsteps_varday);
    }
    else
    {
        cout << "ERROR with the daily variation file" << endl;
    }

    ctx.climate.WDailyMean_all.assign(ctx.time.nbsteps_varday, float());
    ctx.climate.tDailyMean_all.assign(ctx.time.nbsteps_varday, float());
    ctx.climate.VPDDailyMean_all.assign(ctx.time.nbsteps_varday, float());
    ctx.climate.windDailyMean_all.assign(ctx.time.nbsteps_varday, float());

    for (int j = 0; j < ctx.time.nbsteps_varday; j++)
    {
        ctx.climate.WDailyMean_all[j] = ctx.climate.tDailyMean_all[j] = ctx.climate.VPDDailyMean_all[j] = ctx.climate.windDailyMean_all[j] = 0.0;
    }

    ctx.climate.WDailyMean_year = 0.0;

    for (int i = 0; i < ctx.time.nbdays; i++)
    {

        float DailyMeanTemperature_current = 0.0, DailyMeanWindSpeed_current = 0.0, DailyMeanIrradiance_current = 0.0, DailyMeanVapourPressureDeficit_current = 0.0;

        for (int j = 0; j < ctx.time.nbsteps_varday; j++)
        {
            DailyMeanTemperature_current += ctx.climate.varday_T[i * ctx.time.nbsteps_varday + j];
            DailyMeanWindSpeed_current += ctx.climate.varday_WS[i * ctx.time.nbsteps_varday + j];
            DailyMeanIrradiance_current += ctx.climate.varday_light[i * ctx.time.nbsteps_varday + j];
            DailyMeanVapourPressureDeficit_current += ctx.climate.varday_vpd[i * ctx.time.nbsteps_varday + j];

            ctx.climate.WDailyMean_all[j] += ctx.climate.varday_light[i * ctx.time.nbsteps_varday + j];
            ctx.climate.VPDDailyMean_all[j] += ctx.climate.varday_vpd[i * ctx.time.nbsteps_varday + j];
            ctx.climate.tDailyMean_all[j] += ctx.climate.varday_T[i * ctx.time.nbsteps_varday + j];
            ctx.climate.windDailyMean_all[j] += ctx.climate.varday_WS[i * ctx.time.nbsteps_varday + j];
        }

        DailyMeanTemperature_current *= ctx.time.inv_nbsteps_varday;
        DailyMeanWindSpeed_current *= ctx.time.inv_nbsteps_varday;
        DailyMeanIrradiance_current *= ctx.time.inv_nbsteps_varday;
        DailyMeanVapourPressureDeficit_current *= ctx.time.inv_nbsteps_varday;

        ctx.climate.WDailyMean_year += DailyMeanIrradiance_current;

        ctx.climate.DailyMeanTemperature.push_back(DailyMeanTemperature_current);
        ctx.climate.DailyMeanWindSpeed.push_back(DailyMeanWindSpeed_current);
        ctx.climate.DailyMeanIrradiance.push_back(DailyMeanIrradiance_current);
        ctx.climate.DailyMeanVapourPressureDeficit.push_back(DailyMeanVapourPressureDeficit_current);
    }

    for (int j = 0; j < ctx.time.nbsteps_varday; j++)
    {
        ctx.climate.WDailyMean_all[j] /= ctx.time.nbdays;
        ctx.climate.VPDDailyMean_all[j] /= ctx.time.nbdays;
        ctx.climate.tDailyMean_all[j] /= ctx.time.nbdays;
        ctx.climate.windDailyMean_all[j] /= ctx.time.nbdays;
    }

    ctx.climate.WDailyMean_year *= ctx.params.SWtoPPFD / ctx.time.nbdays;

    ctx.climate.tnight = ctx.climate.NightTemperature[0];
    ctx.climate.precip = ctx.climate.Rainfall[0];
    ctx.climate.WSDailyMean = ctx.climate.DailyMeanWindSpeed[0];
    ctx.climate.WDailyMean = ctx.climate.DailyMeanIrradiance[0] * ctx.params.SWtoPPFD;
    ctx.climate.tDailyMean = ctx.climate.DailyMeanTemperature[0];
    ctx.climate.VPDDailyMean = ctx.climate.DailyMeanVapourPressureDeficit[0];
}

//! Global function: This function reads inputs from the environmental variation input file
//! - v.3.0 of the code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=ctx.time.iter and not ctx.time.iterperyear) and change ctx.time.iterperyear by ctx.time.nbiter here.
//! @param ctx simulation context (all model state)
void ReadInputClimate(Context &ctx)
{
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_climate << endl;

    fstream InClim(ctx.fileio.inputfile_climate, ios::in);

    if (InClim)
    {
        // we go through all lines in the input file
        ctx.time.nbdays = 0;
        InClim.getline(ctx.fileio.buffer, 256, '\n');
        // cout << "Header line: " << ctx.fileio.buffer << endl;

        string line;
        while (getline(InClim, line))
        {
            istringstream linestream(line);

            float NightTemperature_current, Rainfall_current;
            linestream >> NightTemperature_current >> Rainfall_current;

            ctx.climate.NightTemperature.push_back(NightTemperature_current);
            ctx.climate.Rainfall.push_back(Rainfall_current);
            ctx.time.nbdays++;
        }

#ifdef WATER
        ctx.time.iterperyear = 365;
#else
        // if (ctx.time.nbdays%365==0) ctx.time.iterperyear=365; //check if this guarantees ctx.time.iterperyear is correctly infered
        ctx.time.iterperyear = 12;
#endif

        ctx.time.timestep = 1.0 / float(ctx.time.iterperyear);
        cout << "Read in climate data for " << ctx.time.nbdays << " days, with" << ctx.time.iterperyear << " iterations per year." << endl;
        // choose average conditions

        ctx.climate.Tnight_year = 0.0;

        for (int i = 0; i < ctx.time.nbdays; i++)
        {
            ctx.climate.Tnight_year += ctx.climate.NightTemperature[i];
        }

        ctx.climate.Tnight_year *= 1.0 / float(ctx.time.nbdays);

        cout << "Successfully read the climate file" << endl;
    }
    else
    {
        cout << "ERROR with the climate file" << endl;
    }
    InClim.close();
}

#else // FULL_CLIMATE

//! Global function: This function reads inputs from the environmental daily variation input file
//! @param ctx simulation context (all model state)
void ReadInputDailyvar(Context &ctx)
{
    // currently very simple reading in, only basic error checking
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_daytimevar << endl;

    fstream InDaily(ctx.fileio.inputfile_daytimevar, ios::in);

    if (InDaily)
    {
        InDaily.getline(ctx.fileio.buffer, 256, '\n'); // read in header

        // we go through all lines in the input file
        ctx.time.nbsteps_varday = 0;
        ctx.time.nbhours_covered = 0.0;

        string line;
        while (getline(InDaily, line))
        {
            istringstream linestream(line);

            float starttime_current, endtime_current, varday_light_current, varday_vpd_current, varday_T_current, varday_WS_current;
            linestream >> starttime_current >> endtime_current >> varday_light_current >> varday_vpd_current >> varday_T_current >> varday_WS_current;
            float nbhours_current;
            if (endtime_current > starttime_current)
                nbhours_current = endtime_current - starttime_current; // just calculate the time window
            else
                nbhours_current = 24.0 - (starttime_current - endtime_current); // i.e. if starttime is 23.5 and endtime is 0.5, etc. Calculate the inverse time window and then subtract from 24

            ctx.time.nbhours_covered += nbhours_current;
            ctx.climate.varday_light.push_back(varday_light_current);
            ctx.climate.varday_vpd.push_back(varday_vpd_current);
            ctx.climate.varday_T.push_back(varday_T_current);
            ctx.climate.varday_WS.push_back(varday_WS_current);
            cout << "at: " << starttime_current << " vardaytime_light: " << ctx.climate.varday_light[ctx.time.nbsteps_varday] << " vardaytime_vpd: " << ctx.climate.varday_vpd[ctx.time.nbsteps_varday] << " vardaytime_T: " << ctx.climate.varday_T[ctx.time.nbsteps_varday] << " vardaytime_WS: " << ctx.climate.varday_WS[ctx.time.nbsteps_varday] << endl;
            ctx.time.nbsteps_varday++;
        }

#ifdef WATER
        totday_light = 0.0;
        for (int i = 0; i < ctx.time.nbsteps_varday; i++)
        {
            totday_light += ctx.climate.varday_light[i];
            cout << i << "\t" << ctx.climate.varday_light[i] << endl;
        }
        cout << totday_light << endl;
        totday_light *= 1.0 / float(ctx.time.nbsteps_varday);
        cout << totday_light << endl;
#endif
        cout << "Read in: " << ctx.time.nbsteps_varday << " timesteps per day, covering " << ctx.time.nbhours_covered << " hours of the day." << endl;
        cout << "Successfully read in daytime variation file" << endl;
        ctx.time.inv_nbsteps_varday = 1.0 / float(ctx.time.nbsteps_varday);
    }
    else
    {
        cout << "ERROR with the daily variation file" << endl;
    }
}

//! Global function: This function reads inputs from the environmental variation input file
//! - v.3.0 of the code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=ctx.time.iter and not ctx.time.iterperyear) and change ctx.time.iterperyear by ctx.time.nbiter here.
//! @param ctx simulation context (all model state)
void ReadInputClimate(Context &ctx)
{
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_climate << endl;

    fstream InClim(ctx.fileio.inputfile_climate, ios::in);

    if (InClim)
    {
        // we go through all lines in the input file
        ctx.time.iterperyear = 0;
        InClim.getline(ctx.fileio.buffer, 256, '\n');
        // cout << "Header line: " << ctx.fileio.buffer << endl;

        string line;
        while (getline(InClim, line))
        {
            istringstream linestream(line);

            float DailyMeanTemperature_current, NightTemperature_current, Rainfall_current, DailyMeanWindSpeed_current, DailyMeanIrradiance_current, DailyMeanVapourPressureDeficit_current;
            linestream >> DailyMeanTemperature_current >> NightTemperature_current >> Rainfall_current >> DailyMeanWindSpeed_current >> DailyMeanIrradiance_current >> DailyMeanVapourPressureDeficit_current;

            ctx.climate.DailyMeanTemperature.push_back(DailyMeanTemperature_current);
            ctx.climate.NightTemperature.push_back(NightTemperature_current);
            ctx.climate.Rainfall.push_back(Rainfall_current);
            ctx.climate.DailyMeanWindSpeed.push_back(DailyMeanWindSpeed_current);
            ctx.climate.DailyMeanIrradiance.push_back(DailyMeanIrradiance_current);
            ctx.climate.DailyMeanVapourPressureDeficit.push_back(DailyMeanVapourPressureDeficit_current);
            // cout <<  DailyMeanTemperature_current << "\t" <<  NightTemperature_current << "\t" <<  Rainfall_current << "\t" <<  DailyMeanWindSpeed_current << "\t" <<  DailyMeanIrradiance_current << "\t" <<   DailyMeanVapourPressureDeficit_current << endl;
            ctx.time.iterperyear++;
        }

        ctx.time.timestep = 1.0 / float(ctx.time.iterperyear);
        cout << "Read in climate data for " << ctx.time.iterperyear << " iterations per year." << endl;
        // choose average conditions
        ctx.climate.WDailyMean_year = ctx.climate.tDailyMean_year = ctx.climate.VPDDailyMean_year = ctx.climate.Tnight_year = 0.0;
#ifdef WATER
        ctx.climate.windDailyMean_year = 0.0;
#endif

        for (int i = 0; i < ctx.time.iterperyear; i++)
        {
            // ctx.climate.DailyMeanVapourPressureDeficit[i]*=1.33; test main_1.33VPD on MESO@LR - 13/05/22

            ctx.climate.WDailyMean_year += ctx.climate.DailyMeanIrradiance[i] * ctx.params.SWtoPPFD;
            ctx.climate.tDailyMean_year += ctx.climate.DailyMeanTemperature[i];
            ctx.climate.VPDDailyMean_year += ctx.climate.DailyMeanVapourPressureDeficit[i];
#ifdef WATER
            ctx.climate.windDailyMean_year += ctx.climate.DailyMeanWindSpeed[i];
#endif
            ctx.climate.Tnight_year += ctx.climate.NightTemperature[i];
        }

        ctx.climate.WDailyMean_year *= 1.0 / float(ctx.time.iterperyear);
        ctx.climate.tDailyMean_year *= 1.0 / float(ctx.time.iterperyear);
        ctx.climate.VPDDailyMean_year *= 1.0 / float(ctx.time.iterperyear);
#ifdef WATER
        ctx.climate.windDailyMean_year *= 1.0 / float(ctx.time.iterperyear);
#endif
        // cout << "ctx.climate.WDailyMean: " << ctx.climate.WDailyMean_year << endl;
        ctx.climate.Tnight_year *= 1.0 / float(ctx.time.iterperyear);

        ctx.climate.tnight = ctx.climate.NightTemperature[ctx.time.iter % ctx.time.iterperyear];
        ctx.climate.precip = ctx.climate.Rainfall[ctx.time.iter % ctx.time.iterperyear];
        ctx.climate.WSDailyMean = ctx.climate.DailyMeanWindSpeed[ctx.time.iter % ctx.time.iterperyear];
        ctx.climate.WDailyMean = ctx.climate.DailyMeanIrradiance[ctx.time.iter % ctx.time.iterperyear] * ctx.params.SWtoPPFD;
        ctx.climate.tDailyMean = ctx.climate.DailyMeanTemperature[ctx.time.iter % ctx.time.iterperyear];
        ctx.climate.VPDDailyMean = ctx.climate.DailyMeanVapourPressureDeficit[ctx.time.iter % ctx.time.iterperyear];

        cout << "Successfully read the climate file" << endl;
    }
    else
    {
        cout << "ERROR with the climate file" << endl;
    }
    InClim.close();
}

#endif // FULL_CLIMATE

#ifdef WATER
//! Global function: This function reads inputs from the soil input file
//! - in v.3.0, all soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each ctx.diag.layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
//! - ctx.soil.Sat_SWC and ctx.soil.Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
//! @param ctx simulation context (all model state)
void ReadInputSoil(Context &ctx)
{
    cout << endl
         << "Reading in file: " << ctx.fileio.inputfile_soil << endl;
    fstream InSoil(ctx.fileio.inputfile_soil, ios::in);

    if (InSoil)
    {
        InSoil.getline(ctx.fileio.buffer, 256, '\n');
        vector<float> layer_thickness, proportion_Silt, proportion_Clay, proportion_Sand; // in ctx.params.m, %, %,%
        vector<float> SOC, DBD, pH, CEC;                                                  // soil organic content, provided in %; dry bulk density, in g cm-3; pH; cation exchange capacity, in cmol kg-1
        SOC.reserve(20);
        DBD.reserve(20);
        pH.reserve(20);
        CEC.reserve(20);

        layer_thickness.reserve(20);
        proportion_Silt.reserve(20);
        proportion_Clay.reserve(20);
        proportion_Sand.reserve(20);
        ctx.soil.nblayers_soil = 0;

        // we go through all lines in the input file
        string line;
        while (getline(InSoil, line))
        {
            istringstream linestream(line);

            if (ctx.opt._WATER_RETENTION_CURVE == 1)
            {
                float thickness_current, proportion_Silt_current, proportion_Clay_current, proportion_Sand_current, SOC_current, DBD_current, pH_current, CEC_current;
                linestream >> thickness_current >> proportion_Silt_current >> proportion_Clay_current >> proportion_Sand_current >> SOC_current >> DBD_current >> pH_current >> CEC_current;
                layer_thickness.push_back(thickness_current);
                proportion_Silt.push_back(proportion_Silt_current);
                proportion_Clay.push_back(proportion_Clay_current);
                proportion_Sand.push_back(proportion_Sand_current);
                SOC_current *= 10; // to convert SOC in % to gC kg-1 soil, to apply Hodnett & Tomasella 2002
                SOC.push_back(SOC_current);
                DBD.push_back(DBD_current);
                pH.push_back(pH_current);
                CEC.push_back(CEC_current);
            }
            else if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {
                float thickness_current, proportion_Silt_current, proportion_Clay_current, proportion_Sand_current;
                linestream >> thickness_current >> proportion_Silt_current >> proportion_Clay_current >> proportion_Sand_current;
                layer_thickness.push_back(thickness_current);
                proportion_Silt.push_back(proportion_Silt_current);
                proportion_Clay.push_back(proportion_Clay_current);
                proportion_Sand.push_back(proportion_Sand_current);
            }

            ctx.soil.nblayers_soil++;
        }
        cout << "Read in: " << ctx.soil.nblayers_soil << " soil layers" << endl;

        ctx.soil.layer_depth.assign(ctx.soil.nblayers_soil, float());
        float cumulative_depth = 0.0;
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            cumulative_depth += layer_thickness[l];
            ctx.soil.layer_depth[l] = cumulative_depth;
        }
        // (added in the header but kept in here) in this version, all soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each ctx.diag.layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
        // (added in the header but kept in here) ctx.soil.Sat_SWC and ctx.soil.Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
        ctx.soil.Sat_SWC.assign(ctx.soil.nblayers_soil, float());
        ctx.soil.Max_SWC.assign(ctx.soil.nblayers_soil, float());
        ctx.soil.FC_SWC.assign(ctx.soil.nblayers_soil, float());

        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {

            if (ctx.opt._WATER_RETENTION_CURVE == 1)
            {
                ctx.soil.Sat_SWC[l] = 0.01 * (81.799 + (0.099 * proportion_Clay[l]) - (31.42 * DBD[l]) + (0.018 * CEC[l]) + (0.451 * pH[l]) - (0.0005 * proportion_Sand[l] * proportion_Clay[l])); // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.ctx.params.m-3
            }
            else if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {
                ctx.soil.Sat_SWC[l] = 0.01 * (40.61 + (0.165 * proportion_Silt[l]) + (0.162 * proportion_Clay[l]) + (0.00137 * proportion_Silt[l] * proportion_Silt[l]) + (0.000018 * proportion_Silt[l] * proportion_Silt[l] * proportion_Clay[l])); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.ctx.params.m-3
            }
            ctx.soil.Max_SWC[l] = ctx.soil.Sat_SWC[l] * ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l]; // in m3
            cout << "ctx.diag.layer " << l << " Vol=" << ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l] << " m3; ctx.soil.Sat_SWC =" << ctx.soil.Sat_SWC[l] << " MAX_SWC =" << ctx.soil.Max_SWC[l] << " m3." << endl;
        }

        ctx.soil.Ksat.assign(ctx.soil.nblayers_soil, float());
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            ctx.soil.Ksat[l] = 0.007055556 * pow(10, (-0.60 - (0.0064 * proportion_Clay[l]) + (0.0126 * proportion_Sand[l]))); // according to Cosby et al. 1984 (the only expression of k_sat reported in Table 2 of Marthews et al. 2014). k_sat is here in mm/s or equivalently in kg/m2/s.
            cout << "ctx.diag.layer " << l << " ctx.soil.Ksat=" << ctx.soil.Ksat[l] << "mm/s or kg/m2/s  " << ctx.soil.Ksat[l] * 9.8 / 18 << endl;
        }

        ctx.soil.Res_SWC.assign(ctx.soil.nblayers_soil, float());
        ctx.soil.Min_SWC.assign(ctx.soil.nblayers_soil, float());

        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {

            if (ctx.opt._WATER_RETENTION_CURVE == 1)
            {
                ctx.soil.Res_SWC[l] = 0.01 * (22.733 - (0.164 * proportion_Sand[l]) + (0.235 * CEC[l]) - (0.831 * pH[l]) + (0.0018 * proportion_Clay[l] * proportion_Clay[l]) + (0.0026 * proportion_Sand[l] * proportion_Clay[l])); // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.ctx.params.m-3
            }
            if (ctx.opt._WATER_RETENTION_CURVE == 0)
            {
                ctx.soil.Res_SWC[l] = 0.01 * fmaxf(0.0, (-2.094 + (0.047 * proportion_Silt[l]) + (0.431 * proportion_Clay[l]) - (0.00827 * proportion_Silt[l] * proportion_Clay[l]))); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
            }
            ctx.soil.Min_SWC[l] = ctx.soil.Res_SWC[l] * ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l]; // in m3

            cout << "ctx.diag.layer " << l << " Vol=" << ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l] << "m3; Res=" << ctx.soil.Res_SWC[l] << " MIN_SWC =" << ctx.soil.Min_SWC[l] << " m3" << endl;
        }

        if (ctx.opt._WATER_RETENTION_CURVE == 1)
        {
            ctx.soil.a_vgm.assign(ctx.soil.nblayers_soil, float());
            ctx.soil.b_vgm.assign(ctx.soil.nblayers_soil, float());
            ctx.soil.c_vgm.assign(ctx.soil.nblayers_soil, float());
            ctx.soil.m_vgm.assign(ctx.soil.nblayers_soil, float());

            for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            {

                float alpha_vgm = 1000 * exp((-2.294 - (3.526 * proportion_Silt[l]) + (2.440 * (0.1 * SOC[l])) - (0.076 * CEC[l]) - (11.331 * pH[l]) + (0.019 * proportion_Silt[l] * proportion_Silt[l])) * 0.01); // RENAMED: ctx.params.alpha → alpha_vgm (van Genuchten-Mualem soil parameter; distinct from global ctx.params.alpha = apparent quantum yield), in MPa-1 (ie after already dividing by rho*g, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                ctx.soil.a_vgm[l] = -1.0 / alpha_vgm;
                float n_vgm = exp((62.986 - (0.833 * proportion_Clay[l]) - (0.529 * (SOC[l] * 0.1)) + (0.593 * pH[l]) + (0.007 * proportion_Clay[l] * proportion_Clay[l]) - (0.014 * proportion_Sand[l] * proportion_Silt[l])) * 0.01); // this is the n parameter of the van Genuchten-Mualem model, dimensionless, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                ctx.soil.m_vgm[l] = 1.0 - 1.0 / n_vgm;
                ctx.soil.b_vgm[l] = 1.0 / ctx.soil.m_vgm[l];
                ctx.soil.c_vgm[l] = 1.0 - ctx.soil.m_vgm[l];

                ctx.soil.FC_SWC[l] = (ctx.soil.Res_SWC[l] + (ctx.soil.Sat_SWC[l] - ctx.soil.Res_SWC[l]) * pow((pow(0.01 * alpha_vgm, 1 / ctx.soil.c_vgm[l]) + 1), -(1 / ctx.soil.b_vgm[l]))) * ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l]; // this is the ctx.diag.layer water content at field capacity, in m3. As in Marthews et al. 2014 (cf. note in Table 2), we used Phi at FC=-10kPa and not -33kPa, following Marshall et al., 1996; Townend et al., 2001; Tomasella and Hodnett, 2004)

                cout << "ctx.diag.layer " << l << " alpha_vgm=" << alpha_vgm << "\t" << "n_vgm=" << n_vgm << " ctx.soil.FC_SWC=" << ctx.soil.FC_SWC[l] << endl;
            }
        }
        else if (ctx.opt._WATER_RETENTION_CURVE == 0)
        {
            ctx.soil.phi_e.assign(ctx.soil.nblayers_soil, float());
            ctx.soil.b.assign(ctx.soil.nblayers_soil, float());
            for (int l = 0; l < ctx.soil.nblayers_soil; l++)
            {
                ctx.soil.phi_e[l] = -0.001 * (0.285 + (0.000733 * proportion_Silt[l] * proportion_Silt[l]) - (0.00013 * proportion_Silt[l] * proportion_Clay[l]) + (0.0000036 * proportion_Silt[l] * proportion_Silt[l] * proportion_Clay[l])); // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. ctx.soil.phi_e is here in MPa.
                ctx.soil.b[l] = exp(1.197 + (0.00417 * proportion_Silt[l]) - (0.0045 * proportion_Clay[l]) + (0.000894 * proportion_Silt[l] * proportion_Clay[l]) - (0.00001 * proportion_Silt[l] * proportion_Silt[l] * proportion_Clay[l]));  // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
                // ctx.soil.phi_e[l]=-0.00000001*pow(10.0,(2.17-(0.0063*proportion_Clay[l])-(0.0158*proportion_Sand[l])))*(1000*9.80665); // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. In MPa.
                // ctx.soil.b[l]=3.10+0.157*proportion_Clay[l]-0.003*proportion_Sand[l]; // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. Dimensionless.

                ctx.soil.FC_SWC[l] = (ctx.soil.Res_SWC[l] + (ctx.soil.Sat_SWC[l] - ctx.soil.Res_SWC[l]) * pow(-0.01 / ctx.soil.phi_e[l], -(1 / ctx.soil.b[l]))) * ctx.grid.sites_per_dcell * ctx.grid.LH * ctx.grid.LH * layer_thickness[l]; // this is the ctx.diag.layer water content at filed capacity, in m3. As in Marthews et al. 2014 (cf. note in Table 2), we used Phi at FC=-10kPa and not -33kPa, following Marshall et al., 1996; Townend et al., 2001; Tomasella and Hodnett, 2004)

                cout << "ctx.diag.layer " << l << " ctx.soil.phi_e=" << ctx.soil.phi_e[l] << "\t" << "b=" << ctx.soil.b[l] << " ctx.soil.FC_SWC=" << ctx.soil.FC_SWC[l] << endl;
            }
        }
        cout << "Successfully read the soil file" << endl;
    }
    else
    {
        cout << "ERROR with the soil file" << endl;
    }
    cout << endl;
}
#endif
