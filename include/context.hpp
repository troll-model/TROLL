#ifndef TROLL_CONTEXT_HPP
#define TROLL_CONTEXT_HPP

#include <fstream>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "species.hpp"
#include "tree.hpp"

//! Input and output file names
struct FileIO
{
    char buffer[256]; //!< Character strings used to read file names, or other features
    char inputfile[256]; //!< Character strings used to read file names, or other features
    char inputfile_daytimevar[256]; //!< Character strings used to read file names, or other features
    char inputfile_climate[256]; //!< Character strings used to read file names, or other features
    char inputfile_soil[256]; //!< Character strings used to read file names, or other features
    char outputinfo[256]; //!< Character strings used to read file names, or other features
    char inputfile_inventory[256]; //!< Character strings used to read file names, or other features
    char inputfile_pointcloud[256]; //!< Character strings used to read file names, or other features
    char inputfile_SWC[256]; //!< Path of the soil water content file (option -w), e.g. the _final_SWC3D.txt output of a previous run
    char inputfile_species[256]; //!< Vector of input files
};

//! Raw command-line arguments, run index (-n) and random seed (-r)
struct InputBuffers
{
    char *bufi = nullptr; //!< Command-line value of -i: path of the global parameter file (points into argv)
    char *bufi_daytimevar = nullptr; //!< Command-line value of -d: path of the half-hourly climate file (points into argv)
    char *bufi_climate = nullptr; //!< Command-line value of -m: path of the daily climate file (points into argv)
    char *bufi_soil = nullptr; //!< Command-line value of -p: path of the soil file (points into argv)
    char *buf = nullptr; //!< Command-line value of -o: prefix of all output file names (points into argv)
    char *bufi_data = nullptr; //!< Command-line value of -f: path of the forest inventory used as initial state (points into argv)
    char *bufi_pointcloud = nullptr; //!< Command-line value of -l: path of the point-cloud (lidar simulation) parameter file (points into argv)
    char *bufi_dataSWC = nullptr; //!< Command-line value of -w: path of the initial soil water content file (points into argv)
    char *bufi_species = nullptr; //!< Command-line value of -s: path of the species parameter file (points into argv)
    int run_index = 0;       //!< Run index from the -n option, used in output file names (formerly the global easympi_rank)
    long long seed = -1;     //!< RNG seed from the -r / --seed option; -1 = not set on the command line
};

//! Output file streams
struct OutputConfig
{
    std::fstream output_info; //!< Basic simulation information
    std::fstream output_basic[4]; //!< Default output streams, always used
    std::fstream output_extended[9]; //!< Extended TROLL outputs, preserved from previous versions, might need further clean-up
    std::fstream output_visual[2]; //!< Outputs for visualization/gif creation, new in v.3.1.2
    std::fstream output_pointcloud; //!< Optional las file output, new in v.3.1.6
#ifdef MIP_Lichstein
    std::fstream output_MIP_eco; //!< Ecosystem-level output of the Lichstein model intercomparison (MIP_Lichstein option)
    std::fstream output_MIP_ind; //!< Individual-level output of the Lichstein model intercomparison (MIP_Lichstein option)
#endif
#ifdef Output_ABC
    std::fstream output_abc[11]; //!< Output streams for Approximate Bayesian Computation
#endif
#ifdef WATER
    std::fstream output[40]; //!< Output files
#endif
#ifdef TRACK_INDIVIDUALS
    std::fstream output_track[3]; //!< Output streams for tracking of trees
#endif
};

//! Model options (flags) read from the global parameter file
struct ModelOptions
{
    bool _NONRANDOM; //!< User control: If _NONRANDOM == 1, the seeds for the random number generators will be kept fixed at 1, for bug fixing
    bool _GPPcrown; //!< User control: This defines an option to compute only GPP from the topmost value of PPFD and GPP, instead of looping within the crown.
    bool _BASICTREEFALL; //!< User control: if defined: treefall is a source of tree death (and if TREEFALL not defined, this is modeled through simple comparison between tree height and a threshold t_Ct, if not defined, treefall is not represented as a separated and independent source of death, but instead, all tree death are due to the deathrate value)
    bool _SEEDTRADEOFF; //!< User control: if defined: the number of seeds produced by each tree is determined by the tree NPP allocated to reproduction and the species seed mass, otherwise the number of seeds is fixed; besides, seedling recruitment in one site is not made by randomly and 'equiprobably' picking one species among the seeds present at that site but the probability of recruitment among the present seeds is proportional to the number of seeds (in s_Seed[site]) time the seed mass of each species
    bool _NDD; //!< User control: if defined, negative density dependant processes affect both the probability of seedling recruitment and the local tree death rate. The term of density-dependance is computed as the sum of conspecific tree basal area divided by their distance to the focal tree within a neighbourhood (circle of radius 15m)
    bool _CROWN_MM; //!< User control: new in v.2.4.1: Michaelis Menten allometry for crowns instead of power law, since v.2.5: not a macro anymore, but set at runtime (little overhead since queried only once at birth and once per timestep), !!!: power law is the default. If Michaelis Menten type allometry is activated, the parameters have to be changed in input sheet accordingly
    bool _sapwood; //!< User control: two ways of parameterising sapwood density: constant thickness (0), Fyllas, but with lower limit (1)
    bool _seedsadditional; //!< User control: excess carbon into seeds? no/yes=(0/1)
    bool _LL_parameterization; //!< User control: two ways for parameterising leaf lifespan: empirical (derived by Sylvain Schmitt, TODO: from which data?), Kikuzawa model (0,1)
    bool _FromInventory; //!< User control: if defined, an additional input file can be provided to start simulations from an existing data set or a simulated data set (5 parameters are needed: x and y coordinates, dbh, species_label, species
    bool _OUTPUT_extended; //!< User control: changed in v.3.1 from _OUTPUT_reduced to _OUTPUT_extended, uses extended set of output streams
    bool _OUTPUT_inventory; //!< User control: added in v4.0 to ease stop and start from inventories.
    int _LA_regulation; //!< User control: updated v.3.1: potentially three ways of parameterising leaf dynamic allocation, but currently using only two ways: no regulation (0), never exceed LAImax, i.e. the maximum LAI under full sunlight (1), adjust LAI to the current light environment (2). To switch between option 1 and 2, only one line is necessary in CalcLAmax()
    int _OUTPUT_pointcloud; //!< User control: ATTENTION! At the moment assumes a little-endian system (most personal computers, but not necessarily server systems), because LAS fles are in little-endian! If == 1, creates a point cloud from a simplified ALS simulation;
    int _SOIL_LAYER_WEIGHT; //!< !< User control: three ways of computing the fraction of transpiration supplied by each soil layer in individual tree total transpiration, and to weight the tree average water potential in the root zone. (0: root biomass only, 1: relative root-to soil conductance, 2: relative estimated maximal transpiration as in Duursma & Medlyn 2012).
    int _WATER_RETENTION_CURVE; //!< !< User control: different water retention curves can be used. So far two are have been implemented: either brooks & Corey (0), either van Genuchten-Mualem (1). To each wtare retention cirve option is associated a different set of pedo-transfer functions, see Table 2 (texture-based Tomasella & Hodnett 1998
};

//! Spatial grid: plot dimensions, vertical layers, soil cells (dcells), borders
struct Grid
{
    int sites; //!< Number of pixels in the scene (cols*rows)
    int cols; //!< Number of columns in the scene
    int rows; //!< Number of rows in the scene
    int nbspp; //!< Number of species
    int length_dcell; //!< Linear size of a dcell (coarser grid of the simulated scene); default since v.3.0 standard
    int linear_nb_dcells; //!< Linear number of dcells
    int sites_per_dcell; //!< Number of sites per dcell
    int nbdcells; //!< Total number of dcells; note that nbdcells = linear_nb_dcells * linear_nb_dcells
    std::vector<int> site_DCELL; //!< Index of the soil cell (dcell) containing each site; size sites
    float i_sites_per_dcell; //!< 1.0/sites_per_dcell, where sites_per_dcell is number of sites per dcell
    int HEIGHT; //!< Maximum height (m)
    int SBORD; //!< Boundary condition, defined as RMAX*cols
    float NV; //!< Number of cells per m (vertical)
    float NH; //!< Number of cells per m (horizontal)
    float LV; //!< LV = 1.0/NV; NV is vertical number of cells per m
    float LH; //!< LH = 1.0/NH; NH is horizontal number of cells per m
    int RMAX; //!< Maximum crown radius (m)
    int dbhmaxincm; //!< Maximum diameter at breast height (DBH) in cm (the convention in forestry; !!! departure from usual SI units)
    int leafdem_resolution; //!< Resolution for leaf demography model !!!MORE_COMMENT
};

//! Time stepping state: current iteration, timestep length, number of iterations and days
struct TimeState
{
    int iterperyear; //!< Number of iterations per year (=12 if monthly timestep, =365 if daily timestep)
    int nbiter; //!< Total number of timesteps
    int iter; //!< Current timestep
    int nbout; //!< Number of outputs
    int freqout; //!< Frequency HDF outputs
    int nbdays; //!< Number of days with explicit daily climate variation provided in input
    int nbsteps_varday; //!< Number of steps for environmental variables
    float inv_nbsteps_varday; //!< The inverse of the step number for environmental variables
    float nbhours_covered; //!< Number of hours of environmental variables in the environmental input file
    float timestep; //!< Duration of one timestep (in years)=1/iterperyear
};

//! Climate forcing: daily series, half-hourly series and current-step values
struct Climate
{
    std::vector<float> DailyMeanTemperature; //!< Daily mean temperature in degree Celsius
    std::vector<float> DailyMeanIrradiance; //!< Daily mean irradiance (W/m^2)
    std::vector<float> DailyMeanVapourPressureDeficit; //!< Daily mean vapour pressure deficit (kPa)
    std::vector<float> NightTemperature; //!< Night temperature in degree Celsius
    std::vector<float> Rainfall; //!< Rainfall (mm)
    std::vector<float> DailyMeanWindSpeed; //!< Wind speed (m/s)

    std::vector<float> varday_light; //!< Light irradiance variation during an average day, since v.2.5: standardized with respect to the mean and summing to 0 (adimensional)
    std::vector<float> varday_vpd; //!< Vapour pressure deficit (VPD) variation during an average day, since v.2.5:  standardized with respect to the mean and summing to 0 (adimensional)
    std::vector<float> varday_T; //!< Temperature variation during an average day, since v.2.5:  standardized with respect to the mean and summing to 0 (adimensional)
    std::vector<float> varday_WS; //!< Wind speed variation during an average day, since v.2.5:  standardized with respect to the mean and summing to 0 (adimensional)

    float tnight; //!< Night mean temperature (degree C)
    float precip; //!< Rainfall  (mm)
    float WSDailyMean; //!< WindSpeed (m/s)
    float WDailyMean; //!< Daily mean irradiance (average for timestep) (micromol PAR photon/m^2/s), used in the photosynthesis part. !!! W/m2 is the common unit of meteorological stations ==> need to be converted using SWtoPPFD
    float tDailyMean; //!< Daily mean temperature (degree C)
    float VPDDailyMean; //!< Daily mean VapourPressureDeficit (kPa)
    float WDailyMean_year; //!< Average WDailyMean per year
    float tDailyMean_year; //!< Average tDailyMean per year
    float VPDDailyMean_year; //!< Average VPDDailyMean per year
    float windDailyMean_year; //!< Average WindSpeed per year
    float Tnight_year; //!< Average tnight per year

    std::vector<float> WDailyMean_all; //!< Mean daily course of irradiance over the whole climate file, per daytime step; used when iter = -1 (LAImax table)
    std::vector<float> VPDDailyMean_all; //!< Mean daily course of VPD over the whole climate file, per daytime step; used when iter = -1 (LAImax table)
    std::vector<float> tDailyMean_all; //!< Mean daily course of air temperature over the whole climate file, per daytime step; used when iter = -1 (LAImax table)
    std::vector<float> windDailyMean_all; //!< Mean daily course of wind speed over the whole climate file, per daytime step; used when iter = -1 (LAImax table)
};

//! Precomputed lookup tables: temperature responses (0.1 degC bins), light and longwave extinction, wind profile, crown sites, LAImax
struct LookupTables
{
    int nbTbins; //!< Number of bins for the temperature lookup tables
    float iTaccuracy; //!< Inverse of accuracy of a temperature bin (e.g. if Taccuracy is 0.1 or 0.5 °C, then iTaccuracy is 10.0 or 2.0, respectively)

    std::vector<float> LookUp_KmT; //!< Michaelis-Menten constant of Rubisco for CO2, including O2 inhibition, per 0.1 degC temperature bin (von Caemmerer 2000); expressed relative to Cair in one of the two Photosyn variants
    std::vector<float> LookUp_GammaT; //!< CO2 compensation point in the absence of dark respiration (Gamma*), per 0.1 degC bin (von Caemmerer 2000); relative to Cair in one of the two Photosyn variants
    std::vector<float> LookUp_VcmaxT; //!< Temperature factor of Vcmax relative to 25 degC, per 0.1 degC bin (Bernacchi et al. 2003; dimensionless)
    std::vector<float> LookUp_JmaxT; //!< Temperature factor of Jmax relative to 25 degC, per 0.1 degC bin (Bernacchi et al. 2003; dimensionless)
    std::vector<float> LookUp_Rleaf; //!< Temperature factor of leaf dark respiration relative to 25 degC, per 0.1 degC bin (Atkin et al. 2015, eq. 1; dimensionless)
    std::vector<float> LookUp_flux_absorption; //!< Mean fraction of top-of-canopy PPFD within a leaf layer, by cumulated LAI above (i) and layer LAI (j), index i + 400 j
    std::vector<float> LookUp_flux; //!< Fraction of top-of-canopy PPFD reaching a given cumulated LAI, exp(-kpar * LAI); index i + 400 j
    std::vector<float> LookUp_ExtinctLW; //!< Extinction factor of thermal (longwave) radiation with cumulated LAI, 0.8 exp(-0.8 LAI) (Leuning et al. 1995, eq. D1)
    std::vector<float> LookUp_VPD; //!< Reduction factor of VPD with cumulated LAI within the canopy (dimensionless)
    std::vector<float> LookUp_T; //!< Decrease in air temperature with cumulated LAI within the canopy, up to 3 degC at LAI >= 7 (degC)
    std::vector<float> LookUp_Rstem; //!< Stem maintenance respiration coefficient per timestep at each 0.1 degC bin, Q10 = 2
    int LookUp_Crown_site[2601]; //!< new in v.2.4: lookup table to fill crown cylinder sequentially from inside to outside, allowing for smooth crown growth

    int nbVPDbins; //!< Number of bins for the VPD lookup tables
    float iVPDaccuracy; //!< Inverse of accuracy of a vpd bin (e.g. if VPDaccuracy is 0.01 kPa, then iVPDaccuracy is 100.0)

    std::vector<std::vector<float>> LookUp_INLR; //!< Isothermal net longwave radiation by temperature bin (0.1 degC) and VPD bin (0.01 kPa), with sky emissivity from Dilley & O'Brien 1998 (W m-2)
    std::vector<float> LookUp_SLOPE; //!< Slope of saturation vapour pressure versus temperature, per 0.1 degC bin (Pa K-1; Jones 2013, eq. 5.15)
    std::vector<float> LookUp_GRADN; //!< Radiative conductance 4 emissivity sigma T^3 / (CPAIR AIRMA), per 0.1 degC bin (mol m-2 s-1; Jones 2013, p. 101)

    int nbHbins; //!< Number of bins for the Wind lookup tables
    float iHaccuracy; //!< Inverse of the height-ratio step of LookUp_Wind (1 / 0.01)

    std::vector<float> LookUp_Wind; //!< Wind speed relative to canopy top, by height ratio h / canopy height in 0.01 steps: logarithmic profile above 0.8 H, exponential (Inoue) below

    std::vector<float> LookUpLAImax;      //!< LAImax per (species, intraspecific variant); filled lazily, see GetLAImaxPrecomputed()
    std::vector<unsigned char> LookUpLAImax_ready; //!< 1 once the matching LookUpLAImax entry has been computed
    int LookUpLAImax_iter = -1;           //!< ctx.time.iter at initialisation (the entries are defined for these conditions)
    float LookUpLAImax_tnight = 0.0f;     //!< ctx.climate.tnight at initialisation
};

//! Soil hydrology: layer properties, water content, potential and fluxes per layer and soil cell (dcell)
struct Soil
{
    int nblayers_soil; //!< Number of soil layers (for water module)
    std::vector<float> layer_depth; //!< Depth of the bottom of each soil layer below the surface (m)

    std::vector<float> Sat_SWC; //!< Saturated volumetric water content of each soil layer (m3 m-3)
    std::vector<float> Max_SWC; //!< Maximum water volume of each soil layer in one dcell (m3)
    std::vector<float> FC_SWC; //!< Water volume at field capacity of each soil layer in one dcell (m3)
    std::vector<float> Res_SWC; //!< Residual volumetric water content of each soil layer (m3 m-3)
    std::vector<float> Min_SWC; //!< Minimum (residual) water volume of each soil layer in one dcell (m3)

    std::vector<float> Ksat; //!< Saturated hydraulic conductivity of each soil layer (mm s-1, i.e. kg m-2 s-1; Cosby et al. 1984)
    std::vector<float> a_vgm; //!< van Genuchten-Mualem parameter -1/alpha of each soil layer (MPa, so that soil_phi3D is in MPa)
    std::vector<float> b_vgm; //!< van Genuchten-Mualem parameter 1/m of each soil layer (dimensionless)
    std::vector<float> c_vgm; //!< van Genuchten-Mualem parameter 1 - m of each soil layer (dimensionless)
    std::vector<float> m_vgm; //!< van Genuchten-Mualem parameter m = 1 - 1/n of each soil layer (dimensionless)
    std::vector<float> phi_e; //!< Air-entry water potential of each soil layer, Campbell model (MPa; Tomasella & Hodnett 1998)
    std::vector<float> b; //!< Exponent of the Campbell water retention curve of each soil layer (dimensionless)

    std::vector<std::vector<float>> SWC3D; //!< Soil water content per layer and dcell (m3)
    std::vector<std::vector<float>> soil_phi3D; //!< Soil water potential per layer and dcell (MPa)
    std::vector<std::vector<float>> Ks; //!< Soil hydraulic conductivity per layer and dcell at the current water content (mm s-1)
    std::vector<std::vector<float>> KsPhi; //!< Product of conductivity and water potential (Campbell model) per layer and dcell; used for the water potential in the root zone
    std::vector<std::vector<float>> Transpiration; //!< Water taken up by transpiration per layer and dcell during the current timestep (same unit as SWC3D)

    std::vector<std::vector<float>> LAI_DCELL; //!< Mean leaf area density per height layer and dcell (m2 m-2 per layer); [0] is the total LAI of the dcell

    std::vector<float> LAI_young; //!< Plot-mean LAI of young leaves per height layer (m2 m-2)
    std::vector<float> LAI_mature; //!< Plot-mean LAI of mature leaves per height layer (m2 m-2)
    std::vector<float> LAI_old; //!< Plot-mean LAI of old leaves per height layer (m2 m-2)

    std::vector<float> Canopy_height_DCELL; //!< Mean canopy height of each dcell (m)
    std::vector<int> HSum_DCELL; //!< Sum of the canopy heights of the sites of each dcell, used to compute Canopy_height_DCELL (m)
    std::vector<float> TopWindSpeed_DCELL; //!< Wind speed at the top of the canopy of each dcell (m s-1)

    std::vector<float> Interception; //!< Rainfall intercepted by the canopy of each dcell during the current timestep (mm of rainfall)
    std::vector<float> Throughfall; //!< Rainfall reaching the soil of each dcell during the current timestep (mm of rainfall)
    std::vector<float> Runoff; //!< Surface runoff of each dcell during the current timestep (mm of rainfall, as Throughfall)
    std::vector<float> Leakage; //!< Drainage below the deepest soil layer of each dcell during the current timestep (m3)
    std::vector<float> Evaporation; //!< Soil evaporation from the top layer of each dcell during the current timestep (m3; note: Interception and Throughfall are in mm)

    float transpiration_1016; //!< Variable used only to compare with eddy-flux tower data.
    float abund_phi_root; //!< Summary statistic output: abundance-weighted phi-root
    float abund10_phi_root; //!< Summary statistic output: abundance-weighted phi-root, including only tree with dbh>10cm
    float agb_phi_root; //!< Summary statistic output: agb-weighted phi-root
};

//! Counters and sums accumulated for outputs
struct Diagnostics
{
    int nblivetrees; //!< For output -- number of live trees at each timestep
    int nbtrees_n10, nbtrees_n30; //!< For output -- number of trees dbh > 10 cm, computed at beginning of each timestep
    int nbdead_n1, nbdead_n10, nbdead_n30; //!< For output -- number of deaths dbh > 1 cm, computed at each timestep
    int nbTreefall1, nbTreefall10, nbTreefall30; //!< For output -- number of treefalls at each timestep (dbh > 1cm), _BASICTREEFALL
    int nbtrees_carbstarv_n1, nbtrees_carbstarv_n10, nbtrees_carbstarv_n30; //!< For output -- number of trees > 1 cm with carbon starvation

    std::vector<int> nbdbh; //!< Number of trees per 1 cm DBH class, for outputs
    std::vector<float> layer; //!< Plot sum of LAI3D per height layer, for outputs (m2 m-2 summed over sites)
};

//! Settings of the simulated lidar point clouds
struct PointCloud
{
    float mean_beam_pc; //!< The mean number of shots per m2 for the point cloud sampling
    float sd_beam_pc; //!< The standard deviation of the shots per m2 for the point cloud sampling
    float klaser_pc; //!< The klight for the point cloud sampling, needs to be multiplied with the transmittance_laser to get effective klaser (i.e. 0.5 * 0.4 = 0.2)
    float transmittance_laser; //!< Transmittance of the laser when hitting something (or probability of continuing)
    int iter_pointcloud_generation; //!< Iteration at which point cloud should be generated
};

//! Random number generator and multivariate distribution of leaf traits (N, P, LMA)
struct RNGState
{
    gsl_rng *gslrand = nullptr; //!< Random number generator
    gsl_matrix *mcov_N_P_LMA = nullptr; //!< Covariance matrix for leaf_properties
    gsl_vector *mu_N_P_LMA = nullptr; //!< Mean values of the distributions and the output vector for the multivariate draw
    gsl_vector *variation_N_P_LMA = nullptr; //!< Mean values of the distributions and the output vector for the multivariate draw
    int covariance_status; //!< Covariance status: if one of N, P, or LMA has zero variation, the Cholesky decomposition fails, we then use no correlation at all
};

//! Seed rain and recruitment state per site and species
struct SpeciesState
{
    std::vector<std::vector<int>> SPECIES_SEEDS; //!< Seeds per site and species: presence (0/1), or number of seeds with the SEEDTRADEOFF option
    std::vector<double> p_seed; //!< Probability of each site receiving a seed from the external seed rain (multinomial weights, uniform)
    std::vector<unsigned int> n_seed; //!< Number of external seeds per site in the current timestep (multinomial draw)

    std::vector<double> p_species; //!< Relative frequency of each species in the external seed rain (multinomial weights)
    std::vector<unsigned int> n_species; //!< Number of external seeds per species in the current timestep (multinomial draw)

    std::vector<int> SPECIES_GERM; //!< Species with at least one seed at the current site, from which the recruit is drawn
    std::vector<float> PROB_S; //!< Per-species seed probabilities; allocated but not read by the current code (candidate for removal)
};

//! Intraspecific trait variation: standard deviations, correlations and precomputed distributions (10000 variants)
struct Intraspecific
{
    float d_intraspecific_height[10000]; //!< distribution of intraspecific values for maximal tree height
    float d_intraspecific_CR[10000]; //!< distribution of intraspecific values for maximal crown radius
    float d_intraspecific_CD[10000]; //!< distribution of intraspecific values for maximal crown depth
    float d_intraspecific_P[10000]; //!< distribution of intraspecific values for leaf phosphorus
    float d_intraspecific_N[10000]; //!< distribution of intraspecific values for leaf nitrogen
    float d_intraspecific_LMA[10000]; //!< distribution of intraspecific values for leaf mass per area (LMA)
    float d_intraspecific_wsg[10000]; //!< distribution of intraspecific values for wood density (wsg)
    float d_intraspecific_dbhmax[10000]; //!< distribution of intraspecific values for maximal trunk diameter
    float d_intraspecific_leafarea[10000]; //!< distribution of intraspecific values for leaf area (LA)
    float d_intraspecific_tlp[10000]; //!< distribution of intraspecific values for turgor loss point (TLP)

    // Standard deviations and correlations for intraspecific variation
    float sigma_height; //!< Standard deviation of intraspecific variation in maximal height
    float sigma_CR; //!< Standard deviation of intraspecific variation in crown radius
    float sigma_CD; //!< Standard deviation of intraspecific variation in crown depth
    float sigma_P; //!< Standard deviation of intraspecific variation in leaf phosphorus concentration
    float sigma_N; //!< Standard deviation of intraspecific variation in leaf nitrogen concentration
    float sigma_LMA; //!< Standard deviation of intraspecific variation in leaf mass per area (LMA)
    float sigma_wsg; //!< Standard deviation of intraspecific variation in wood specific gravity (wsg)
    float sigma_dbhmax; //!< Standard deviation of intraspecific variation in maximal trunk diameter
    float sigma_leafarea; //!< Standard deviation of intraspecific variation in leaf area (LA)
    float sigma_tlp; //!< Standard deviation of intraspecific variation in turgor loss point (TLP)

    float corr_CR_height; //!< Intraspecific correlation between crown radius and maximal height
    float corr_N_P; //!< Intraspecific correlation between leaf nitrogen and leaf phosphorus
    float corr_N_LMA; //!< Intraspecific correlation between leaf nitrogen and LMA
    float corr_P_LMA; //!< Intraspecific correlation between leaf phosphorus and LMA
    float cov_N_P; //!< Intraspecific covariance between leaf nitrogen and leaf phosphorus
    float cov_N_LMA; //!< Intraspecific covariance between leaf nitrogen and LMA
    float cov_P_LMA; //!< Intraspecific covariance between leaf phosphorus and LMA
};

//! Crown geometry parameters
struct CrownGeometry
{
    float crown_gap_fraction; //!< Fraction of gaps in a crown, new in v.2.4.1: adds gaps according to a predetermined gap_fraction in tree crown. Gaps inserted regularly (i.e. every 4 pixels for a gap fraction of 0.25). In v.2.5:  gap fraction is now given as a target fraction and, in filling up a tree's crown, gaps are inserted when fraction_filled is too high: more accurate approximation of gap fraction (e.g., gap fractions > 0.5 are now possible). Also in v.2.5: intraspecific variation in crown radius affects gap fraction: trees with larger crowns also have more gaps, so large crowns get more open and smaller crowns denser, while keeping the same leaf area. !!!IMPROVE_IN_FUTURE
    float shape_crown; //!< Percentage of crown radius that reaches the crown top and defines how steep the crown slope is; only required if CROWN_UMBRELLA is activated, but added to main code so that input sheet does not have to be modified.
    float Rndd, deltaR, deltaD; //!< Distance beyond which negative density dependence (NDD) effect is not accounted anymore

    float BAtot; //!< !!!UPDATE

    int extent_visual; //!< Extent of visualization outputs, bounded by extent of simulation
    int mincol_visual, maxcol_visual; //!< ! Global variables: computed from extent, defining the two visualization outputs
    int minrow_visual, maxrow_visual; //!< ! Global variables: computed from extent, defining the two visualization outputs
    int minrow_visual_slice, maxrow_visual_slice; //!< ! Global variables: computed from extent, defining the two visualization outputs
};

//! Global model parameters read from the global parameter file
struct SimParams
{
    // Light and photosynthesis
    float p_nonvert; //!< Ratio of non-vertical incident light
    float SWtoPPFD; //!< Conversion factor for shortwave irradiance measured in W/m2 to PPFD in micromol of PAR (micromol/s/m^2, as used in the Farquhar model). Around 2.0-2.5 in the tropics. Data at Nouragues (comparing photon count and irradiance) give a value: 2.27. Depends on cloudiness (in non-cloudy areas, the fraction of PAR in irradiance arriving on the ground may be much lower). This is typically equals to ca 0.5*4.57, where 0.5 stands for the fact that ca 50% of the total solar energy reaching the Earth’s surface corresponds to PAR, and 4.57 umol J-1 is basically equal to 10^6/(Emean*Avogadro number) where Emean is the average energy of a photon across PAR wavelength range: Ephoton=h*c/wavelength).
    float PPFDtoSW; //!< Conversion factor from PPFD in micromol of PAR to shortwave irradiance measured in J/ Inverse of SWtoPPFD.
    float klight; //!< Light absorption rate or extinction coefficient used in Beer-Lambert law to compute light within the canopy
    float kpar; //!< New in v.2.5: effective light absorption rate or extinction coefficient used in Beer-Lambert law to compute light within the canopy; kpar = klight * absorptance_leaves
    float phi; //!< True quantum yield; previously was apparent quantum yield (in micromol C/micromol photon): quantum yield multiplied by leaf absorptance. Quantum yield is often provided per absorbed light flux, so one should multiply incident PPFD by leaf absorptance (Poorter et al American Journal of Botany). For some authors, it should be species-dependent or environment dependent, but these options are not implemented here (see eg Domingues et al 2014 Plant Ecology & Diversity). As of v2.5: phi is the quantum yield, without multiplication with absorptance (see absorptance_leaves)
    float theta; //!< Parameter of the Farquhar model set to 0.7 in this version. For some authors, it should be species-dependent or environment dependent, but these options are not implemented here
    float absorptance_leaves; //!< Absorptance of leaves (close to 0.91 for tropical tree species)
    float alpha; //!< Apparent quantum yield to electron transport in mol e-/mol photons, equal to the true quantum yield multiplied by leaf absorbance
    float g1; //!< G1 parameter of Medlyn et al's model of stomatal conductance. v230: defined as a global parameter shared by species, instead of a Species class variable. !!!UPDATE A species-specific value of g1 (cf Lin et al. 2015 NCC or Wu et al. 2019 GCB) is needed to simulate functional shift through a regeneration. Different values of g1 across PFT are also used by Xu et al. 2016 New Phytologist using ED2+SPA applied on tropical dry forests
    float g0; //!< Minimum leaf conductance, in mmolH20 m-2 s-1, parameter of stomatal conductance model. Defined as a global parameter shared by species, in absence of a clear variation with other traits (see Duursma et al. 2018 New Phytologist, Slot et al. 2021 New Phytologist, METRADICA data)

    // CO2 / atmosphere
    float Cair; //!< Atmospheric CO_2 concentration (in ppm). If in the future CO_2 is allowed to vary, Cair should have the same status as other climatic variables
    float iCair; //!< Inverse of Cair
    float PRESS; //!< Atmospheric CO_2 pressure (in kPa). If in the future PRESS should vary, it should have the same status as other climatic variables

    // Phenology drought response
    float pheno_a0; //!< Threshold beyond which a change in t_pheno_factor (that controls old leaves shedding) is triggered, in proportion of leaf TLP
    float pheno_b0; //!< Threshold beyond which a change in t_pheno_factor (that controls old leaves shedding) is triggered, in proportion of tree height
    float pheno_delta; //!< Amplitude of change in t_pheno_factor per timestep

    // Tree allometry
    float H0; //!< Initial tree height (m)
    float DBH0; //!< Initial tree DBH (m)
    float CD0; //!< Initial tree crown depth (m)
    float CD_a; //!< Allometric parameter crown depth intercept
    float CD_b; //!< Allometric parameter crown depth slope
    float CR_a; //!< Allometric parameter crown radius log scale intercept (translates into factor on regular scale)
    float CR_b; //!< Allometric parameter crown radius log scale slope (translates into exponent on regular scale)
    float CR_min; //!< Allometric parameter minimum crown radius (in m)
    float dens; //!< Initial crown leaf density (in m^2/m^3)
    float fallocwood; //!< Fraction of biomass allocated to above ground wood (branch turnover+stem)
    float falloccanopy; //!< Fraction of biomass allocated to canopy (leaves + reproductive organs + twigs)
    float vC; //!< Variance of treefall threshold

    // Mortality
    float m; //!< Basal death rate
    float m1; //!< Slope of the dependence between death rate and wood density (wsg); usually a negative constant (v.2.2)
    float p_tfsecondary; //!< Probability that a death due to a treefall is a treefall itself (v.2.4.0)
    float hurt_decay; //!< "healing factor" rate at which t_hurt, the tree-level negative impact of a treefall on a tree, declines each timestep (v.2.4.0)

    // Seed dispersal
    float Cseedrain; //!< Constant used to scale total seed rain per hectare across species
    float nbs0; //!< Number of seeds produced and dispersed by each mature tree when SEEDTRADEOFF is not defined
};

//! Spatial fields: cumulated leaf area density and treefall damage
struct SimFields
{
    std::vector<std::vector<float>> LAI3D; //!< Leaf area index above each height (cumulated from the canopy top), per height layer and site including a border of SBORD sites (m2 m-2)
    std::vector<unsigned short> Thurt[3]; //!< Treefall field
};

//! All state of a TROLL simulation, passed by reference to the model functions (replaces the global variables of TROLL 4.0)
struct Context
{
    FileIO fileio; //!< Input and output file names
    InputBuffers buffers; //!< Raw command-line arguments, run index and seed
    ModelOptions opt; //!< Model options read from the global parameter file
    Grid grid; //!< Spatial grid: plot size, layers, soil cells (dcells)
    TimeState time; //!< Time stepping: iteration, timestep, days
    Climate climate; //!< Climate forcing: daily and half-hourly series
    LookupTables lookup; //!< Precomputed lookup tables (temperature responses, light, wind, LAImax)
    Soil soil; //!< Soil hydrology per layer and dcell
    Diagnostics diag; //!< Diagnostics accumulated for outputs
    PointCloud pc; //!< Point-cloud (lidar) simulation settings
    RNGState rng; //!< Random number generator and multivariate trait distribution
    SpeciesState species; //!< Seed rain and recruitment state
    Intraspecific intra; //!< Intraspecific trait variation: standard deviations, correlations and precomputed distributions
    CrownGeometry crown; //!< Crown geometry parameters
    SimParams params; //!< Global model parameters (light, photosynthesis, allometry, mortality, dispersal)
    SimFields field; //!< Spatial fields: leaf area density and treefall
    std::vector<Species> S; //!< All species
    std::vector<Tree> T;    //!< All tree sites
#ifdef WATER
    std::vector<RespGPPResult> resp_gpp; //!< Per-site fluxes precomputed in parallel each timestep (see Evolution)
#endif
    OutputConfig out;   // last: large fstream objects, accessed rarely
};

#endif