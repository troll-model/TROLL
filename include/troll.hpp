//////////////////////////////////////////////////////////////////////////////////
/*! \mainpage
 *  \section TROLL
 *
 *  Individual-based forest dynamics simulator
 *
 * History
 * -------
 *    - Version 1: Jerome Chave (Sept 1997 to Sept 2001; Chave 1999, 2001)
 *    - Version 2.0: Jerome Chave (March 2011)
 *    - Version 2.1 & 2.2: Isabelle Marechaux & Jerome Chave (Oct 2013 to May 2016)
 *    - Version 2.3: Isabelle Marechaux, Fabian Fischer, Jerome Chave (Oct 2016 to March 2017; Marechaux & Chave 2017)
 *    - Version 2.4 & 2.5: Fabian Fischer (Feb 2018 to May 2020)
 *    - Version 3.0: Fabian Fischer, Isabelle Marechaux, Jerome Chave (Jan 2021)
 *    - Version 3.1: Isabelle Maréchaux (May-June 2021)
 *    - Version 4.0: Isabelle Maréchaux (Feb 2023)
 *
 * Compiling options
 * -----------------
 *    - Compilation test on osx and linux plateforms
 *    - Compilation command: _g++ main.cpp -O3 -o TROLL.out -lgsl -lgslcblas -Wall_
 *    - Code profiling: _g++ -O3 -Wall -o troll main.cpp -lgsl -lgslcblas -lm -g -pg_
 *
 * Bibliography
 * -----------
 *    - Chave, J. (1999). Study of structural, successional and spatial patterns in tropical rain forests using TROLL, a spatially explicit forest model. Ecological modelling, 124(2-3), 233-254.
 *    - Chave, J. (2001). Spatial patterns and persistence of woody plant species in ecological communities. The American Naturalist, 157(1), 51-65.
 *    - Maréchaux, I., & Chave, J. (2017). An individual‐based forest model to jointly simulate carbon and tree diversity in Amazonia: description and applications. Ecological Monographs, 87(4), 632-664.
 *
 * Comments
 * -------
 *     - TODO: include iomanip, typedef, be careful with parameter_names initialisation
 *     - As of v.2.4, GNU scientific library (gsl) is needed (current version is 2.6)
 *        -# to work, the gsl_linalg.h header is needed, so versions above 1.7 are the minimum
 *        -# on osx, type "brew install gsl"
 *     - Basic uniform sampling functions (for more information: c:
 *        -# gsl_rng_uniform(gsl_rand) which samples range [0,1)
 *        -# gsl_rng_uniform_pos(gsl_rand) which samples range (0,1) to avoid 0 for log transformations
 *        -# gls_rng_uniform_int(gsl_rand, unsigned long int n) samples long integers from 0 to n-1
 */
//////////////////////////////////////////////////////////////////////////////////

#include "troll_defines.hpp"

// LIBRAIRIES
#include <cstdio>
#include <iostream>
#include <iomanip> //new IM (for setprecision)
#include <fstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <random>

#ifdef MPI
#include "mpi.h"
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

using namespace std;

#include "context.hpp"

// FileIO and InputBuffers migrated to ctx.fileio / ctx.buffers
// FILE OUTPUT STREAMS migrated to ctx.out (OutputConfig in context.hpp)

// USER CONTROLS migrated to ctx.opt:
// _NONRANDOM, _GPPcrown, _BASICTREEFALL, _SEEDTRADEOFF, _NDD, _CROWN_MM,
// _OUTPUT_extended, _OUTPUT_inventory, _FromInventory, _sapwood,
// _seedsadditional, _LL_parameterization, _LA_regulation,
// _OUTPUT_pointcloud, _SOIL_LAYER_WEIGHT, _WATER_RETENTION_CURVE

// GLOBAL PARAMETERS OF THE SIMULATION
// (sites, cols, rows, nbspp migrated to ctx.grid)

// RNG and covariance machinery migrated to ctx.rng:
// gslrand, mcov_N_P_LMA, mu_N_P_LMA, variation_N_P_LMA, covariance_status

// DCELL/Grid variables migrated to ctx.grid:
// length_dcell, linear_nb_dcells, sites_per_dcell, nbdcells,
// site_DCELL, i_sites_per_dcell, HEIGHT, dbhmaxincm, RMAX, SBORD,
// leafdem_resolution, NV, NH, LV, LH

// p_nonvert, Cseedrain, nbs0, Cair, iCair, PRESS migrated to ctx.params

// CrownGeometry migrated to ctx.crown:
// crown_gap_fraction, shape_crown, Rndd, deltaR, deltaD, BAtot,
// extent_visual, mincol_visual, maxcol_visual, minrow_visual,
// maxrow_visual, minrow_visual_slice, maxrow_visual_slice

//! ENVIRONMENTAL VARIABLES

// Climate vectors migrated to ctx.climate:
// varday_light, varday_vpd, varday_T, varday_WS,
// DailyMeanTemperature, DailyMeanIrradiance,
// DailyMeanVapourPressureDeficit, NightTemperature,
// Rainfall, DailyMeanWindSpeed

// LOOKUP TABLES migrated to ctx.lookup:
// nbTbins, iTaccuracy, LookUp_KmT/GammaT/VcmaxT/JmaxT/Rleaf,
// LookUp_flux_absorption/flux/ExtinctLW/VPD/T/Rstem/Crown_site,
// nbVPDbins, iVPDaccuracy, LookUp_INLR/SLOPE/GRADN,
// nbHbins, iHaccuracy, LookUp_Wind, LookUpLAImax

// ENVIRONMENTAL VARIABLES migrated to ctx.climate:
// tnight, precip, WSDailyMean, WDailyMean, tDailyMean, VPDDailyMean,
// WDailyMean_year, tDailyMean_year, VPDDailyMean_year,
// windDailyMean_year, Tnight_year,
// WDailyMean_all, VPDDailyMean_all, tDailyMean_all, windDailyMean_all

// Simulation parameters migrated to ctx.params:
// SWtoPPFD, PPFDtoSW, klight, kpar, phi, theta, absorptance_leaves, g1, g0,
// pheno_a0, pheno_b0, pheno_delta, alpha, vC, H0, DBH0, CD0,
// fallocwood, falloccanopy, dens, CD_a, CD_b, CR_a, CR_b, CR_min,
// p_tfsecondary, hurt_decay, m, m1

// Intraspecific variation parameters migrated to ctx.intra:
// sigma_height/CR/CD/P/N/LMA/wsg/dbhmax/leafarea/tlp,
// corr_CR_height/N_P/N_LMA/P_LMA, cov_N_P/N_LMA/P_LMA

// Intraspecific arrays migrated to ctx.intra:
// d_intraspecific_height/CR/CD/P/N/LMA/wsg/dbhmax/leafarea/tlp

// LookUpLAImax migrated to ctx.lookup

// Simulation fields migrated to ctx.field:
// LAI3D, Thurt[3]

// Soil globals migrated to ctx.soil:
// nblayers_soil, layer_depth, Sat_SWC, Max_SWC, FC_SWC, Res_SWC, Min_SWC,
// Ksat, a_vgm, b_vgm, c_vgm, m_vgm, phi_e, b, SWC3D, soil_phi3D, Ks, KsPhi,
// LAI_DCELL, LAI_young, LAI_mature, LAI_old, Canopy_height_DCELL, HSum_DCELL,
// TopWindSpeed_DCELL, Interception, Throughfall, Runoff, Leakage, Evaporation,
// Transpiration, transpiration_1016, abund_phi_root, abund10_phi_root, agb_phi_root

// SEED DENSITIES migrated to ctx.species:
// SPECIES_SEEDS, p_seed, n_seed, p_species, n_species,
// SPECIES_GERM, PROB_S

// PointCloud globals migrated to ctx.pc:
// mean_beam_pc, sd_beam_pc, klaser_pc, transmittance_laser, iter_pointcloud_generation

#ifdef Output_ABC
// Output creation for Approximate Bayesian Computation and Summary Statistics calculation. ABC produces large files, so metrics are pre-calculated within TROLL
int chmchange_iter = 43; //!< Global variable: ad-hoc variable. Output at step 43 iterations before the final output, to simulate a suite of successive LiDAR measurements at Nouragues, from March 2012 to October 2015, so 3 * 12 + 7 = 43 !!!UPDATE
// the margins of the ABC analysis, can be smaller than the actually simulated plot to deal with edge effects

int margin;       //!< Global variable: for ABC -- margin used to deal with edge effects
int row_start;    //!< Global variable: for ABC -- row start after margin clipping
int row_end;      //!< Global variable: for ABC -- row end after margin clipping
int col_start;    //!< Global variable: for ABC -- col start after margin clipping
int col_end;      //!< Global variable: for ABC -- col end after margin clipping
int sites_abc;    //!< Global variable: for ABC -- number of sites after margin clipping
int nbvisited;    //!< Global variable: for ABC -- number of visited sites !!!UPDATE
int patch_size;   //!< Global variable: for ABC -- patch size !!!UPDATE
float isites_abc; //!< Global variable: for ABC -- inverse of number of sites after margin clipping

int *chm_field_previous(0);     //!< Global vector: for ABC -- canopy height model (chm) calculated at previous step !!!UPDATE
int *chm_field_current(0);      //!< Global vector: for ABC -- canopy height model (chm) calculated at current step !!!UPDATE
int *chm_field_previous_ALS(0); //!< Global vector: for ABC -- canopy height model (chm) inferred from Airborne Lidar Scanning (ALS) at previous step
int *chm_field_current_ALS(0);  //!< Global vector: for ABC -- canopy height model (chm) inferred from Airborne Lidar Scanning (ALS) at current step
int *chm_field_changes(0);      //!< Global vector: for ABC -- change in canopy height as calculated
int *chm_field_changes_ALS(0);  //!< Global vector: for ABC -- change in canopy height as observed
float **transmittance_simulatedALS(0);
float **transmittance_direct(0);
int **transmittance_simulatedALS_sampling(0);
float GPP_MA[120] = {0.0};        //!< Global variable: moving average of GPP across 120 iterations (originally 120 months, could also be made dependent on iterperyear)
float Litterfall_MA[120] = {0.0}; //!< Global variable: moving average of litterfall across 120 iterations (originally 120 months, could also be made dependent on iterperyear)
float Mortality_MA[120] = {0.0};  //!< Global variable: moving average of mortality across 120 iterations (originally 120 months, could also be made dependent on iterperyear)
float Treefall_MA[120] = {0.0};   //!< Global variable: moving average of treefall across 120 iterations (originally 120 months, could also be made dependent on iterperyear)
#endif

// DIAGNOSTICS OUTPUT STATISTICS migrated to ctx.diag:
// nblivetrees, nbtrees_n10/n30, nbdead_n1/n10/n30,
// nbTreefall1/10/30, nbtrees_carbstarv_n1/n10/n30, nbdbh, layer
// (Output_ABC globals omitted: Output_ABC is #undef'd)

#ifdef CHECK_CARBON
float carbon_assimilated_total; //!< Global variable: --  !!!UPDATE
float carbon_net_total;         //!< Global variable: --  !!!UPDATE
float carbon_stored_leaves;     //!< Global variable: --  !!!UPDATE
float carbon_stored_trunk;      //!< Global variable: --  !!!UPDATE
float carbon_stored_free;       //!< Global variable: --  !!!UPDATE
#endif

// PROCESSORS FOR MPI OPTION
extern int mpi_rank;     //!< Global variable: processor rank (Message Passing Interface)
extern int mpi_size;     //!< Global variable: processor number (Message Passing Interface)
extern int easympi_rank; //!< Global variable: processor rank (easy Message Passing Interface)

// MPI PROCEDURES
#ifdef MPI
unsigned short **LAIc[2];                                        //!< Global 3D field: sharing LAI across boundaries
void MPI_ShareSeed(unsigned char **, int);                       //!< Global MPI function: Communication of border fields in the parallel version of the code
void MPI_ShareField(Context &ctx, unsigned short **, unsigned short ***, int); //!< Global MPI function: Communication of fields
void MPI_ShareTreefall(unsigned short **, int);                  //!< Global MPI function: Communication of treefalls
#endif

#ifdef CROWN_UMBRELLA
// TREESHAPE CALCULATIONS
// All these calculations are currently defined as non-member functions, but could potentially be converted to class Tree
// See function description for more details
void GetPPFDabove(Context &ctx, int height, int site, float noinput, float (&PPFD)[2]); //!< Global function: PPFD retrieval for function leafarea_max()
#ifdef WATER
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[6]); //!< Global function: calculates the canopy environment // RENAMED: dens → dens_layer
#else
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[4]); //!< Global function: calculates the canopy environment // RENAMED: dens → dens_layer
#endif
void AddCrownVolumeLayer(Context &ctx, int row_center, int col_center, float height, float CR, float CD, int crownvolume[70]); //!< Global function: calculates packing densities
void UpdateLAI3D(Context &ctx, int height, int site, float dens_layer, float &LA_cumulated);                                    //!< Global function: update of LAI3D field, called by CalcLAI() // RENAMED: dens → dens_layer
#ifdef CHM_SPIKEFREE
void UpdateCHMvector(int height, int site, float noinput, vector<int> &chm); //!< Global function: remove outliers in canopy height model (CHM); vector option
void UpdateCHM(int height, int site, float noinput, int *chm);               //!< Global function: remove outliers in canopy height model (CHM)
#endif
void OutputCrownSliced(Context &ctx, int height, int site, int row_slice, vector<float> &output_statistics);      //!< Global function: write a slice of a crown to file
void KeepFloatAsIs(float input, float &output, float CD, float height, int layer_fromtop);          //!< Global function: dummy function when no modification is needed
void KeepIntAsIs(int input, int &output, float CD, float height, int layer_fromtop);                //!< Global function: dummy function when no modification is needed
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);             //!< Global function: a modifying function that converts LAI to the density of a specific layer, using the GetDensity functions
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);   //!< Global function: a modifying function that converts LAI to percentage LAI in and above the current layer, using the GetDensity functions; can be used to directly allocate LAI without looping over LAI3D field; new in v.3.1
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens_layer); //!< Global function: deduces within-crown densities from LAI // RENAMED: dens → dens_layer
void GetDensityUniform(float LAI, float CD, float &dens_layer);                                           //!< Global function: deduces within-crown density from LAI // RENAMED: dens → dens_layer
int GetCrownIntarea(float radius);                                                                  //!< Global function: converts floating point crown area into integer value, imposing lower and upper limits
float GetRadiusSlope(Context &ctx, float CR, float crown_extent, float crown_position);                           //!< Global function: linear decrease of crown radius
float GetRadiusCylinder(float CR, float crown_extent, float crown_position);                        //!< Global function: not currently used, but returns the input radius
template <typename I, typename O, typename F>
void CircleAreaUpdateCrownStatistic_template(Context &ctx, int row_center, int col_center, int pos_start, int pos_end, float fraction_filled_target, float &fraction_filled_actual, int height_layer, I CrownStatistic_input, O &CrownStatistic_output, F UpdateCrownStatistic)
{
    for (int i = pos_start; i < pos_end; i++)
    {
        if (fraction_filled_actual > fraction_filled_target)
        {
            fraction_filled_actual = (fraction_filled_actual * float(i)) / (float(i) + 1.0);
        }
        else
        {
            fraction_filled_actual = (fraction_filled_actual * float(i) + 1.0) / (float(i) + 1.0);

            int site_relative = ctx.lookup.LookUp_Crown_site[i];
            int row, col;
            row = row_center + site_relative / 51 - 25;
            col = col_center + site_relative % 51 - 25;

            if (row >= 0 && row < ctx.grid.rows && col >= 0 && col < ctx.grid.cols)
            {
                int site = col + ctx.grid.cols * row;
                UpdateCrownStatistic(height_layer, site, CrownStatistic_input, CrownStatistic_output);
            }
        }
    }
}
//!< Global function: Template function called by LoopLayerUpdateCrownStatistic_template

template <typename G, typename I, typename O, typename M, typename F>
void LoopLayerUpdateCrownStatistic_template(
    Context &ctx, // context (parameters)
    int row_center, int col_center, // voxel center
    float height, float CR, float CD, // tree dimensions
    float fraction_filled_target, 
    int shell_fromtop, 
    G GetRadiusLayer, 
    I CrownStatistic_input, 
    O &CrownStatistic_output, 
    M ModifyCrownStatistic_input, 
    F UpdateCrownStatistic_output
) {
    int crown_top = int(height);

    // we start out with 0 actually filled voxels. As a result, the first voxel will always be filled
    float fraction_filled_actual = 0.0;
    if (CD <= 3.0)
    {
        I CrownStatistic_input_modified;
        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_modified, CD, height, shell_fromtop);

        int crown_intarea_previous = 0;
        int crown_intarea = GetCrownIntarea(CR);
        int layer_cylinder = crown_top - shell_fromtop;
        CircleAreaUpdateCrownStatistic_template(ctx, row_center, col_center, crown_intarea_previous, crown_intarea, fraction_filled_target, fraction_filled_actual, layer_cylinder, CrownStatistic_input_modified, CrownStatistic_output, UpdateCrownStatistic_output);
    }
    else
    {
        // This function computes the extent of the crown at every height ctx.diag.layer, given a specific function
        // it separates out the innermost sector (a slowly increasing cylinder), and the surrounding parts of the crown
        // first the metrics with respect to the internal crown structure (i.e. z coordinate with respect to crown base)
        float crownshell_base = height - CD + 2.0;                   // lower reference point for the crown slope function is two layers up from the crown base
        float crownshell_extent = height - crownshell_base;          // this is the extent from the "base ctx.diag.layer" to the top
        float crownshell_extent_toplayer = floor(crownshell_extent); // this is the extent to the lower limit of the toplayer
        // then we translate the crown coordinates into discretised variables with respect to the absolute location in the voxel field, as needed for location in the voxel field, with layers defined from top to bottom
        int height_innermost = crown_top - shell_fromtop;
        int height_toplayer = int(crownshell_base + crownshell_extent_toplayer) - shell_fromtop;
        int height_baselayer = int(crownshell_base + 1.0) - shell_fromtop;

        // now calculate the two modifications of the input statistic
        I CrownStatistic_input_innermost;
        I CrownStatistic_input_outer;

        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_innermost, CD, height, shell_fromtop);
        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_outer, CD, crownshell_base, shell_fromtop);

        // now do calculations
        // first the inner crown shell section that grows dynamically
        int crown_intarea_previous = 0;

        float radius_innermost = GetRadiusLayer(CR, crownshell_extent, crownshell_extent_toplayer);
        int crown_intarea_innermost = GetCrownIntarea(radius_innermost);
        CircleAreaUpdateCrownStatistic_template(ctx, row_center, col_center, crown_intarea_previous, crown_intarea_innermost, fraction_filled_target, fraction_filled_actual, height_innermost, CrownStatistic_input_innermost, CrownStatistic_output, UpdateCrownStatistic_output);
        crown_intarea_previous = crown_intarea_innermost;

        // now loop through the outer crown shell cylinders
        for (int h_outer = height_toplayer; h_outer >= height_baselayer; h_outer--)
        {
            // calculating the radius of the current ctx.diag.layer depending on the respective slopes, to be replaced by function
            // float radius_height = CR - crown_slope * (h_outer - height_baselayer);    // for the lowest ctx.diag.layer, i.e. h == height_baselayer, radius = t_CR
            int extent_layerouter = h_outer - height_baselayer;
            float radius_height = GetRadiusLayer(CR, crownshell_extent, extent_layerouter);
            int crown_intarea = GetCrownIntarea(radius_height);
            CircleAreaUpdateCrownStatistic_template(ctx, row_center, col_center, crown_intarea_previous, crown_intarea, fraction_filled_target, fraction_filled_actual, h_outer, CrownStatistic_input_outer, CrownStatistic_output, UpdateCrownStatistic_output);
            crown_intarea_previous = crown_intarea;
        }
    }
}
//!< Global function: Template function to loop across crown layer and update
#endif

// Fields computation
float CalcVcmaxm(float lma, float nmass, float pmass);             //!< Returns Vcmaxm, in micromol C g-1 s-1
float CalcRdark(float lma, float nmass, float pmass, float vcmax); //!< Returns Rdark, in micromol C m-2 s-1
#ifdef G0
float QUAD(float A, float B, float C, int ROOT); //!< Returns solution of quadratic equation (larger root of ROOT=1, smaller root if ROOT=-1)
#endif
void CalcLAI(Context &ctx, Tree &tree); 

// GLOBAL FUNCTIONS
void AssignSpeciesParam(Species &S, const std::string &name, const std::string &value); //!< Global function: assign a species parameter by name
void AssignValuePointcloud(Context &ctx, std::string parameter_name, std::string parameter_value); //!< Global function: assign a pointcloud parameter by name
void ReadInputGeneral(Context &ctx);                                //!< Global function: read in global parameter sheet
void ReadInputSpecies(Context &ctx);                                //!< Global function: read in species parameter sheet
void ReadInputDailyvar(Context &ctx);                               //!< Global function: read in daily climate variation
void ReadInputClimate(Context &ctx);                                //!< Global function: read in climate data
void ReadInputSoil(Context &ctx);                                   //!< Global function: read in soil data
void Initialise(Context &ctx);                                      //!< Global function: initialisation with bare ground conditions
void InitialiseIntraspecific(Context &ctx);                         //!< Global function: initialise intraspecific variables
void InitialiseLookUpLAImax(Context &ctx);                          //!< Global function: initialise lookup table for LAImax
void InitialiseLookUpTables(Context &ctx);                          //!< Global function: initialise lookup tables
void InitialiseOutputStreams(Context &ctx);                         //!< Global function: initialisation of output streams
void ReadInputInventory(Context &ctx);                              //!< Global function: updated in v.3.1: initialisation from inventories
void ReadInputPointcloud(Context &ctx);                             //!< Global function: introduced in v.3.1.6: reads the parameter sheet for the point cloud simulation
void AllocMem(Context &ctx);                                        //!< Global function: Field dynamic memory allocation
void Evolution(Context &ctx);                                       //!< Global function: Evolution at each timestep
void UpdateSeeds(Context &ctx);                                     //!< Global function: Compute field Seed
void UpdateField(Context &ctx);                                     //!< Global function: Update all fields
void TriggerTreefall(Context &ctx);                                 //!< Global function: Treefall gap formation; v.2.4
void TriggerTreefallSecondary(Context &ctx);                        //!< Global function: Secondary treefall gap formation
void FillSeed(Context &ctx, int col, int row, int spp);                           //!< Global function: update SPECIES_SEEDS field; v.2.5
void RecruitTree(Context &ctx);                                     //!< Global function: tree germination module; v.2.5
void Average(Context &ctx);                                                        //!< Global function: output of the global averages every timestep
void OutputField(Context &ctx);                                                    //!< Global function: output of the field variables every timestep
void OutputSnapshot(Context &ctx, fstream &output, bool header, float dbh_limit); //!< Global function: output snapshots of the scene at one point in time
void OutputLAI(Context &ctx, fstream &output_transmLAI3D);                        //!< Global function: writes the whole 3D LAI voxel field to file
void OutputCHM(Context &ctx, fstream &output_CHM);                                //!< Global function: Outputs CHM
void OutputVisual(Context &ctx);                                                   //!< Global function: Output function for visualization purposes
void CloseOutputs(Context &ctx);
void FreeMem(Context &ctx);

void ExportPointcloud(float mean_beam, float sd_beam, float klaser, float transmittance_laser, fstream &output_pointcloud); //!< Global function: point cloud output, v.3.1.6; kept separately from other output functions, as we write to a dedicated external file format (.las)

// HELPER FUNCTIONS
int GetTimeofyear(Context &ctx);                              //!< Helper function, new in v.3.1: converts current iteration into time of year, also works backwards (negative iterations)
float CalcHeightBaseline(float &ah, float &hmax, float &dbh); //!< Helper function: calculates mean predicted height from allometry
float CalcCDBaseline(Context &ctx, float &height);            //!< Helper function, new in v.3.1: calculates mean predicted crown radius from allometry
float CalcCRBaseline(Context &ctx, float &dbh);               //!< Helper function new in v.3.1: calculates mean predicted crown diameter from allometry
int CalcIntabsorb(float absorb_prev, float absorb_delta);     //!< Helper function: returns index for LookUpTables of absorbed flux (considering leaves above and within voxel)
int CalcIntabsorb(float absorb_prev);                         //!< Helper function: returns index for LookUpTables of absorbed flux (only considering leaves above voxel)

// DATA READING FUNCTIONS
// update v.3.1: already defined here to enable calling from functions of class Tree()
template <typename N>
void SetParameter(string &parameter_name, string &parameter_value, N &parameter, N parameter_min, N parameter_max, N parameter_default, bool quiet);
void SetParameter(string &parameter_name, string &parameter_value, string &parameter, string parameter_default, bool quiet);

#ifdef Output_ABC
// these are functions needed to create ABC output
void InitialiseABC(Context &ctx);                                                                                                                                                                                                                                                                                                                                                                                              //!< Global ABC function: initialise ABC conditions
void UpdateMovingAveragesABC(Context &ctx);                                                                                                                                                                                                                                                                                                                                                                                //!< Global ABC function: yearly statistics
void UpdateDBHtrackingABC(Context &ctx);                                                                                                                                                                                                                                                                                                                                                                                   //!< Global ABC function: update DBH function for ABC routines
void UpdateTransmittanceCHM_ABC(int mean_beam, float sd_beam, float klaser, float transmittance_laser);                                                                                                                                                                                                                                                                                                                    //!< Global function: calculating the TROLL transmittance field from simulated LiDAR
void OutputABCWriteHeaders(fstream &output_traitconservation, fstream &output_field, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_transmittance, fstream &output_transmittance_ALS, fstream &output_LAIfield, fstream &output_LAIfield_ALS, fstream &output_chmpotential, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass); //!< Global function: write headers for ABC outputs
void OutputABCConservationTraits(Context &ctx, fstream &output_traitconservation);                                                                                                                                                                                                                                                                                                                                         //!< Global ABC function: assesses the conservation of traits between input and output in TROLL
void OutputABC_ground(Context &ctx, fstream &output_field);                                                                                                                                                                                                                                                                                                                                                                //!< Global ABC function: returns ABC ground outputs
void OutputABC_species(Context &ctx, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass);                                                                                                                                                                                                                                                       //!< Global ABC function: returns ABC outputs for species
void OutputABC_CHM(Context &ctx, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_chmpotential);                                                                                                                                                                                                                                                                                                              //!< Global function: returns ABC outputs for canopy height model (CHM)
void OutputABC_transmittance(Context &ctx, fstream &output_transmittance, fstream &output_transmittance_ALS);                                                                                                                                                                                                                                                                                                              //!< Global ABC function: ABC outputs
void OutputABC(Context &ctx);                                                                                                                                                                                                                                                                                                                                                                                              //!< Global ABC function: output general ABC statistics
void OutputCHM(Context &ctx, fstream &output_CHM);                                                                                                                                                                                                                                                                                                                                                                        //!< Global ABC function: creates a Canopy Height Model and LAD profile
#endif

#ifdef TRACK_INDIVIDUALS
void TrackingData_andOutput(Context &ctx); //!< Global function: tree level tracking of key variables
#endif

// Species and Tree class declarations moved to include/species.hpp and include/tree.hpp
// (included transitively via context.hpp)

//// UTIL Templates

template <typename N>
bool AreEqual(N value1, N value2, int precision)
{
    return abs(value1 - value2) < pow(10, -precision); // make sure that namespace std is activated, otherwise abs is only defined for integers
}

inline string GetParameter(string &parameter_name, vector<string> &parameter_names, vector<string> &parameter_values)
{
    unsigned int parameter_index = 0;
    while (parameter_index < parameter_names.size() && parameter_names[parameter_index] != parameter_name)
    {
        parameter_index++;
    }

    /* if (parameter_name=="s_name") {
         for (int i=0; i<=parameter_names.size(); i=i+1) cout << parameter_names[i] << "=" << parameter_values[i] << "\t";
         cout << endl;
         cout << parameter_names[parameter_index] << " ="  << parameter_values[parameter_index] << endl;
         cout << parameter_names.size() << "\t"  << parameter_values.size() << endl;

     }*/
    string parameter_value;
    if (parameter_index < parameter_names.size())
        parameter_value = parameter_values[parameter_index];
    else
    {
        parameter_value = ""; // should not happen, potentially throw warning
        cout << "Problem finding parameter " << parameter_name << " please check input string in parameter definitions." << endl;
    }
    return (parameter_value);
}

template <typename N>
void CompareParameters(N &parameter_value, N parameter_min, N parameter_max, N parameter_default, string &parameter_name, vector<string> &parameter_names, vector<string> &parameter_values, bool quiet)
{
    // function could be used to systematically compare empirically provided values and computed values
    // since this could result in a lot of computation and output during initialisation, we rather print out a Snapshot of the initial configuration, to be compared to input
    N parameter_value_emp;
    string parameter_value_emp_string = GetParameter(parameter_name, parameter_names, parameter_values);
    SetParameter(parameter_name, parameter_value_emp_string, parameter_value_emp, parameter_min, parameter_max, parameter_default, quiet);
    if (parameter_value_emp != parameter_default)
    {
        int precision = 4;
        bool equal = AreEqual(parameter_value, parameter_value_emp, precision);
        if (!equal)
            cout << "Supplied parameter " << parameter_name << " does not correspond to calculated value at precision " << precision << endl;
    }
}
