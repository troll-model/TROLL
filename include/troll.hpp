//////////////////////////////////////////////////////////////////////////////////
/*! \page history History and bibliography
 *  \section troll_history TROLL
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
 * Building
 * --------
 *    - See README.md: TROLL is built with CMake since 4.1.0.
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
unsigned short **LAIc[2]; //!< Global 3D field: sharing LAI across boundaries
/**
 * @brief Communication of border fields in the parallel version of the code.
 * @details Only if the MPI option has been enabled.
 */
void MPI_ShareSeed(unsigned char **, int);
/**
 * @brief Communication of fields.
 */
void MPI_ShareField(Context &ctx, unsigned short **, unsigned short ***, int);
/**
 * @brief Communication of treefalls.
 */
void MPI_ShareTreefall(unsigned short **, int);
#endif

#ifdef CROWN_UMBRELLA
// TREESHAPE CALCULATIONS
// All these calculations are currently defined as non-member functions, but could potentially be converted to class Tree
// See function description for more details
/**
 * @brief PPFD retrieval for function leafarea_max().
 * @param ctx simulation context (all model state)
 * @param height tree height or height layer (m)
 * @param site site index
 * @param noinput unused (common signature of crown statistics)
 * @param ppfd_CA [in,out] accumulator: [0] PPFD summed over the crown area, [1] number of sites summed
 */
void GetPPFDabove(Context &ctx, int height, int site, float noinput, float (&ppfd_CA)[2]);
#ifdef WATER
/**
 * @brief Calculates the canopy environment.
 * @param ctx simulation context (all model state)
 * @param height tree height or height layer (m)
 * @param site site index
 * @param dens_layer leaf area density of the crown in the layer (m2 m-2 per 1 m layer); [out] for GetDensity*
 * @param canopy_environment_cumulated [in,out] accumulators over the crown area: light, VPD, temperature decrease, leaf area (and with WATER: incident light, longwave extinction)
 * @details This function adds to the environmental variables provided in canopy_environment_cumulated: the PPFD, VPD, Tmp and leafarea_layer retrieval function for Fluxh().
 */
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[6]); // RENAMED: dens → dens_layer
#else
/**
 * @brief Calculates the canopy environment.
 * @param ctx simulation context (all model state)
 * @param height tree height or height layer (m)
 * @param site site index
 * @param dens_layer leaf area density of the crown in the layer (m2 m-2 per 1 m layer); [out] for GetDensity*
 * @param canopy_environment_cumulated [in,out] accumulators over the crown area: light, VPD, temperature decrease, leaf area (and with WATER: incident light, longwave extinction)
 * @details This function adds to the environmental variables provided in canopy_environment_cumulated: the PPFD, VPD, Tmp and leafarea_layer retrieval function for Fluxh().
 */
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[4]); // RENAMED: dens → dens_layer
#endif
/**
 * @brief Calculates packing densities.
 */
void AddCrownVolumeLayer(Context &ctx, int row_center, int col_center, float height, float CR, float CD, int crownvolume[70]);
/**
 * @brief Update of LAI3D field, called by CalcLAI().
 * @param ctx           simulation context
 * @param height        height layer (m)
 * @param site          site index
 * @param dens_layer    leaf area density of the crown in the layer (m2 m-2 per 1 m layer)
 * @param LA_cumulated  [in,out] leaf area added to the LAI3D field so far (m2)
 */
void UpdateLAI3D(Context &ctx, int height, int site, float dens_layer, float &LA_cumulated); // RENAMED: dens → dens_layer
#ifdef CHM_SPIKEFREE
/**
 * @brief Remove outliers in canopy height model (CHM); vector option.
 */
void UpdateCHMvector(int height, int site, float noinput, vector<int> &chm);
/**
 * @brief Remove outliers in canopy height model (CHM).
 */
void UpdateCHM(int height, int site, float noinput, int *chm);
#endif
/**
 * @brief Write a slice of a crown to file.
 */
void OutputCrownSliced(Context &ctx, int height, int site, int row_slice, vector<float> &output_statistics);
/**
 * @brief Dummy function when no modification is needed.
 */
void KeepFloatAsIs(float input, float &output, float CD, float height, int layer_fromtop);
/**
 * @brief Dummy function when no modification is needed.
 */
void KeepIntAsIs(int input, int &output, float CD, float height, int layer_fromtop);
/**
 * @brief A modifying function that converts LAI to the density of a specific layer, using the GetDensity functions.
 * @details Modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results. LAI is the input, dens_layer the output.
 */
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);
/**
 * @brief A modifying function that converts LAI to percentage LAI in and above the current layer, using the GetDensity functions; can be used to directly allocate LAI without looping over LAI3D field; new in v.3.1.
 * @details Modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results. LAI is the input, dens_layer the output.
 */
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);
/**
 * @brief Deduces within-crown densities from LAI.
 * @param LAI            leaf area index of the crown (m2 m-2)
 * @param CD             crown depth (m)
 * @param dens_top       [out] leaf area density of the top layer
 * @param dens_belowtop  [out] leaf area density of the layer below the top
 * @param dens_layer     [out] leaf area density of the lower layers
 */
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens_layer); // RENAMED: dens → dens_layer
/**
 * @brief Deduces within-crown density from LAI.
 * @param LAI leaf area index of the crown (m2 m-2)
 * @param CD crown depth (m)
 * @param dens_layer leaf area density of the crown in the layer (m2 m-2 per 1 m layer); [out] for GetDensity*
 */
void GetDensityUniform(float LAI, float CD, float &dens_layer); // RENAMED: dens → dens_layer
/**
 * @brief Converts floating point crown area into integer value, imposing lower and upper limits.
 */
int GetCrownIntarea(float radius);
/**
 * @brief Linear decrease of crown radius.
 */
float GetRadiusSlope(Context &ctx, float CR, float crown_extent, float crown_position);
/**
 * @brief Not currently used, but returns the input radius.
 */
float GetRadiusCylinder(float CR, float crown_extent, float crown_position);
/**
 * @brief Template function called by LoopLayerUpdateCrownStatistic_template.
 *
 * Fill one crown layer site by site, in the precomputed spiral order (LookUp_Crown_site), skipping
 * sites so that the filled fraction tracks fraction_filled_target, and apply UpdateCrownStatistic
 * to each filled site.
 * @param ctx                     simulation context
 * @param row_center              row of the crown centre
 * @param col_center              column of the crown centre
 * @param pos_start               first position in the spiral order
 * @param pos_end                 position after the last one
 * @param fraction_filled_target  target fraction of the crown area filled with leaves (1 - crown gap fraction)
 * @param fraction_filled_actual  [in,out] fraction filled so far
 * @param height_layer            height layer (m)
 * @param CrownStatistic_input    value passed to UpdateCrownStatistic for each site (e.g. leaf area density of the layer)
 * @param CrownStatistic_output   [in,out] statistic accumulated over the sites
 * @param UpdateCrownStatistic    function applied to each site (e.g. UpdateLAI3D, GetCanopyEnvironment)
 */
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
/**
 * @brief Template function to loop across crown layer and update.
 *
 * Apply a crown statistic to every layer of a crown: the crown is processed as shells from the top,
 * each layer filled by CircleAreaUpdateCrownStatistic_template with the radius given by GetRadiusLayer.
 * @param ctx                          simulation context
 * @param row_center                   row of the crown centre
 * @param col_center                   column of the crown centre
 * @param height                       tree height (m)
 * @param CR                           crown radius (m)
 * @param CD                           crown depth (m)
 * @param fraction_filled_target       target fraction of the crown area filled with leaves
 * @param shell_fromtop                offset of the processed shell below the crown top (layers)
 * @param GetRadiusLayer               radius of a layer from the crown shape (GetRadiusSlope or GetRadiusCylinder)
 * @param CrownStatistic_input         crown-level input (e.g. crown LAI)
 * @param CrownStatistic_output        [in,out] statistic accumulated over the crown
 * @param ModifyCrownStatistic_input   converts the input for each layer (e.g. LAI2dens, KeepFloatAsIs)
 * @param UpdateCrownStatistic_output  function applied to each site of each layer
 */
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
#endif

/**
 * @brief Returns Vcmaxm, in micromol C g-1 s-1.
 */
float CalcVcmaxm(float LMA, float Nmass, float Pmass);
/**
 * @brief Returns Rdark, in micromol C m-2 s-1.
 */
float CalcRdark(float lma, float nmass, float pmass, float vcmax);
#ifdef G0
/**
 * @brief Returns solution of quadratic equation (larger root of ROOT=1, smaller root if ROOT=-1).
 */
float QUAD(float A, float B, float C, int ROOT);
#endif
/**
 * @brief Fills LAI3D within the crown of a tree.
 * @details Called by UpdateField. Modified in v.2.3: additional contribution to voxels that are not fully occupied by the tree crown. !!!: this does not calculate ctx.field.LAI3D directly, this only calculates the density in each voxel belonging to a tree. The final LAI field is calculated outside of the class Tree. Modified in v.2.4 and v.2.5: introducing an alternative crown shape, "umbrella"-like, inspired by previous shell models and similar to the crown shapes in the PPA. If activated, crowns contain three layers of vegetation that, once the crown goes beyond 3m in depth, will bend downwards on the edges with a linear slope. Since v.2.5 all loops (CalcLAI, Fluxh, leafarea_max) are executed through the same template. This allows to implement other crown shapes in the future and ensures that modifications are carried through across the code.
 * @param ctx simulation context (all model state)
 * @param tree reference to tree object
 */
void CalcLAI(Context &ctx, Tree &tree);
float GetCrownAreaFilled(float crown_area, float fraction_filled_target);

// GLOBAL FUNCTIONS
/**
 * @brief Assign a species parameter by name.
 * @details Global function: This function provides limits and defaults for species-specific parameters. Register the species-file columns (name, type, allowed range) used to read species parameters.
 * @param S [in,out] species whose parameter is set
 * @param name column name in the species file
 * @param value value read for this species (text)
 */
void AssignSpeciesParam(Species &S, const std::string &name, const std::string &value);
/**
 * @brief Assign a pointcloud parameter by name.
 * @param ctx simulation context (all model state)
 * @param parameter_name parameter name
 * @param parameter_value value read from the file (text)
 */
void AssignValuePointcloud(Context &ctx, std::string parameter_name, std::string parameter_value);
/**
 * @brief Read in global parameter sheet.
 * @details Global function: This function reads inputs from the general input file. Also sets up initial values for some parameters computed from input parameters.
 * - !!!: TO IMPLEMENT: automatic stop of program under error to avoid running empty simulations on cluster etc.
 * - !!!: TO IMPLEMENT: separate error messages (i.e. no input file provided, empty input file, etc.)
 * @param ctx simulation context (all model state)
 */
void ReadInputGeneral(Context &ctx);
/**
 * @brief Read in species parameter sheet.
 * @details Global function: This function reads inputs from the species input file.
 * @param ctx simulation context (all model state)
 */
void ReadInputSpecies(Context &ctx);
/**
 * @brief Read in daily climate variation.
 * @details Global function: This function reads inputs from the environmental daily variation input file.
 * @param ctx simulation context (all model state)
 */
void ReadInputDailyvar(Context &ctx);
/**
 * @brief Read in climate data.
 * @details Global function: This function reads inputs from the environmental variation input file.
 * - v.3.0 of the code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=ctx.time.iter and not ctx.time.iterperyear) and change ctx.time.iterperyear by ctx.time.nbiter here.
 * @param ctx simulation context (all model state)
 */
void ReadInputClimate(Context &ctx);
/**
 * @brief Read in soil data.
 * @details Global function: This function reads inputs from the soil input file.
 * - in v.3.0, all soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each ctx.diag.layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (ctx.soil.Sat_SWC, ctx.soil.Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
 * - ctx.soil.Sat_SWC and ctx.soil.Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
 * @param ctx simulation context (all model state)
 */
void ReadInputSoil(Context &ctx);
/**
 * @brief Initialisation with bare ground conditions.
 * @param ctx simulation context (all model state)
 */
void Initialise(Context &ctx);
/**
 * @brief Initialise intraspecific variables.
 */
void InitialiseIntraspecific(Context &ctx);
/**
 * @brief Initialise lookup table for LAImax.
 */
void InitialiseLookUpLAImax(Context &ctx);
/**
 * @brief Initialise lookup tables.
 */
void InitialiseLookUpTables(Context &ctx);
/**
 * @brief Initialisation of output streams and intraspecific variables.
 * @param ctx simulation context (all model state)
 */
void InitialiseOutputStreams(Context &ctx);
/**
 * @brief Updated in v.3.1: initialisation from inventories.
 * @param ctx simulation context (all model state)
 */
void ReadInputInventory(Context &ctx);
/**
 * @brief Introduced in v.3.1.6: reads the parameter sheet for the point cloud simulation.
 * @param ctx simulation context (all model state)
 */
void ReadInputPointcloud(Context &ctx);
/**
 * @brief Field dynamic memory allocation.
 * @param ctx simulation context (all model state)
 */
void AllocMem(Context &ctx);
/**
 * @brief Evolution at each timestep.
 */
void Evolution(Context &ctx);
/**
 * @brief Compute field Seed.
 */
void UpdateSeeds(Context &ctx);
/**
 * @brief Compute seed dispersal.
 * This routine implements the reproduction stage (nbs seeds are produced per tree) and dispersal of the trees
 * - Dispersal is equiprobable in all direction and normally distributed (for 2D normal distribution, the absolute distance follows a Rayleigh distribution)
 * - Reproduction only occurs for mature trees
 * - New v.2.1 threshold of maturity is defined as a size threshold (and not age as before), following Wright et al 2005 JTE
 */
void DisperseSeed(Context &ctx, Tree &tree);
/**
 * @brief Update SPECIES_SEEDS field; v.2.5.
 */
void FillSeed(Context &ctx, int col, int row, int spp);
/**
 * @brief Update all fields.
 * This is an important function for TROLL. It includes many of the operations.
 * Sets the iteration environment.
 * @note the current structure of code suppose that environment is periodic (a period = a year),
 * if one wants to input a variable climate, with interannual variation and climate change
 * along the simulation, a full climatic input needs to be input
 * (ie number of columns=ctx.time.iter and not ctx.time.iterperyear)
 * and change ctx.time.iterperyear by ctx.time.nbiter here.
 */
void UpdateField(Context &ctx);
/**
 * @brief Treefall gap formation; v.2.4.
 */
void TriggerTreefall(Context &ctx);
/**
 * @brief Secondary treefall gap formation.
 */
void TriggerTreefallSecondary(Context &ctx);
/**
 * @brief Tree germination module; v.2.5.
 */
void RecruitTree(Context &ctx);
/**
 * @brief Output of the global averages every timestep.
 */
void Average(Context &ctx);
/**
 * @brief Output of the field variables every timestep.
 */
void OutputField(Context &ctx);
/**
 * @brief Output snapshots of the scene at one point in time.
 */
void OutputSnapshot(Context &ctx, fstream &output, bool header, float dbh_limit);
/**
 * @brief Writes the whole 3D LAI voxel field to file.
 */
void OutputLAI(Context &ctx, fstream &output_transmLAI3D);
/**
 * @brief Outputs CHM.
 */
void OutputCHM(Context &ctx, fstream &output_CHM);
/**
 * @brief Output function for visualization purposes.
 */
void OutputVisual(Context &ctx);
#ifdef WATER
/**
 * @brief Write tree variables to file.
 */
void OutputTreeStandard(Context &ctx, Tree &tree, std::ostream &output);
/**
 * @brief Write tree variables to standard output.
 */
void OutputTreeStandard(Context &ctx, Tree &tree);
#endif
/**
 * @brief Compute DBH histograms.
 */
void histdbh(Context &ctx, Tree &tree);
/**
 * @brief Compute average and update species output variables.
 */
void Average(Context &ctx, Tree &tree);
void CloseOutputs(Context &ctx);
/**
 * @brief Free dynamic memory.
 * @param ctx simulation context (all model state)
 */
void FreeMem(Context &ctx);

/**
 * @brief Point cloud output, v.3.1.6; kept separately from other output functions, as we write to a dedicated external file format (.las).
 * @param mean_beam            mean number of laser beams per site
 * @param sd_beam              standard deviation of the number of beams per site
 * @param klaser               extinction coefficient of the laser in foliage
 * @param transmittance_laser  transmittance of the laser through empty voxels
 * @param output_pointcloud    output stream of the point cloud
 */
void ExportPointcloud(float mean_beam, float sd_beam, float klaser, float transmittance_laser, fstream &output_pointcloud);

// HELPER FUNCTIONS
/**
 * @brief Helper function, new in v.3.1: converts current iteration into time of year, also works backwards (negative iterations).
 */
int GetTimeofyear(Context &ctx);
/**
 * @brief Helper function: calculates mean predicted height from allometry.
 */
float CalcHeightBaseline(float &ah, float &hmax, float &dbh);
/**
 * @brief Helper function, new in v.3.1: calculates mean predicted crown radius from allometry.
 */
float CalcCDBaseline(Context &ctx, float &height);
/**
 * @brief Helper function new in v.3.1: calculates mean predicted crown diameter from allometry.
 */
float CalcCRBaseline(Context &ctx, float &dbh);
/**
 * @brief Helper function: returns index for LookUpTables of absorbed flux (considering leaves above and within voxel).
 * @details Upper bound on LAI within one voxel and above voxel (beyond 9.95 and 19.95, none of the environmental variables should change), needed for LookUp tables.
 */
int CalcIntabsorb(float absorb_prev, float absorb_delta);
/**
 * @brief Helper function: returns index for LookUpTables of absorbed flux (only considering leaves above voxel).
 * @details Upper bound on LAI within one voxel and above voxel, needed for LookUp tables, for cases where there won't be any absorption within the voxel (absorb_delta = 0.0).
 */
int CalcIntabsorb(float absorb_prev);

// DATA READING FUNCTIONS
// update v.3.1: already defined here to enable calling from functions of class Tree()
/**
 * @brief See param_registry.hpp: numeric parameter with range check.
 * @param parameter_name parameter name
 * @param parameter_value value read from the file (text)
 * @param parameter [out] variable receiving the value
 * @param parameter_min minimum allowed value
 * @param parameter_max maximum allowed value
 * @param parameter_default default value used when the value is invalid
 * @param quiet suppress console messages if true
 */
template <typename N>
void SetParameter(string &parameter_name, string &parameter_value, N &parameter, N parameter_min, N parameter_max, N parameter_default, bool quiet);
/**
 * @brief See param_registry.hpp: string parameter.
 * @param parameter_name parameter name
 * @param parameter_value value read from the file (text)
 * @param parameter [out] variable receiving the value
 * @param parameter_default default value used when the value is empty
 * @param quiet suppress console messages if true
 */
void SetParameter(string &parameter_name, string &parameter_value, string &parameter, string parameter_default, bool quiet);

#ifdef Output_ABC
// these are functions needed to create ABC output
/**
 * @brief Global ABC function: initialise ABC conditions.
 * @details Determine area to be taken into consideration for ABC and margin not to be taken into consideration (artefacts along border).
 * @param ctx simulation context (all model state)
 */
void InitialiseABC(Context &ctx);
/**
 * @brief Global ABC function: yearly statistics.
 */
void UpdateMovingAveragesABC(Context &ctx);
/**
 * @brief Global ABC function: update DBH function for ABC routines.
 */
void UpdateDBHtrackingABC(Context &ctx);
/**
 * @brief Global function: calculating the TROLL transmittance field from simulated LiDAR.
 */
void UpdateTransmittanceCHM_ABC(int mean_beam, float sd_beam, float klaser, float transmittance_laser);
/**
 * @brief Global function: write headers for ABC outputs.
 */
void OutputABCWriteHeaders(fstream &output_traitconservation, fstream &output_field, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_transmittance, fstream &output_transmittance_ALS, fstream &output_LAIfield, fstream &output_LAIfield_ALS, fstream &output_chmpotential, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass);
/**
 * @brief Global ABC function: assesses the conservation of traits between input and output in TROLL.
 */
void OutputABCConservationTraits(Context &ctx, fstream &output_traitconservation);
/**
 * @brief Global ABC function: returns ABC ground outputs.
 */
void OutputABC_ground(Context &ctx, fstream &output_field);
/**
 * @brief Global ABC function: returns ABC outputs for species.
 */
void OutputABC_species(Context &ctx, fstream &output_species, fstream &output_species10, fstream &output_traits, fstream &output_traits10, fstream &output_biomass);
/**
 * @brief Global function: returns ABC outputs for canopy height model (CHM).
 */
void OutputABC_CHM(Context &ctx, fstream &output_CHM, fstream &output_CHM_ALS, fstream &output_chmpotential);
/**
 * @brief Global ABC function: ABC outputs.
 */
void OutputABC_transmittance(Context &ctx, fstream &output_transmittance, fstream &output_transmittance_ALS);
/**
 * @brief Global ABC function: output general ABC statistics.
 */
void OutputABC(Context &ctx);
/**
 * @brief Global ABC function: creates a Canopy Height Model and LAD profile.
 */
void OutputCHM(Context &ctx, fstream &output_CHM);
#endif

#ifdef TRACK_INDIVIDUALS
/**
 * @brief Global function: tree level tracking of key variables.
 */
void TrackingData_andOutput(Context &ctx);
/**
 * @brief Global function: initialize tree tracking variables.
 */
float StartTracking(Context &ctx, Tree &tree);
#endif

// Species and Tree class declarations moved to include/species.hpp and include/tree.hpp
// (included transitively via context.hpp)

//// UTIL Templates

/**
 * @brief True if value1 and value2 are equal to the given number of decimal places.
 * @param value1 first value
 * @param value2 second value
 * @param precision number of decimal places compared
 */
template <typename N>
bool AreEqual(N value1, N value2, int precision)
{
    return abs(value1 - value2) < pow(10, -precision); // make sure that namespace std is activated, otherwise abs is only defined for integers
}

/**
 * @brief Value of parameter_name in the parallel lists of names and values read from a parameter file ("" if absent).
 * @param parameter_name parameter name
 * @param parameter_names parameter names read from a parameter file
 * @param parameter_values values read, in the same order as parameter_names
 */
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

/**
 * @brief Report how a parameter value compares with its allowed range and default (diagnostics when reading inputs).
 * @param parameter_value value read from the file (text)
 * @param parameter_min minimum allowed value
 * @param parameter_max maximum allowed value
 * @param parameter_default default value used when the value is invalid
 * @param parameter_name parameter name
 * @param parameter_names parameter names read from a parameter file
 * @param parameter_values values read, in the same order as parameter_names
 * @param quiet suppress console messages if true
 */
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
