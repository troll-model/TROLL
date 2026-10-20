#ifndef TROLL_TREE_HPP
#define TROLL_TREE_HPP

#include <vector>
#include <fstream>
#include <string>
#include "troll_defines.hpp"

// Forward declaration — full definition in context.hpp
struct Context;

#ifdef WATER
//! Global structure: For GPPleaf FUNCTION
struct leafFluxes
{
    float carbon_flux; //!< Net leaf carbon assimilation (umol CO2 m-2 s-1)
    float water_flux;  //!< Leaf water flux: transpiration from FluxesLeaf (umol H2O m-2 s-1); stomatal conductance to CO2 from Photosyn (mol m-2 s-1)
};
#endif

//! Daily carbon and water fluxes of one tree, as computed by Tree::ComputeRespGPP().
//! Computing them reads the tree and shared fields but writes nothing, so it can run
//! for all trees in parallel; Tree::ApplyRespGPP() then stores them in tree order.
struct RespGPPResult
{
    float GPP = 0.0f;                //!< Gross primary production of the tree in the timestep (gC, as Tree::t_GPP)
    float transpiration = 0.0f;      //!< Transpiration of the tree in the timestep (as Tree::t_transpiration)
    float transpiration_1016 = 0.0f; //!< contribution to ctx.soil.transpiration_1016
    float Rstem = 0.0f;              //!< Stem respiration of the tree in the timestep (gC)
    float Rnight = 0.0f;             //!< Leaf respiration at night in the timestep (gC)
};

// ###########################################
//! Tree class
// ############################################

class Tree
{
public:
    int t_site;              //!< Geolocation of the tree
    int t_NPPneg;            //!< Diagnostic variable: number of consecutive timesteps with NPP<0; v.2.2
    int t_CrownDisplacement; //!< Displacement of the crown center with respect to the stem.
#ifdef MIP_Lichstein
    int t_inInventory; //!< 1 if the tree was read from the input inventory (MIP_Lichstein outputs)
#endif

    float t_age;       //!< Tree age, also indicates whether tree is alive (live trees are such that t_age > 0.0)
    float t_hmax;      //!< Allometric parameter, not real maximum
    float t_ah;        //!< Allometric parameter, for consistency with t_hmax also an individual parameter; v.2.4
    float t_dbh;       //!< Diameter at breast height (m) beware: this scales with NH, the horizontal size of voxels
    float t_dbhmature; //!< Reproductive size threshold; v.2.3
    float t_dbhmax;    //!< Maximum diameter at breast height (dbh), as estimated from field data
    float t_height;    //!< Total tree height (m) beware: this scales with NV, the vertical size of voxels, renamed v.3.1 for convenience
    float t_CD;        //!< crown depth (m) beware: this scales with NV, the vertical size of voxels, renamed v.3.1 for convenience
    float t_CR;        //!< crown radius (m) beware: this scales with NH, the horizontal size of voxels, renamed v.3.1 for convenience
    float t_Ct;        //!< flexural force threshold, _BASICTREEFALL
    float t_GPP;       //!< Gross primary productivity of the tree (gC/timestep)
    float t_NPP;       //!< Net primary productivity of the tree (gC/timestep)
    float t_Rday;      //!< Daytime leaf respiration of the tree (gC/timestep)
    float t_Rnight;    //!< Nighttime leaf respiration of the tree (gC/timestep)
    float t_Rstem;     //!< Stem respiration (gC/timestep)
    float t_LA;        //!< Total crown leaf area (m^2); v.2.2, renamed in v.3.1 for convenience
    float t_youngLA;   //!< Total young leaf area (m^2); v.2.2
    float t_matureLA;  //!< Total mature leaf area (m^2); v.2.2
    float t_oldLA;     //!< Total old leaf area (m^2); v.2.2
    float t_LAI;       //!< Total leaf area index (m^2/m^2); v.2.5
    float t_litter;    //!< Tree litterfall at each timestep, in dry mass (g); v.2.2

    std::vector<float> t_NDDfield; //!< Tree field useful when option _NDD is activated

#ifdef WATER
    float t_root_depth;                //!< Tree rooting depth (m)
    float t_phi_root;                  //!< Soil water potential in the root zone (MPa)
    std::vector<float> t_root_biomass;      //!< Tree root biomass in each soil layer, in g
    std::vector<float> t_soil_layer_weight; //!< Soil layer weight (to compute t_phi_root), if different from root biomass
    float t_WSF;                       //!< Tree water stress factor for stomatal conductance, unitless, between 0 and 1
    float t_WSF_A;                     //!< Tree water stress factor for photosynthetic capacities, unitless, between 0 and 1
    float t_transpiration;             //!< Amount of water taken up from the soil and transpired at each timestep
    float t_g1_0;                      //!< Medlyn stomatal slope g1 of the tree from wood specific gravity, g1 = 6.53 - 3.97 wsg (Lin et al. 2015) (kPa^0.5)
    float t_g1;                        //!< Medlyn g1 used for stomatal conductance (kPa^0.5): t_g1_0 x WSF_s (stomatal water stress factor, set in Water_availability; GMD Part 1 Eq. 36)

#ifdef PHENO_DROUGHT
    int t_Ndays_dry;       //!< Days with root-zone water potential below the species turgor loss point (PHENO_DROUGHT)
    int t_Ndays_wet;       //!< Counterpart of t_Ndays_dry for days above the turgor loss point (PHENO_DROUGHT)
    float t_pheno_factor;  //!< Multiplicative factor (<= 1) on the residence time of old leaves under drought; f0 in GMD Part 2, Sect. 2
#endif

#endif

    // new in v.2.4.0: traits defined at individual level
    float t_Pmass; //!< Leaf phosphorus content, defined at tree scale (g/g)
    float t_Nmass; //!< Leaf nitrogen content, defined at tree scale (g/g))
    float t_LMA;   //!< Leaf mass per area (LMA), defined at tree scale (g/m^2)
#ifdef WATER
    float t_leafarea;   //!< Leaf area (LA), defined at tree scale
    float t_tlp;        //!< Leaf turgor loss point (TLP), defined at tree scale
    float t_wleaf;      //!< Leaf width (in m)
    float t_phi_lethal; //!< plant water potential at drought-induced death (MPa)
    float t_itlp;       //!< inverse of t_tlp
    float t_b;          //!< parameter used to compute the water stress factor (WSF) for stomatal limitation
#endif
    float t_wsg;   //!< Wood specific gravity, defined at tree scale (g/cm^3)
    float t_Rdark; //!< Dark respiration rate at PPFD = 0 (micromol C/m^2/s)
    float t_Vcmax; //!< Maximal rate of carboxylation, on an area basis (micromolC/m^2/s)
    float t_Jmax;  //!< Maximal rate of electron transport, on an area basis (micromol/m^2/s)
    float t_leaflifespan;  //!< Average leaf lifespan (months)
    float t_lambda_young;  //!< Residence time in the young leaf class; v.2.5
    float t_lambda_mature; //!< Residence time in the mature leaf class; v.2.5
    float t_lambda_old;    //!< Residence time in the old leaf class; v.2.5

    float t_mult_height; //!< Intraspecific multiplier for height; v.2.4.0
    float t_mult_CR;     //!< Intraspecific multiplier for crown radius; v.2.4.0
    float t_mult_CD;     //!< Intraspecific multiplier for crown depth; v.2.4.0
    float t_mult_P;      //!< Intraspecific multiplier for leaf phosphorus concentration; v.2.4.0
    float t_mult_N;      //!< Intraspecific multiplier for leaf nitrogen concentration; v.2.4.0
    float t_mult_LMA;    //!< Intraspecific multiplier for leaf mass per area; v.2.4.0
    float t_mult_dbhmax; //!< Intraspecific multiplier for dbh maximum; v.2.4.0
    float t_dev_wsg;     //!< Intraspecific absolution deviation for wood specific gravity; v.2.4.0
#ifdef WATER
    float t_mult_leafarea; //!< Intraspecific multiplier for leaf area
    float t_mult_tlp;      //!< Intraspecific multiplier for turgor loss point
#endif

    float t_LAImax;          //!< Maximal LAI
    float t_LAmax;           //!< Maximal leaf area
    float t_carbon_storage;  //!< Persistent C storage pool.
    float t_carbon_biometry; //!< Temporary biometry C storage pool.
    float t_fraction_filled; //!< Filled fraction of tree crown (1-gap fraction); v.2.5.0

    int t_multiplier_seed; //!< Not documented  !!!UPDATE
    float t_sapwood_area;  //!< Sapwood area !!!UPDATE units, description

#ifdef CROWN_UMBRELLA
    float t_Crown_Slope; //!< Not documented  !!!UPDATE
#endif

    unsigned short
        t_from_Data, //!< Indicator: tree was born through initialisation (1) or through simulation routine (0)
        t_sp_lab,    //!< Stores the species label.
        t_hurt;      //!< Stories the treefall index

#ifdef Output_ABC
    float t_dbh_previous; //!< Not documented  !!!UPDATE
#endif

#ifdef TRACK_INDIVIDUALS
    int t_month_born; //!< Month of birth (TRACK_INDIVIDUALS; note: tree.cpp uses t_timeofyear_born, so this option does not compile)

    int t_seedsproduced; //!< Seeds produced in the timestep (TRACK_INDIVIDUALS)
    int t_seedsproduced_sumyear; //!< Seeds produced in the current year (TRACK_INDIVIDUALS)
    int t_time_carbonstarvation; //!< Number of timesteps with carbon starvation (TRACK_INDIVIDUALS)
    int t_time_carbonstarvation_year; //!< Number of timesteps with carbon starvation in the current year (TRACK_INDIVIDUALS)
    float t_GPP_sumyear; //!< Sum of GPP over the current year (TRACK_INDIVIDUALS)
    float t_NPP_sumyear; //!< Sum of NPP over the current year (TRACK_INDIVIDUALS)
    float t_GPPsquared_sumyear; //!< Sum of squared GPP over the current year (for the variance) (TRACK_INDIVIDUALS)
    float t_NPPsquared_sumyear; //!< Sum of squared NPP over the current year (for the variance) (TRACK_INDIVIDUALS)
    float t_Rday_sumyear; //!< Sum of daytime leaf respiration over the current year (TRACK_INDIVIDUALS)
    float t_Rnight_sumyear; //!< Sum of night-time leaf respiration over the current year (TRACK_INDIVIDUALS)
    float t_Rstem_sumyear; //!< Sum of stem respiration over the current year (TRACK_INDIVIDUALS)
    float t_LAIabove_effavgyear; //!< Mean effective LAI above the tree over the current year (TRACK_INDIVIDUALS)
    float t_carbon_storage_avgyear; //!< Mean carbon storage over the current year (TRACK_INDIVIDUALS)
    float t_dbh_tracked; //!< DBH at the start of the tracking period (m) (TRACK_INDIVIDUALS)
    float t_height_tracked; //!< Height at the start of the tracking period (m) (TRACK_INDIVIDUALS)
    float t_CR_tracked; //!< Crown radius at the start of the tracking period (m) (TRACK_INDIVIDUALS)
    float t_agb_tracked; //!< Aboveground biomass at the start of the tracking period (TRACK_INDIVIDUALS)
    float t_LAIcum; //!< LAI above the tree, cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_LAIeffcum; //!< Effective LAI above the tree, cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_GPPcum; //!< GPP cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_NPPcum; //!< NPP cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_LAIsquared_cum; //!< Squared LAI cumulated over the tracking period (for the variance) (TRACK_INDIVIDUALS)
    float t_LAIeffsquared_cum; //!< Squared effective LAI cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_GPPsquared_cum; //!< Squared GPP cumulated over the tracking period (TRACK_INDIVIDUALS)
    float t_NPPsquared_cum; //!< Squared NPP cumulated over the tracking period (TRACK_INDIVIDUALS)
#endif
    // Constructor
    //! Constructor: an empty site (no tree), with per-tree vectors sized for this simulation
    //! @param ctx simulation context (all model state)
    Tree(Context &ctx);
    void Birth(Context &ctx, int, int);
    int BirthFromInventory(Context &ctx, int site, std::vector<std::string> &parameter_names, std::vector<std::string> &parameter_values, int &nb_speciesrandom);
    void Death(Context &ctx);
    void Growth(Context &ctx, const RespGPPResult *precomputed = nullptr);
#ifdef WATER
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer, float &PPFD_incident, float &ExtinctLW);
#else
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer);
#endif

#ifdef WATER
    void Water_availability(Context &ctx);
    void Water_uptake(Context &ctx);
#endif

    void CalcRespGPP(Context &ctx); //!< Carbon and water fluxes of the tree in the timestep: GPP, transpiration, stem and night-time respiration (WATER: ComputeRespGPP then ApplyRespGPP)
#ifdef WATER
    RespGPPResult ComputeRespGPP(Context &ctx); //!< Fluxes for this timestep; modifies neither the tree nor ctx
    void ApplyRespGPP(Context &ctx, const RespGPPResult &r); //!< Store fluxes in the tree; add to ctx.soil.transpiration_1016
#endif
    void CalcNPP(Context &ctx); //!< Net primary production of the timestep, from GPP and respiration
    void UpdateLeafDynamics(Context &ctx);
    void UpdateTreeBiometry(Context &ctx);
    void UpdateVolumeDensity(); //!< Crown leaf area index t_LAI: leaf area over the crown area filled with leaves (gaps excluded)
    void DisperseSeed(Context &ctx);
    void Treefall(Context &ctx, float angle);
    void Update(Context &ctx, bool fluxes_precomputed = false);
    void Average(Context &ctx); //!< Add the tree to the per-species sums of its species (sumstats outputs)
    void CalcLAI(Context &ctx);
    void histdbh(Context &ctx); //!< Add the tree to the DBH histogram ctx.diag.nbdbh (1 cm classes, rounded down)
    float CalcJmaxm(Context &ctx); //!< Maximum electron transport rate per leaf mass from the tree's leaf traits
    float CalcAGB(Context &ctx); //!< Aboveground biomass of the tree (kg)
    float CalcIncrementDBH(Context &ctx, float delta_agb); //!< DBH increment corresponding to an aboveground biomass increment delta_agb
    float CalcCarbonStorageMax(Context &ctx); //!< Maximum non-structural carbon storage: 5% of the tree carbon (half of NSC), x1.25 for roots (gC)
    float CalcCt(Context &ctx);

    void CalcLeafLifespan(Context &ctx);
    void InitialiseLeafPools(Context &ctx); //!< Initial partition of leaf area into young, mature and old leaves
    void UpdateSapwoodArea(Context &ctx, float ddbh);
    void UpdateHeight(Context &ctx); //!< Height from DBH (Michaelis-Menten allometry, GMD Part 1 Eq. 16)
    void UpdateCR(Context &ctx); //!< Crown radius from DBH (GMD Part 1 Eq. 17)
    void UpdateCD(Context &ctx); //!< Crown depth from height (GMD Part 1 Eq. 18)
    float DeathRateNDD(Context &ctx, float dbh, float carbon_starv, float ndd); //!< Death rate including negative density dependence (NDD option)
#ifdef WATER
    float DeathRate(Context &ctx, float dbh, float carbon_starv, float phi_root); //!< Death rate of the tree per timestep (WATER version)
    //! Leaf energy balance coupled with photosynthesis and stomatal conductance, for one leaf
    //! layer and one half-hourly step: iterates leaf temperature to convergence (GMD Part 1, Sect. 2.5).
    //! @param ctx        simulation context (parameters, lookup tables)
    //! @param PPFD       photosynthetic photon flux density on the leaves of the layer (umol photons m-2 s-1)
    //! @param VPDa       air vapour pressure deficit at the layer (kPa)
    //! @param Ta         air temperature at the layer (degC)
    //! @param WIND       wind speed at the layer (m s-1)
    //! @param ExtinctLW  extinction factor of net longwave radiation at the layer (LookUp_ExtinctLW, dimensionless)
    //! @param PPFDtop    PPFD at the top of the canopy (umol photons m-2 s-1)
    //! @param PPFDinc    PPFD incident on the layer, i.e. above it (umol photons m-2 s-1)
    //! @param Tatop      air temperature at the top of the canopy (degC), for the net longwave lookup
    //! @param VPDatop    vapour pressure deficit at the top of the canopy (kPa), for the net longwave lookup
    //! @return carbon_flux: net assimilation (umol CO2 m-2 s-1); water_flux: transpiration (umol H2O m-2 s-1)
    leafFluxes FluxesLeaf(Context &ctx, float PPFD, float VPDa, float Ta, float WIND, float ExtinctLW, float PPFDtop, float PPFDinc, float Tatop, float VPDatop);
    //! Daily mean leaf fluxes of one crown layer: FluxesLeaf over the daytime half-hourly steps.
    //! The climate arguments are relative to the canopy-top half-hourly series, not absolute values.
    //! @param ctx        simulation context (climate series)
    //! @param PPFD       light at the layer relative to the canopy top (fraction; multiplied by the half-hourly irradiance and SWtoPPFD)
    //! @param VPD        VPD at the layer relative to the canopy top (factor)
    //! @param T          decrease in air temperature from the canopy top to the layer (degC)
    //! @param W          wind speed at the layer relative to the canopy top (factor)
    //! @param ExtinctLW  extinction factor of net longwave radiation at the layer (dimensionless)
    //! @param PPFDinc    light incident on the layer relative to the canopy top (fraction)
    //! @param wf_1016    [out] transpiration of the layer for the 10-16 cm DBH diagnostic (transpiration1016 output)
    //! @return carbon_flux and water_flux averaged over the daytime steps (units as FluxesLeaf)
    leafFluxes dailyFluxesLeaf(Context &ctx, float PPFD, float VPD, float T, float W, float ExtinctLW, float PPFDinc, float &wf_1016);
    leafFluxes Photosyn(Context &ctx, float PPFD, float TLEAF, float CS, float DS); //!< Leaf photosynthesis and stomatal conductance (see the definition for parameters)
    void OutputTreeStandard(Context &ctx, std::fstream &output); //!< Write the tree's state to a snapshot or site output stream
    void OutputTreeStandard(Context &ctx); //!< Write the tree's state to standard output
#else
    float DeathRate(Context &ctx, float dbh, float carbon_starv); //!< Death rate of the tree per timestep (version without water)
    float GPPleaf(Context &ctx, float, float, float);
    float dailyGPPleaf(Context &ctx, float, float, float);
    float dailyGPPcrown(Context &ctx, float, float, float, float);
#endif
    float Rdayleaf(Context &ctx, float T);
    float dailyRdayleaf(Context &ctx, float T); //!< Daytime leaf respiration averaged over the half-hourly steps of a day
    void CalcLAImax(Context &ctx);
    void CalcLAmax(Context &ctx, float &LAIexperienced_eff, float &LAmax);
    float predLeafLifespanKikuzawa(Context &ctx);
    float GetCrownAreaFilled(float crown_area); //!< Crown area actually filled with leaves, given the crown gap fraction (m2)
#ifdef TRACK_INDIVIDUALS
    float StartTracking(Context &ctx);
#endif
};

#endif
