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
    /**
     * @brief Constructor: an empty site (no tree).
     *
     * Per-tree vectors are sized for this simulation.
     *
     * @param ctx Simulation context (global model state).
     */
    Tree(Context &ctx);

    /**
     * @brief Actions needed at birth for a tree.
     *
     * More precisely, the actions performed when a stem enters the > 1 cm
     * trunk diameter class.
     *
     * @param ctx   Simulation context (global model state).
     * @param nume  Species index of the recruit.
     * @param site0 Site where the tree is recruited.
     */
    void Birth(Context &ctx, int nume, int site0);

    /**
     * @brief Initialize a tree from field data (input inventory).
     *
     * Modelled for compatibility with Tree::Birth; for comments regarding
     * allometries and t_LA, see Tree::Growth.
     *
     * @param ctx              Simulation context (global model state).
     * @param site             Site of the tree.
     * @param parameter_names  Column names of the inventory file.
     * @param parameter_values Values of this tree's row.
     * @param nb_speciesrandom [in,out] Number of trees whose species was not
     *                                  found and was drawn at random.
     * @return 1 on success, 0 otherwise (tree not read from the inventory).
     */
    int BirthFromInventory(Context &ctx, int site, std::vector<std::string> &parameter_names, std::vector<std::string> &parameter_values, int &nb_speciesrandom);

    /**
     * @brief Record the death of a tree and empty the site.
     *
     * Records the basic properties of a tree that has died and empties the
     * variables at the tree site. It does not apply the destructor to the
     * tree object.
     *
     * @param ctx Simulation context (global model state).
     */
    void Death(Context &ctx);

    /**
     * @brief Tree growth during one timestep.
     *
     * Updates age, computes GPP and respiration, then allocates carbon to
     * leaves, storage and growth (dbh, height, crown size).
     *
     * @param ctx Simulation context (global model state).
     */
    void Growth(Context &ctx);

#ifdef WATER
    /**
     * @brief Mean light flux, VPD and temperature received by one crown layer.
     *
     * Computes the light (PPFD), VPD and air temperature received by a tree
     * crown layer at height h. The values are weighted by the leaf area in the
     * layer. The climate quantities are now local variables, not tree
     * variables.
     *
     * @param ctx            Simulation context (global model state).
     * @param h              Height layer (m).
     * @param PPFD           [out] Light at the layer relative to the canopy top (fraction).
     * @param VPD            [out] VPD at the layer relative to the canopy top (factor).
     * @param Tmp            [out] Decrease in air temperature from the canopy top (degC).
     * @param leafarea_layer [out] Leaf area of the crown in the layer (m2).
     * @param PPFD_incident  [out] Light incident on the layer relative to the canopy top (fraction).
     * @param ExtinctLW      [out] Extinction factor of net longwave radiation (dimensionless).
     */
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer, float &PPFD_incident, float &ExtinctLW);
#else
    /**
     * @brief Mean light flux, VPD and temperature received by one crown layer.
     *
     * Computes the light (PPFD), VPD and air temperature received by a tree
     * crown layer at height h. The values are weighted by the leaf area in the
     * layer.
     *
     * @param ctx            Simulation context (global model state).
     * @param h              Height layer (m).
     * @param PPFD           [out] Light at the layer relative to the canopy top (fraction).
     * @param VPD            [out] VPD at the layer relative to the canopy top (factor).
     * @param Tmp            [out] Decrease in air temperature from the canopy top (degC).
     * @param leafarea_layer [out] Leaf area of the crown in the layer (m2).
     */
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer);
#endif

#ifdef WATER
    /**
     * @brief Compute the tree's water availability in the root zone.
     *
     * Updates the root-zone soil water potential t_phi_root and the water
     * stress factors t_WSF / t_WSF_A from the soil water conditions.
     *
     * @param ctx Simulation context (global model state).
     */
    void Water_availability(Context &ctx);

    /**
     * @brief Add the tree's contribution to the stand transpiration field.
     *
     * Adds up each tree's contribution to the soil transpiration field, i.e.
     * the water removed from the soil through the tree's transpiration.
     * Similar to CalcLAI, which adds each tree's contribution to the LAI3D
     * field. Called by UpdateField.
     *
     * @param ctx Simulation context (global model state).
     */
    void Water_uptake(Context &ctx);
#endif

    /**
     * @brief Calculate the tree's GPP and carbon/water fluxes.
     *
     * Carbon and water fluxes of the tree in the timestep: GPP, transpiration,
     * stem and night-time respiration. Since v.2.4 t_GPP and t_Rday are
     * updated in a separate function. Under WATER, Rday is computed within
     * Photosyn/FluxesLeaf and is included in the GPP.
     *
     * @param ctx Simulation context (global model state).
     */
    void CalcRespGPP(Context &ctx);

    /**
     * @brief Calculate the tree's net primary production.
     *
     * Net primary production of the timestep, from GPP and respiration.
     *
     * @param ctx Simulation context (global model state).
     */
    void CalcNPP(Context &ctx);

    /**
     * @brief Update leaf dynamics and the carbon allocation to leaves.
     *
     * Allocates NPP to leaves within the three leaf-age classes (young,
     * mature, old); only the old leaves generate litterfall, and the dynamics
     * of the leaf cycle is generated by the dynamics of NPP.
     *
     * @param ctx Simulation context (global model state).
     */
    void UpdateLeafDynamics(Context &ctx);

    /**
     * @brief Update tree biometry, including allometry.
     *
     * Converts the carbon allocated to growth into dbh, height and crown
     * increments. Standalone function since v.2.3.0.
     *
     * @param ctx Simulation context (global model state).
     */
    void UpdateTreeBiometry(Context &ctx);

    /**
     * @brief Update the crown leaf area density.
     *
     * Crown leaf area index t_LAI: leaf area over the crown area filled with
     * leaves (gaps excluded).
     */
    void UpdateVolumeDensity();

    /**
     * @brief Create a treefall.
     *
     * Tree falling routine (formerly FallTree, _BASICTREEFALL). Creates a
     * treefall (but no longer a treefall probability). Takes the fall angle as
     * an argument and can be used for primary, secondary treefalls, forestry
     * or other disturbances.
     *
     * @param ctx   Simulation context (global model state).
     * @param angle Direction of the fall (radians).
     */
    void Treefall(Context &ctx, float angle);

    /**
     * @brief Update the tree for one timestep (growth or death).
     *
     * Calls the appropriate DeathRate module; if the death condition is met,
     * Tree::Death() is called, otherwise Tree::Growth().
     *
     * @param ctx Simulation context (global model state).
     */
    void Update(Context &ctx);

    /**
     * @brief Maximum electron transport rate per leaf mass.
     *
     * From the tree's leaf traits.
     *
     * @param ctx Simulation context (global model state).
     * @return Maximum electron transport rate per leaf mass.
     */
    float CalcJmaxm(Context &ctx);

    /**
     * @brief Aboveground biomass of the tree.
     *
     * @param ctx Simulation context (global model state).
     * @return Aboveground biomass of the tree (kg).
     */
    float CalcAGB(Context &ctx);

    /**
     * @brief DBH increment corresponding to an aboveground biomass increment.
     *
     * @param ctx       Simulation context (global model state).
     * @param delta_agb Aboveground biomass increment (kg).
     * @return DBH increment (m).
     */
    float CalcIncrementDBH(Context &ctx, float delta_agb);

    /**
     * @brief Maximum non-structural carbon storage.
     *
     * 5% of the tree carbon (half of NSC), x1.25 for roots.
     *
     * @param ctx Simulation context (global model state).
     * @return Maximum carbon storage (gC).
     */
    float CalcCarbonStorageMax(Context &ctx);

    /**
     * @brief Compute the treefall threshold.
     *
     * Used when the _BASICTREEFALL option is activated. Slightly updated in
     * v.3.1.
     *
     * @param ctx Simulation context (global model state).
     * @return Flexural force threshold.
     */
    float CalcCt(Context &ctx);

    /**
     * @brief Determine the leaf lifespan.
     *
     * Either from an empirical function of leaf traits or from the Kikuzawa
     * model.
     *
     * @param ctx Simulation context (global model state).
     */
    void CalcLeafLifespan(Context &ctx);

    /**
     * @brief Initial partition of leaf area into the age classes.
     *
     * Partitions the leaf area into young, mature and old leaves.
     *
     * @param ctx Simulation context (global model state).
     */
    void InitialiseLeafPools(Context &ctx);

    /**
     * @brief Update the sapwood area, limited by the dbh increase.
     *
     * Options: fixed percentage (ctx.opt._sapwood == 0) or based on the tree's
     * leaf area (ctx.opt._sapwood > 0).
     *
     * @param ctx  Simulation context (global model state).
     * @param ddbh DBH increase (m).
     */
    void UpdateSapwoodArea(Context &ctx, float ddbh);

    /**
     * @brief Update the tree height.
     *
     * Height from DBH (Michaelis-Menten allometry, GMD Part 1 Eq. 16).
     *
     * @param ctx Simulation context (global model state).
     */
    void UpdateHeight(Context &ctx);

    /**
     * @brief Update the crown radius.
     *
     * Crown radius from DBH (GMD Part 1 Eq. 17).
     *
     * @param ctx Simulation context (global model state).
     */
    void UpdateCR(Context &ctx);

    /**
     * @brief Update the crown depth.
     *
     * Crown depth from height (GMD Part 1 Eq. 18).
     *
     * @param ctx Simulation context (global model state).
     */
    void UpdateCD(Context &ctx);

    /**
     * @brief Death rate including negative density dependence (NDD).
     *
     * The function returns a dimensionless probability that the tree dies
     * during the current simulation step.  It combines a baseline mortality
     * based on the tree diameter with a correction factor representing
     * competition from neighbouring trees (NDD).
     *
     * @param ctx           Simulation context (global model state).
     * @param dbh           Diameter at breast height (meters).
     * @param carbon_starv  Number of timesteps the tree has suffered carbon starvation.
     * @param ndd           Negative density-dependence index of the tree's neighbourhood.
     *
     * @return The death rate (unit-less, typically a probability per step).
     */
    float DeathRateNDD(Context &ctx,
                       float dbh,
                       float carbon_starv,
                       float ndd);
#ifdef WATER
    /**
     * @brief Death rate of the tree per timestep (WATER version).
     *
     * @param ctx          Simulation context (global model state).
     * @param dbh          Diameter at breast height (m).
     * @param carbon_starv Number of timesteps with carbon starvation.
     * @param phi_root     Water potential in the root zone (MPa).
     * @return The death rate of the tree per timestep.
     */
    float DeathRate(Context &ctx, float dbh, float carbon_starv, float phi_root);

    /**
     * @brief Leaf energy balance coupled with photosynthesis and stomatal conductance.
     *
     * For one leaf layer and one half-hourly step: iterates leaf temperature
     * to convergence (GMD Part 1, Sect. 2.5).
     *
     * @param ctx       Simulation context (parameters, lookup tables).
     * @param PPFD      Photosynthetic photon flux density on the leaves of the layer (umol photons m-2 s-1).
     * @param VPDa      Air vapour pressure deficit at the layer (kPa).
     * @param Ta        Air temperature at the layer (degC).
     * @param WIND      Wind speed at the layer (m s-1).
     * @param ExtinctLW Extinction factor of net longwave radiation at the layer (LookUp_ExtinctLW, dimensionless).
     * @param PPFDtop   PPFD at the top of the canopy (umol photons m-2 s-1).
     * @param PPFDinc   PPFD incident on the layer, i.e. above it (umol photons m-2 s-1).
     * @param Tatop     Air temperature at the top of the canopy (degC), for the net longwave lookup.
     * @param VPDatop   Vapour pressure deficit at the top of the canopy (kPa), for the net longwave lookup.
     * @return carbon_flux: net assimilation (umol CO2 m-2 s-1); water_flux:
     *         transpiration (umol H2O m-2 s-1).
     */
    leafFluxes FluxesLeaf(Context &ctx, float PPFD, float VPDa, float Ta, float WIND, float ExtinctLW, float PPFDtop, float PPFDinc, float Tatop, float VPDatop);

    /**
     * @brief Daily mean leaf fluxes of one crown layer.
     *
     * FluxesLeaf over the daytime half-hourly steps. The climate arguments are
     * relative to the canopy-top half-hourly series, not absolute values.
     *
     * @param ctx       Simulation context (climate series).
     * @param PPFD      Light at the layer relative to the canopy top (fraction; multiplied by the half-hourly irradiance and SWtoPPFD).
     * @param VPD       VPD at the layer relative to the canopy top (factor).
     * @param T         Decrease in air temperature from the canopy top to the layer (degC).
     * @param W         Wind speed at the layer relative to the canopy top (factor).
     * @param ExtinctLW Extinction factor of net longwave radiation at the layer (dimensionless).
     * @param PPFDinc   Light incident on the layer relative to the canopy top (fraction).
     * @param wf_1016   [out] Transpiration of the layer for the 10-16 cm DBH diagnostic (transpiration1016 output).
     * @return carbon_flux and water_flux averaged over the daytime steps (units as FluxesLeaf).
     */
    leafFluxes dailyFluxesLeaf(Context &ctx, float PPFD, float VPD, float T, float W, float ExtinctLW, float PPFDinc, float &wf_1016);

    /**
     * @brief Leaf photosynthesis and stomatal conductance.
     *
     * Farquhar-von Caemmerer-Berry photosynthesis coupled with Medlyn stomatal
     * conductance at the leaf surface (GMD Part 1, Sect. 2.5.1, Eqs. 24-25).
     *
     * @param ctx   Simulation context.
     * @param PPFD  Photosynthetic photon flux density on the leaf (umol photons m-2 s-1).
     * @param TLEAF Leaf temperature (degC).
     * @param CS    CO2 concentration at the leaf surface (umol mol-1).
     * @param DS    Vapour pressure deficit at the leaf surface (kPa).
     * @return carbon_flux: net assimilation (umol CO2 m-2 s-1); water_flux:
     *         stomatal conductance to CO2 (mol m-2 s-1).
     */
    leafFluxes Photosyn(Context &ctx, float PPFD, float TLEAF, float CS, float DS);
#else
    /**
     * @brief Death rate of the tree per timestep (version without water).
     *
     * @param ctx          Simulation context (global model state).
     * @param dbh          Diameter at breast height (m).
     * @param carbon_starv Number of timesteps with carbon starvation.
     * @return The death rate of the tree per timestep.
     */
    float DeathRate(Context &ctx, float dbh, float carbon_starv);

    /**
     * @brief Leaf gross primary productivity.
     *
     * Primary productivity (net assimilation, i.e. after leaf respiration is
     * at least accounted through Autotrophic respiration) per unit leaf area,
     * from the Farquhar-von Caemmerer-Berry model and the Medlyn et al. 2011
     * stomatal conductance model.
     *
     * @param ctx  Simulation context (global model state).
     * @param PPFD Light on the leaf (non-WATER version).
     * @param VPD  Vapour pressure deficit (kPa).
     * @param T    Air temperature (degC).
     * @return Assimilation per unit leaf area (micromol CO2 m-2 s-1).
     */
    float GPPleaf(Context &ctx, float PPFD, float VPD, float T);

    /**
     * @brief Daily leaf GPP averaged over the half-hourly steps.
     *
     * Daily primary productivity (assimilation; computed from Tree::GPPleaf)
     * averaged across the daily fluctuations in climatic conditions, per unit
     * leaf area.
     *
     * @param ctx  Simulation context (global model state).
     * @param PPFD Light on the leaf, relative to the canopy top (non-WATER version).
     * @param VPD  VPD relative to the canopy top.
     * @param T    Temperature decrease from the canopy top (degC).
     * @return Daily assimilation per unit leaf area (micromol CO2 m-2 s-1).
     */
    float dailyGPPleaf(Context &ctx, float PPFD, float VPD, float T);

    /**
     * @brief Whole-crown daily GPP.
     *
     * Faster, whole-crown GPP calculation, similar to Tree::dailyGPPleaf.
     *
     * @param ctx  Simulation context (global model state).
     * @param PPFD Light on the crown, relative to the canopy top (non-WATER version).
     * @param VPD  VPD relative to the canopy top.
     * @param T    Temperature decrease from the canopy top (degC).
     * @param LAI  Leaf area index of the crown (m2 m-2).
     * @return The whole-crown daily GPP.
     */
    float dailyGPPcrown(Context &ctx, float PPFD, float VPD, float T, float LAI);
#endif

    /**
     * @brief Leaf daytime respiration.
     *
     * @param ctx Simulation context (global model state).
     * @param T  Leaf temperature (degC).
     * @return Daytime leaf respiration (micromol CO2 m-2 s-1).
     */
    float Rdayleaf(Context &ctx, float T);

    /**
     * @brief Daytime leaf respiration averaged over the half-hourly steps.
     *
     * @param ctx Simulation context (global model state).
     * @param T  Factor applied to the half-hourly temperatures (varday_T).
     * @return Daytime leaf respiration averaged over a day.
     */
    float dailyRdayleaf(Context &ctx, float T);

    /**
     * @brief Determine the maximum LAI of the tree.
     *
     * LAImax (maximum LAI) is assumed to lie between 0 and 10, and is such
     * that adding more leaves results in a net loss of carbon. The range is
     * narrowed down by bisection (sequentially halving the possible range).
     * Uses actually absorbed PPFD instead of incident PPFD.
     *
     * @param ctx Simulation context (global model state).
     */
    void CalcLAImax(Context &ctx);

    /**
     * @brief Determine the maximum leaf area for the current environment.
     *
     * Similar to CalcLAImax, but directly calculates the maximum leaf area for
     * the current light environment the tree experiences rather than for a
     * theoretical day.
     *
     * @param ctx                 Simulation context (global model state).
     * @param LAIexperienced_eff [out] Effective leaf area index experienced at the crown top (m2 m-2).
     * @param LAmax              [out] Maximum leaf area of the tree (m2).
     */
    void CalcLAmax(Context &ctx, float &LAIexperienced_eff, float &LAmax);

    /**
     * @brief Leaf lifespan from the Kikuzawa model.
     *
     * Based on Kikuzawa, K. (1991). A cost-benefit analysis of leaf habit and
     * leaf longevity of trees and their geographical pattern, The American
     * Naturalist 138(5):1250-1263.
     *
     * @param ctx Simulation context (global model state).
     * @return Optimal (timestep-scale) leaf lifespan.
     */
    float predLeafLifespanKikuzawa(Context &ctx);
};

#endif
