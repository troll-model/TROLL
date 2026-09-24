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
    float carbon_flux;
    float water_flux;
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
    int t_inInventory;
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
    float t_g1_0;
    float t_g1;

#ifdef PHENO_DROUGHT
    int t_Ndays_dry;
    int t_Ndays_wet;
    float t_pheno_factor;
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
    int t_month_born; //!< Not documented  !!!UPDATE

    int t_seedsproduced;
    int t_seedsproduced_sumyear;
    int t_time_carbonstarvation;
    int t_time_carbonstarvation_year;

    float t_GPP_sumyear;
    float t_NPP_sumyear;
    float t_GPPsquared_sumyear;
    float t_NPPsquared_sumyear;
    float t_Rday_sumyear;
    float t_Rnight_sumyear;
    float t_Rstem_sumyear;
    float t_LAIabove_effavgyear;
    float t_carbon_storage_avgyear;

    float t_dbh_tracked;
    float t_height_tracked;
    float t_CR_tracked;
    float t_agb_tracked;

    float t_LAIcum;
    float t_LAIeffcum;
    float t_GPPcum;
    float t_NPPcum;
    float t_LAIsquared_cum;
    float t_LAIeffsquared_cum;
    float t_GPPsquared_cum;
    float t_NPPsquared_cum;

#endif
    // Constructor
    Tree(Context &ctx);

    void Birth(Context &ctx, int, int);
    int BirthFromInventory(Context &ctx, int site, std::vector<std::string> &parameter_names, std::vector<std::string> &parameter_values, int &nb_speciesrandom);
    void Death(Context &ctx);
    void Growth(Context &ctx);
#ifdef WATER
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer, float &PPFD_incident, float &ExtinctLW);
#else
    void Fluxh(Context &ctx, int h, float &PPFD, float &VPD, float &Tmp, float &leafarea_layer);
#endif

#ifdef WATER
    void Water_availability(Context &ctx);
    void Water_uptake(Context &ctx);
#endif

    void CalcRespGPP(Context &ctx);
    void CalcNPP(Context &ctx);
    void UpdateLeafDynamics(Context &ctx);
    void UpdateTreeBiometry(Context &ctx);
    void UpdateVolumeDensity();

    void Treefall(Context &ctx, float angle);
    void Update(Context &ctx);

    float CalcJmaxm(Context &ctx);
    float CalcAGB(Context &ctx);
    float CalcIncrementDBH(Context &ctx, float delta_agb);
    float CalcCarbonStorageMax(Context &ctx);
    float CalcCt(Context &ctx);

    void CalcLeafLifespan(Context &ctx);
    void InitialiseLeafPools(Context &ctx);
    void UpdateSapwoodArea(Context &ctx, float ddbh);
    void UpdateHeight(Context &ctx);
    void UpdateCR(Context &ctx);
    void UpdateCD(Context &ctx);

    float DeathRateNDD(Context &ctx, float, float, float);
#ifdef WATER
    float DeathRate(Context &ctx, float, float, float);
    leafFluxes FluxesLeaf(Context &ctx, float, float, float, float, float, float, float, float, float);
    leafFluxes dailyFluxesLeaf(Context &ctx, float, float, float, float, float, float, float &);
    leafFluxes Photosyn(Context &ctx, float, float, float, float);
#else
    float DeathRate(Context &ctx, float, float);
    float GPPleaf(Context &ctx, float, float, float);
    float dailyGPPleaf(Context &ctx, float, float, float);
    float dailyGPPcrown(Context &ctx, float, float, float, float);
#endif
    float Rdayleaf(Context &ctx, float T);
    float dailyRdayleaf(Context &ctx, float T);

    void CalcLAImax(Context &ctx);
    void CalcLAmax(Context &ctx, float &LAIexperienced_eff, float &LAmax);
    float predLeafLifespanKikuzawa(Context &ctx);
};

#endif
