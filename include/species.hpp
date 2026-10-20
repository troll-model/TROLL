#ifndef TROLL_SPECIES_HPP
#define TROLL_SPECIES_HPP

#include <string>
#include <vector>
#include "troll_defines.hpp"

// Forward declaration — full definition in context.hpp
struct Context;

// ###########################################
//! Species class
// ############################################

class Species
{

public:
    int s_nbind;   //!< Number of individuals per species
    int s_nbind10; //!< Number of individuals > 10cm per species
    int s_nbind30; //!< Number of individuals > 30cm per species
    int s_nbext;   //!< Total number of incoming seeds in the simulated plot at each timestep (seed rain) -- v.2.2
    std::string s_name; //!< Species name, updated to string in v.3.0 (higher flexibility)

#ifdef LCP_alternative
#else
    float s_LCP; //!< Light compensation point  (micromol photon/m^2/s)
#endif
    float s_ds;           //!< Mean seed dispersal distance (m)
    float s_dbhmax;       //!< Maximal dbh (m) !! data usually report dbh in cm
    float s_hmax;         //!< Maximal height (m)
    float s_LMA;          //!< Leaf mass per area (g/m^2)
    float s_Nmass;        //!< Leaf nitrogen concentration (g/g) v.2.01
    float s_Pmass;        //!< Leaf phosphorous concentration (g/g) v.2.01
    float s_wsg;          //!< Wood specific gravity (g/cm^3)
    float s_ah;           //!< Parameter for allometric height-dbh equation
    float s_regionalfreq; //!< Regional frequency; v.3.0 !!!UPDATE
    float s_leafarea;     //!< Drymass; v.3.0  !!!UPDATE
    float s_seedmass;     //!< Seed mass (g); See Baraloto & Forget 2007 dataset v.2.3; deprecated in v.2.2, but still necessary for SEEDTRADEOFF
    float s_iseedmass;    //!< Inverse of seed mass (1/g), v.2.3
    // float s_output_field[12];         // scalar output fields, deprecated since v.3.1, replaced by actual sumstats for readability/code accessibility
    // Per-species sums for sumstats_species.txt and sumstats.txt, reset each timestep
    float s_sum1; //!< Living trees of the species (output sums, before per-hectare scaling)
    float s_sum10; //!< Trees with DBH >= 10 cm
    float s_sum30; //!< Trees with DBH >= 30 cm
    float s_ba; //!< Basal area, all trees (m2)
    float s_ba10; //!< Basal area, trees with DBH >= 10 cm (m2)
    float s_agb; //!< Aboveground biomass (kg)
    float s_gpp; //!< GPP in the timestep (gC)
    float s_npp; //!< NPP in the timestep (gC)
    float s_rday; //!< Daytime leaf respiration in the timestep (gC)
    float s_rnight; //!< Night-time leaf respiration in the timestep (gC)
    float s_rstem; //!< Stem respiration in the timestep (gC)
    float s_litterfall; //!< Leaf litterfall in the timestep (g dry mass)

    float s_tlp; //!< Leaf water potential at turgor loss point (MPa); defined for consistency when WATER is deactivated
#ifdef WATER     // Some of these parameters may include intraspecific variability, as in v.2.4.1.
    // float s_g1;
    // float s_dWSF;
    // float s_phi_lethal,           //!< plant water potential at drought-induced death (MPa)
    // s_itlp,         //!< inverse of s_tlp
    // s_b;            //!< parameter used to compute the water stress factor (WSF) for stomatal limitation
#endif

#ifdef Output_ABC
    float s_dbhmax_realized; //!< Maximum dbh found across individuals from this species
#endif


    //! Constructor of the Species class
    Species();

    void Init(Context &ctx);
};

#endif
