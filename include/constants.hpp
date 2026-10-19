#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <vector>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_double.h>

// Global constants (e.g. PI and derivatives, physical constant etc)
# define PI 3.141592654  //!< Global constant: Pi
# define twoPi 6.2831853071 //!< Global constant: 2.0*Pi
# define Pis2 1.570796327 //!< Global constant: Pi/2.0
# define iPi 0.3183099 //!< Global constant: 1.0/Pi

// Constants used only when WATER is defined
#define ABSZERO -273.15    //!< Minimal temperature (in °C)
#define AIRMA 0.02896279  //!< molelcular mass of air (kg/mol)
#define CPAIR 1010.0      //!< heat capacity of air (J/kg/K)
#define EMLEAF 0.95       //!< emissivity of thermal radiation by leaf
#define SIGMA 5.67e-8     //!< Stefan-Boltzmann constant (W m-2 K-4)
#define H2OLV0 2.501e6    //!< latent heat H2O (J/kg)
#define H2OMW 18.0e-3     //!< mol mass H2O (kg/mol)
#define RCONST 8.314      //!< universal gas constant (J/mol/K)
#define DHEAT 21.5e-6     //!< molecular diffusivity to heat (m2/s)
#define GBHGBC 1.32       //!< Ratio of Gbh:Gbc (value consistent with Vezy et al. 2018, double-checked by R. Vezy (cf. mail 21/10/2020))
#define GBVGBH 1.075      //!< Ratio of Gbw:Gbh (relative diffusivities of heat and water vapor in air) (cf. equ. A7 in Medlyn et al. 2007; equ. E5 in Leuning et al. 1995)
#define GSVGSC 1.57       //!< Ratio of Gsw:Gsc (relative diffusivities of CO2 and water vapor in air)
#define LSQ 0.15          //!< Light spectral quality: effective spectral quality of light
#define DAYRESP 0.4       //!< Inhibition of leaf dark respiration in day light by ca. 40%, cf. Atkin et al. 2000. Note: should be given in the input file probably

// Single-precision versions for the leaf-flux kernel (Tree::FluxesLeaf, Tree::Photosyn, QUAD).
// All kernel variables are float; with the double constants above, every expression was
// evaluated in double and converted back, which cost ~40% of simulation time (TROLL 4.1.1).
// Keep the kernel in float: use these constants and f-suffixed literals there
// (the kernel files are compiled with -Wdouble-promotion to catch regressions).
constexpr float ABSZERO_F = -273.15f;
constexpr float AIRMA_F = 0.02896279f;
constexpr float CPAIR_F = 1010.0f;
constexpr float H2OLV0_F = 2.501e6f;
constexpr float H2OMW_F = 18.0e-3f;
constexpr float RCONST_F = 8.314f;
constexpr float DHEAT_F = 21.5e-6f;
constexpr float GBHGBC_F = 1.32f;
constexpr float GBVGBH_F = 1.075f;
constexpr float GSVGSC_F = 1.57f;
constexpr float DAYRESP_F = 0.4f;
constexpr float CPAIR_AIRMA_F = CPAIR_F * AIRMA_F;   //!< heat capacity of air per mole (J/mol/K)
constexpr float GBH_FREE_F = 56.2341f * DHEAT_F;     //!< free-convection coefficient for boundary-layer conductance

// the following parameter values could be provided in input
#define Rtot_by_Rabove 1.5 //!< Multiplier coefficient used to roughly account for belowground respiration (coarse and fine root respiration)
#define Rtotleaf_by_Rdark 1.4 //!< Multiplier coefficient used to roughly account for leaf respiration during the day
#define alphaInoue 3.0 //!< alpha coefficient in the model of wind speed profile by Inoue (1963). Values provided in Raupach et al. 1991, Shaw and Pereira (1982), Shuttleworth et al. (1989) (see Rau et al's manuscript)
#define MeteoStation_Height 55.0 //!< Height of the meteo station where the input climate data are taken from. Here supposed to be above the canopy.

#endif // CONSTANTS_HPP