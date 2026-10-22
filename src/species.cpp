#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

// #############################################
//  Species constructor
// #############################################
Species::Species()
{
    s_nbind = 0;
    s_nbind10 = 0;
    s_nbind30 = 0;
};

// ######################################
//! Species class function:
//! Initialise Species vector
// ######################################

void Species::Init(Context &ctx)
{

    // !!!: WARNING! seed volume is provided instead of seed mass
    // We assume a conversion factor of 1.0 from wet volume to wet mass (~density of water)
    // We assume a conversion factor of 0.4 from wet mass to dry mass (~40% of the seed are water)
    s_seedmass *= 0.4;
    s_iseedmass = 1.0 / s_seedmass;
    s_ds = 40.0; // !!!UPDATE

    if (ctx.opt._SEEDTRADEOFF)
        s_nbext = (int(s_regionalfreq * ctx.params.Cseedrain * s_iseedmass) + 1);
    else
        s_nbext = int(s_regionalfreq * ctx.params.Cseedrain * (ctx.grid.sites * ctx.grid.LH * ctx.grid.LH / 10000));

    s_nbind = 0;

#ifdef LCP_alternative
#else
    // Computation of the light compensation point from dark respiration and the quantum yield ctx.params.phi
    // By definition, Rdark is in micromolC/ctx.params.m^2/s and it is used in the Species::NPP() routine
    float Vcmax = CalcVcmaxm(s_LMA, s_Nmass, s_Pmass) * s_LMA;
    s_LCP = (Rtot_by_Rabove * Rtotleaf_by_Rdark * CalcRdark(s_LMA, s_Nmass, s_Pmass, Vcmax)) / ctx.params.phi;
#endif

    s_sum1 = s_sum10 = s_sum30 = s_ba = s_ba10 = s_agb = s_gpp = s_npp = s_rday = s_rnight = s_rstem = s_litterfall = 0.0; // new in v.3.1

    // s_g1=-3.97*s_wsg+6.53;                      // from Lin et al. 2015 Nature Climate Change
    // s_dWSF=1.0/(-0.00395145-0.3626778*s_tlp);         // this is the denominator of WSF= phi_sc - s_tlp, where phi_sc (for 'stomatal closure') corresponds to the onset of water stress, derived from the relationship between leaf ctx.params.phi at 50%loss of stomatal closure (ie. WSF~0.5 considering only stomatal responses to water stress) and turgor loss point (drawn from Bartlett et al. 2016 PNAS), and given the shape of WSF factor here assumed
    // s_phi_lethal= -0.9842 + 3.1795*s_tlp;       // Inferred from data provided in Bartlett et al. 2016 PNAS
    // s_itlp=1/s_tlp;
    // s_b=-2.23*s_itlp; // this results from the following assumption: (i) ctx.params.phi(TLP)=0.97*ctx.params.phi(gs90), from Martin-StPaul et al. 2017 Ecology letters; (ii) WSF=0.1 at ctx.params.phi(gs90)
}
