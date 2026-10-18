#ifndef TROLL_DEFINES_HPP
#define TROLL_DEFINES_HPP

// Compile-time model options.
//
// Every option can be changed without editing this file, from CMake:
//   cmake -DTROLL_DEFINES="Output_ABC;TROLL_NO_WATER" ..
// * Options ON by default are switched OFF by defining TROLL_NO_<OPTION>.
// * Options OFF by default are switched ON by defining <OPTION>.
// CI builds every supported combination listed in .github/workflows/build.yml.
// The former MPI option was removed: it had not been maintained since v.2.2 and no longer compiled.

//! new in v.3.1.4: new way of checking whether light environment supports tree birth, dependent on individual seed/seedling's intraspecific variation >
#ifndef TROLL_NO_LCP_alternative
#define LCP_alternative
#endif

// test with positive g0 in stomatal conductance model
#ifndef TROLL_NO_G0
#define G0
#endif

//! small additions and quick fixes that need to be checked/discussed very carefully
#ifndef TROLL_NO_FF_todiscuss
#define FF_todiscuss
#endif

//! new in v.3.0: If defined, an explicit water cycle is added, with an explicit belowground space. The horizontal resolution of the soil field is currently set by DCELL
#ifndef TROLL_NO_WATER
#define WATER
#endif

#ifndef TROLL_NO_PHENO_DROUGHT
#define PHENO_DROUGHT
#endif

//! new in v.2.4.1, modified in v.2.5: If activated, crowns are assumed to grow cylindrical until they reach 3m in depth, then the center of the crown will experience quicker height growth than the outer parts of the crown, thus leading to an umbrella like crown shape (with slope depending on crown depth and width). Depending on the slope parameter, trees will look more like cylinders or more like cones. Contrary to previous versions, however, the crown will not fill its shape underneath the first three layers. The crown will simply become umbrella-like, with a dense layer that spans out from its stem, and empty space underneath. This allows for a very simple computation of just three crown layers, while retaining their spread across the crown depth, enough realism for LiDAR derived CHMs and for eventual non vertical light penetration. New in v.2.5: Entirely modelled through one template that simulates loops across a crown layer/shell and can be applied to any calculation that requires to calculate crown properties (CalcLAI, Fluxh, etc.). Furthermore, enhanced flexibility for future incorporation of other crown shape functions (e.g. spherical, etc.)
#ifndef TROLL_NO_CROWN_UMBRELLA
#define CROWN_UMBRELLA
#endif

//! new in v.2.4.1, modified in v.2.5: a tree's leaf density changes from top to bottom. Broadly based on Kitajima 2004 (Ann Bot), we assume that trees allocate 50% of their leaves to the first meter of their crown, 25% to the second meter layer, and the rest evenly spread across the crown. In case of the umbrella-like crowns, layers correspond to depth from canopy top. In case of a gap in the top layers, there is no shift of leaf density to a lower layer. This ensures that the gaps are real gaps. New in v.2.5: activated within a single function (LAI2dens) as part of the crown template scheme. Small effect in this case, however, since all leaves are concentrated in three shells anyways, with a shift from (0.33,0.33,0.33) to (0.5, 0.25,0.25) likely no strong effect. If implemented for cylindric crowns with full depth, this might, however, have a strong effect, as upper leaves are less impacted by crown overlap than lower leaves.
// LAI_gradient: OFF by default (define LAI_gradient to enable)

//! new in v.2.4.1, refined in v.2.5: PARAMETERIZATION/OUTPUT TOOL, inclusion of ABC routines for comprehensive parameter inference with TROLL
// Output_ABC: OFF by default (define Output_ABC to enable)

//! new in v.2.5: PARAMETERIZATION/OUTPUT TOOL. Since our empirical CHMs are typically derived with a spike-free/pit-free algorithm (e.g. from LAStools), we should also remove spikes/pits from the TROLL-generated CHM. Here, this is simply implemented as considering all crown gaps as filled. A more sophisticated version could use a simulated lidar point cloud and then apply the spike-free/pit-free algorithms or compare the raw CHMs of simulated lidar and actual lidar. Since we are, however, unlikely to capture the microstructure of crowns well (i.e. crown gaps are likely to be smaller, but more frequent than the simulated gaps), the current approach is probably a good compromise.
#ifndef TROLL_NO_CHM_SPIKEFREE
#define CHM_SPIKEFREE
#endif

//! new in v.2.5: DIAGNOSTIC TOOL, individual-based tracking of trees for PBA (process-based analytics, an accounting-like monitoring of model behavior through individual states and processes). PLEASE BE AWARE: onset of tracking is currently hardcoded to comprise trees born in year 501 in forest regeneration. Also: Filesizes can get VERY VERY LARGE if more trees are included.
// TRACK_INDIVIDUALS: OFF by default (define TRACK_INDIVIDUALS to enable)

//! new in v.2.5: DIAGNOSTIC TOOL, checking of carbon budgets, could potentially be extended for nutrient budget checking in the future. The idea is to keep track of carbon stocks and carbon fluxes every timestep to see whether there are any deviations from expectations - to do so, differences between stocks are computed at each timestep, and can be compared to the gross and net assimilation of carbon
// CHECK_CARBON: OFF by default (define CHECK_CARBON to enable)

#ifndef TROLL_NO_FULL_CLIMATE
#define FULL_CLIMATE
#endif

//! includes specific developments and outputs needed for the MIP experiment led by Jeremy Lichstein.
// MIP_Lichstein: OFF by default (define MIP_Lichstein to enable)

//! cache Look Up Tables (LUT) calculations
// CACHE_LUT: OFF by default (define CACHE_LUT to enable)

#endif // TROLL_DEFINES_HPP
