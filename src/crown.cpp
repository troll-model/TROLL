#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"

////////////////////////////////////////////////
/////// GLOBAL FUNCTIONS ///////////////////////
////////////////////////////////////////////////

#ifdef G0

// Solves the quadratic equation. ROOT: finds smaller (-1), larger (1) root.
float QUAD(float A, float B, float C, int ROOT)
{
    float quad = 0.0;
    float discriminant = B * B - 4.0 * A * C;
    // ROOT=-1: smaller root; ROOT=1: larger root
    if (discriminant < 0)
    {
        cerr << "IMAGINARY ROOTS IN QUADRATIC" << endl;
        quad = 0.0;
    }
    if (A == 0.0)
    {
        if (B == 0.0)
            quad = 0.0;
        else
            quad = -C / B;
    }
    else
        quad = (-B + ROOT * sqrt(discriminant)) / (2.0 * A);
    return quad;
}

#endif

//! Global functions to compute Vcmax (on a mass basis) and leaf dark respiration (on a leaf area basis)
// newIM : to avoid repeating same functions in Tree and Species class. Was needed for the Species class to correctly initialise Rdark and compute LCP

float CalcVcmaxm(float LMA, float Nmass, float Pmass)
{
    float SLA = 10000.0 / LMA;                                                                                                                               // in cm2 g-1
    float Vcmaxm = pow(10.0, fminf((-1.56 + 0.43 * log10(Nmass * 1000.0) + 0.37 * log10(SLA)), (-0.80 + 0.45 * log10(Pmass * 1000.0) + 0.25 * log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA). Nmass and Pmass are given in g g-1, but should be in mg g-1 in equ 2 in Domingues et al. 2010, hence the mutiplication by 1000.
    return (Vcmaxm);                                                                                                                                         //  in micromol C g-1 s-1
}
float CalcRdark(float LMA, float Nmass, float Pmass, float Vcmax)
{
    float Parea = Pmass * LMA;                                                                         // in g ctx.params.m-2
    float Narea = Nmass * LMA;                                                                         // in g ctx.params.m-2
    float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2)); // in micromolC ctx.params.m-2 s-1 //since v.2.5: correction of Atkin et al. 2015 New phytologist formula. The original formula was based on mean-centered values (cf. Atkin et al. 2015, and the clarification/correction published afterwards), unfortunately only one formula (PFT-specific formula) with absolute values is provided and this is the one used here, cf. corrigendum TableS4): https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.13253&file=nph13253-sup-0001-SupInfo.pdf Warning: Vcmax should be provided on an area basis.
    return (Rdark);
}

#ifdef CROWN_UMBRELLA

// Global function: linear decrease of crown radius
float GetRadiusSlope(Context &ctx, float CR, float crown_extent, float crown_position)
{
    float crown_slope = CR * (1.0 - ctx.crown.shape_crown) / crown_extent;
    float radius = CR - crown_slope * float(crown_position);
    return (radius);
}

// Global function: not currently used, but returns the input radius
float GetRadiusCylinder(float CR, float crown_extent, float crown_position)
{
    return (CR);
}

// Global function: converts floating point crown area into integer value, imposing lower and upper limits
int GetCrownIntarea(float crown_radius)
{
    // crown area
    float crown_area = PI * crown_radius * crown_radius;
    int crown_intarea = int(crown_area);      // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea, 1);    // minimum area of crown (1)
    crown_intarea = min(crown_intarea, 1963); // maximum area of crown (radius 25), int(3.14*25*25)
    return (crown_intarea);
}

// Global function: deduces within-crown densities from LAI with a gradient from 50% in top ctx.diag.layer to 25% in belowtop and 25% in all shells underneath (1 ctx.diag.layer for umbrella-like shape)
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens_layer) // RENAMED: ctx.params.dens → dens_layer (was shadowing global ctx.params.dens/ctx.params.ctx.params.dens)
{
    if (CD < 2.0)
    {
        dens_top = dens_belowtop = dens_layer = LAI / CD;
    }
    else if (CD < 3.0)
    {
        dens_top = 0.5 * LAI;
        dens_belowtop = dens_layer = 0.5 * LAI / (CD - 1.0);
    }
    else
    {
        dens_top = 0.5 * LAI;
        dens_belowtop = 0.25 * LAI;
#ifdef CROWN_UMBRELLA
        dens_layer = 0.25 * LAI;
#else
        dens_layer = 0.25 * LAI / (CD - 2.0);
#endif
    }
}

// Global function: deduces within-crown density from LAI, assuming uniform leaf distribution
void GetDensityUniform(float LAI, float CD, float &dens_layer) // RENAMED: ctx.params.dens → dens_layer
{
#ifdef CROWN_UMBRELLA
    float crownshells_limit = fminf(CD, 3.0);
    dens_layer = LAI / crownshells_limit;
#else
    dens_layer = LAI / CD;
#endif
}

// Global function: dummy function when no modification is needed
void KeepFloatAsIs(float input, float &output, float CD, float height, int layer_fromtop)
{
    output = input;
}

void KeepIntAsIs(int input, int &output, float CD, float height, int layer_fromtop)
{
    output = input;
}

// Global function: a modifying function that converts LAI to the density of a specific ctx.diag.layer, using the GetDensity functions
//! - modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results
//! - LAI is the input, dens_layer the output
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop)
{
    int crown_top = int(height);
    int crown_base = int(height - CD);
    float dens_top, dens_belowtop, dens_lowerlayer; // RENAMED: dens → dens_lowerlayer (base-layer crown density from GetDensitiesGradient; distinct from output param dens_layer)
#ifdef LAI_gradient
    GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens_lowerlayer);
#else
    GetDensityUniform(LAI, CD, dens_lowerlayer);
    dens_top = dens_belowtop = dens_lowerlayer;
#endif

    if (CD < 3.0 && crown_top == crown_base)
    {
        dens_layer = LAI; /* full LAI allocation */
    }
    else if (CD < 3.0 && (crown_top - layer_fromtop == crown_base))
    {
        dens_layer = LAI;
    }
    else
    {
        float fraction_layer = height - floor(height); /* this is the fraction that each ctx.diag.layer apart from the topmost ctx.diag.layer will extend into the voxel above */
        if (layer_fromtop == 0)
            dens_layer = dens_top * fraction_layer;
        else if (layer_fromtop == 1)
            dens_layer = dens_top + dens_belowtop * fraction_layer;
        else if (layer_fromtop == 2)
            dens_layer = dens_top + dens_belowtop + dens_lowerlayer * fraction_layer;
        else
            dens_layer = LAI;
    }
}

// Global function: a modifying function that converts LAI to the density of a specific ctx.diag.layer, using the GetDensity functions
//! - modifier for GPP calculation where we need the leaves per ctx.diag.layer to weight our results
//! - LAI is the input, dens_layer the output
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop)
{

    int crown_top = int(height);
    int crown_base = int(height - CD);
    float dens_top, dens_belowtop, dens_lowerlayer; // RENAMED: dens → dens_lowerlayer (base-layer crown density from GetDensitiesGradient; distinct from output param dens_layer)
#ifdef LAI_gradient
    GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens_lowerlayer);
#else
    GetDensityUniform(LAI, CD, dens_lowerlayer);
    dens_top = dens_belowtop = dens_lowerlayer;
#endif

    if (CD < 3.0 && crown_top == crown_base)
    {
        dens_layer = dens_top * CD;
    }
    else if (CD < 3.0 && (crown_top - layer_fromtop == crown_base))
    {
        float fraction_belowbase = float(crown_base + 1) - (height - CD);
        dens_layer = dens_lowerlayer * fraction_belowbase;
    }
    else
    {
        float fraction_layer = height - floor(height);         // this is the fraction that each ctx.diag.layer apart from the topmost ctx.diag.layer will extend into the voxel above
        float fraction_layer_fromabove = 1.0 - fraction_layer; // the inverse of the fraction above

        if (layer_fromtop == 0)
            dens_layer = dens_top * fraction_layer;
        else if (layer_fromtop == 1)
            dens_layer = dens_top * fraction_layer_fromabove + dens_belowtop * fraction_layer;
        else if (layer_fromtop == 2)
            dens_layer = dens_belowtop * fraction_layer_fromabove + dens_lowerlayer * fraction_layer;
        else
            dens_layer = dens_lowerlayer * fraction_layer_fromabove;
    }
}

// Global function: update of ctx.field.LAI3D field, called by CalcLAI()
void UpdateLAI3D(Context &ctx, int height, int site, float dens_layer, float &LA_cumulated) // RENAMED: ctx.params.dens → dens_layer
{
    ctx.field.LAI3D[height][site + ctx.grid.SBORD] += dens_layer;
    LA_cumulated += dens_layer;
}

// Global function: remove outliers in canopy height model (CHM); vector option
#ifdef CHM_SPIKEFREE
void UpdateCHMvector(int height, int site, float noinput, vector<int> &chm)
{
    if (chm[site] < height)
        chm[site] = height;
}
// Global function: remove outliers in canopy height model (CHM)
void UpdateCHM(int height, int site, float noinput, int *chm)
{
    if (chm[site] < height)
        chm[site] = height;
}
#endif

void OutputCrownSliced(Context &ctx, int height, int site, int row_slice, vector<float> &output_statistics)
{
    int row_current = site / ctx.grid.cols;
    int col_current = site % ctx.grid.cols;
    if (row_current == row_slice && col_current >= ctx.crown.mincol_visual && col_current < ctx.crown.maxcol_visual)
    {
        ctx.out.output_visual[1] << ctx.time.iter << "\t" << row_current << "\t" << col_current << "\t" << height;
        for (size_t i = 0; i < output_statistics.size(); i++)
        {
            ctx.out.output_visual[1] << "\t" << output_statistics[i];
        } // we only output tree parts that fall in the current slice extent
        ctx.out.output_visual[1] << endl;
    }
};

// Global function: PPFD retrieval for function leafarea_max()
void GetPPFDabove(Context &ctx, int height, int site, float noinput, float (&ppfd_CA)[2])
{
    // First get voxel field densities
    float absorb_prev = ctx.field.LAI3D[height + 1][site + ctx.grid.SBORD];
    int intabsorb = CalcIntabsorb(absorb_prev);

    // Obtain PPFD for the voxel, and also record the circled area
    ppfd_CA[0] += ctx.climate.WDailyMean * ctx.lookup.LookUp_flux[intabsorb];
    ppfd_CA[1] += 1.0; // add area
}

// Global function: calculates the canopy environment
//! - this function adds to the environmental variables provided in canopy_environment_cumulated
//! - the PPFD, VPD, Tmp and leafarea_layer retrieval function for Fluxh()
#ifdef WATER
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[6]) // RENAMED: ctx.params.dens → dens_layer
{
#else
void GetCanopyEnvironment(Context &ctx, int height, int site, float dens_layer, float (&canopy_environment_cumulated)[4]) // RENAMED: ctx.params.dens → dens_layer
{
#endif
    // first get voxel field densities
    float absorb_prev = ctx.field.LAI3D[height + 1][site + ctx.grid.SBORD];
    float absorb_curr = ctx.field.LAI3D[height][site + ctx.grid.SBORD];
    float absorb_delta = absorb_curr - absorb_prev;
    if (absorb_delta < 0.0)
        absorb_delta = 0.0; // eliminate rounding errors
    int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);
#ifdef WATER
    int intincident = CalcIntabsorb(absorb_prev);
#endif

#ifdef FULL_CLIMATE
    // Obtain PPFD, VPD and T for the voxel
    float PPFD_voxel = ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD_voxel = ctx.lookup.LookUp_VPD[intabsorb];
    float T_voxel = ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float PPFD_voxel_incident = ctx.lookup.LookUp_flux[intincident];
    float ExtinctLW_voxel = ctx.lookup.LookUp_ExtinctLW[intincident];
#endif // WATER

#else // FULL_CLIMATE
      // Obtain PPFD, VPD and T for the voxel
    float PPFD_voxel = ctx.climate.WDailyMean * ctx.lookup.LookUp_flux_absorption[intabsorb];
    float VPD_voxel = ctx.climate.VPDDailyMean * ctx.lookup.LookUp_VPD[intabsorb];
    float T_voxel = ctx.climate.tDailyMean - ctx.lookup.LookUp_T[intabsorb];
#ifdef WATER
    float PPFD_voxel_incident = ctx.climate.WDailyMean * ctx.lookup.LookUp_flux[intincident];
    float ExtinctLW_voxel = ctx.lookup.LookUp_ExtinctLW[intincident];
#endif // WATER

#endif // FULL_CLIMATE

    // Add the three variables up, weighted by leaf density inside voxel
    canopy_environment_cumulated[0] += dens_layer;
    canopy_environment_cumulated[1] += PPFD_voxel * dens_layer;
    canopy_environment_cumulated[2] += VPD_voxel * dens_layer;
    canopy_environment_cumulated[3] += T_voxel * dens_layer;
#ifdef WATER
    canopy_environment_cumulated[4] += PPFD_voxel_incident * dens_layer;
    canopy_environment_cumulated[5] += ExtinctLW_voxel * dens_layer;
#endif
}

// Global function: calculates packing densities
void AddCrownVolumeLayer(Context &ctx, int row_center, int col_center, float height, float CR, float CD, int crownvolume[70])
{
    int crown_top = int(height);
    int crown_base = int(height - CD);
    if (CD <= 3.0)
    {
        // for the smallest crowns it is simply cylinder rings being filled up
        int crown_intarea = GetCrownIntarea(CR);
        for (int h = crown_top; h >= crown_base; h--)
            crownvolume[h] += crown_intarea;
    }
    else
    {
        // For the rest of the crown, we go through different crown shells. We separate out the innermost sector (a slowly increasing cylinder), and the surrounding parts of the crown
        // first the metrics with respect to the internal crown structure (i.e. z coordinate with respect to crown base)
        float crownshell_base = height - CD + 2.0;                   // lower reference point for the crown slope function is two layers up from the crown base
        float crownshell_extent = height - crownshell_base;          // this is the extent from the "base ctx.diag.layer" to the top
        float crownshell_extent_toplayer = floor(crownshell_extent); // this is the extent to the lower limit of the toplayer
        // then we translate the crown coordinates into discretised variables with respect to the absolute location in the voxel field, as needed for location in the voxel field, with layers defined from top to bottom
        int shell_fromtop = 0;
        int height_innermost = crown_top - shell_fromtop;
        int height_toplayer = int(crownshell_base + crownshell_extent_toplayer) - shell_fromtop;
        int height_baselayer = int(crownshell_base + 1.0) - shell_fromtop;

        // now do calculations
        // first the inner crown shell section that grows dynamically
        float radius_innermost = GetRadiusSlope(ctx, CR, crownshell_extent, crownshell_extent_toplayer);
        int crown_intarea_innermost = GetCrownIntarea(radius_innermost);
        for (int h = height_innermost; h >= crown_base; h--)
        {
            crownvolume[h] += crown_intarea_innermost;
        }
        // now loop through the outer crown shell cylinders
        for (int h_outer = height_toplayer; h_outer >= crown_base; h_outer--)
        {
            // calculating the radius of the current ctx.diag.layer depending on the respective slopes, to be replaced by function
            // float radius_height = CR - crown_slope * (h_outer - height_baselayer);    // for the lowest ctx.diag.layer, i.e. h == height_baselayer, radius = t_CR
            int extent_layerouter = max(h_outer - height_baselayer, 0); // we also fill up underneath the baselayer
            float radius_height = GetRadiusSlope(ctx, CR, crownshell_extent, extent_layerouter);
            int crown_intarea = GetCrownIntarea(radius_height);

            crownvolume[h_outer] += (crown_intarea - crown_intarea_innermost);
        }
    }
}
#endif
