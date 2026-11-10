#include "troll.hpp"
#include "constants.hpp"
#include "context.hpp"


// ######################################
//  Global function: Evolution at each ctx.time.timestep
// ######################################
void Evolution(Context &ctx)
{
#ifdef CHECK_CARBON
    if (ctx.time.iter == 0)
    {
        carbon_assimilated_total = 0.0;
        carbon_net_total = 0.0;
    }
#endif

#ifdef WATER
    ctx.soil.transpiration_1016 = 0.0;
#endif

    UpdateField(ctx); // Update light fields and seed banks
    ctx.diag.nbtrees_n10 = ctx.diag.nbtrees_n30 = ctx.diag.nbdead_n1 = ctx.diag.nbdead_n10 = ctx.diag.nbdead_n30 = 0;
    ctx.diag.nbtrees_carbstarv_n1 = ctx.diag.nbtrees_carbstarv_n10 = ctx.diag.nbtrees_carbstarv_n30 = 0;

#ifdef Output_ABC
    nbdead_n10_abc = 0;
#endif
    if (ctx.opt._BASICTREEFALL)
    {
        // secondary treefalls are triggered first, since they have been caused in the previous iteration
        TriggerTreefallSecondary(ctx); // Compute and distribute Treefall events, caused by treefalls in the previous iteration
        TriggerTreefall(ctx);          // Compute and distribute Treefall events, caused by wind drag
    }

    for (int site = 0; site < ctx.grid.sites; site++)
    {
        //**** Tree evolution: Growth or death ****
        ctx.T[site].Update(ctx);
    }

    for (int d = 0; d < ctx.grid.nbdcells; d++)
    {
        for (int l = 0; l < ctx.soil.nblayers_soil; l++)
        {
            ctx.soil.Transpiration[l][d] = 0.0;
        }
    }
    for (int site = 0; site < ctx.grid.sites; site++)
    {
        ctx.T[site].Water_uptake(ctx);
    }

    // Update trees
    Average(ctx); //! Compute averages for outputs
    if (ctx.opt._OUTPUT_extended)
        OutputField(ctx); //! Output the statistics
}


// Helper function
int GetTimeofyear(Context &ctx)
{
    // new function to derive time of year, extended to negative iterations (-1 would be treated as last iteration of previous year)
    int timeofyear;
    if (ctx.time.iter < 0)
        timeofyear = ctx.time.iterperyear - abs(ctx.time.iter) % ctx.time.iterperyear;
    else
        timeofyear = ctx.time.iter % ctx.time.iterperyear;
    return (timeofyear);
};

