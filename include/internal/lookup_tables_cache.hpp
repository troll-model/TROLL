#ifndef LOOKUPTABLES_CACHE_HPP
#define LOOKUPTABLES_CACHE_HPP

struct Context;

bool LoadLookUpTablesFromCache(Context &ctx); //!< Read the lookup tables from the cache (CACHE_LUT); false if unavailable
bool SaveLookUpTablesToCache(Context &ctx);   //!< Write the lookup tables to the cache (CACHE_LUT); false on failure

#endif // LOOKUPTABLES_CACHE_HPP
