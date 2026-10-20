#include <iterator>
#include <string>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include "lookup_tables.hpp"
#include "internal/lookup_tables_cache.hpp"

#include <fstream>
#include <iostream>
#include <filesystem> // C++17 feature

//! The tables depend on the input parameters, so the cache file name contains a
//! hash of the input files and of the compile-time options that change the layout.
// (A single fixed file name silently reused tables computed for other inputs.)
//! @param ctx simulation context (input file names)
static std::string LutCacheFile(const Context &ctx)
{
    std::uint64_t hash = 1469598103934665603ULL; // FNV-1a
    auto mix = [&hash](const char *data, std::size_t n)
    {
        for (std::size_t i = 0; i < n; i++)
        {
            hash ^= (unsigned char)data[i];
            hash *= 1099511628211ULL;
        }
    };
    const char *files[] = {ctx.fileio.inputfile, ctx.fileio.inputfile_species, ctx.fileio.inputfile_climate,
                           ctx.fileio.inputfile_daytimevar, ctx.fileio.inputfile_soil};
    for (const char *file : files)
    {
        std::ifstream in(file, std::ios::binary);
        std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        mix(content.data(), content.size());
        mix("|", 1);
    }
#ifdef WATER
    mix("WATER", 5);
#endif
    char name[64];
    std::snprintf(name, sizeof(name), "./.troll_cache/lut_%016llx.bin", (unsigned long long)hash);
    return name;
}
//! Format tag at the start of the cache file; caches without it are rejected
static const char LUT_CACHE_MAGIC[8] = {'T', 'R', 'O', 'L', 'L', 'L', 'T', '2'}; // bump when the format changes

// --- Save LUTs to file -------------------------------------------------------

//! Write the lookup tables to the cache file (CACHE_LUT option); false on failure
//! @param ctx simulation context (all model state)
bool SaveLookUpTablesToCache(Context &ctx)
{

    // Create the directory if it doesn't exist
    std::filesystem::path cacheDir = ".troll_cache";
    if (!std::filesystem::exists(cacheDir))
    {
        if (!std::filesystem::create_directory(cacheDir))
        {
            std::cerr << "Failed to create directory: " << cacheDir << "\n";
            return false;
        }
    }

    const std::string cache_file = LutCacheFile(ctx);
    std::ofstream out(cache_file, std::ios::binary);
    if (!out)
    {
        std::cerr << "Cannot open LUT cache for writing: " << cache_file << "\n";
        return false;
    }

    // 0) Format header: caches written by older versions lack the accuracy
    //    scalars below and are rejected (and recomputed) on load.
    out.write(LUT_CACHE_MAGIC, sizeof(LUT_CACHE_MAGIC));

    // 1) Scalar sizes and bin accuracies (the latter were previously not cached,
    //    so after a cache load every lookup used bin 0)
    out.write((char *)&ctx.lookup.nbTbins, sizeof(ctx.lookup.nbTbins));
    out.write((char *)&ctx.lookup.nbVPDbins, sizeof(ctx.lookup.nbVPDbins));
    out.write((char *)&ctx.lookup.nbHbins, sizeof(ctx.lookup.nbHbins));
    out.write((char *)&ctx.lookup.iTaccuracy, sizeof(ctx.lookup.iTaccuracy));
    out.write((char *)&ctx.lookup.iVPDaccuracy, sizeof(ctx.lookup.iVPDaccuracy));
    out.write((char *)&ctx.lookup.iHaccuracy, sizeof(ctx.lookup.iHaccuracy));

    // 2) Arrays of size ctx.lookup.nbTbins
    out.write((char *)ctx.lookup.LookUp_KmT.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_GammaT.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_VcmaxT.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_JmaxT.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_Rleaf.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_Rstem.data(), ctx.lookup.nbTbins * sizeof(float));

#ifdef WATER
    // 3) WATER-dependent LUTs
    out.write((char *)ctx.lookup.LookUp_SLOPE.data(), ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_GRADN.data(), ctx.lookup.nbTbins * sizeof(float));

    // INLR: ctx.lookup.nbTbins × ctx.lookup.nbVPDbins
    for (int i = 0; i < ctx.lookup.nbTbins; ++i)
        out.write((char *)ctx.lookup.LookUp_INLR[i].data(), ctx.lookup.nbVPDbins * sizeof(float));

    // Wind LUT
    out.write((char *)ctx.lookup.LookUp_Wind.data(), ctx.lookup.nbHbins * sizeof(float));
#endif

    // 4) Flux LUTs: all size 80000
    out.write((char *)ctx.lookup.LookUp_flux_absorption.data(), 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_flux.data(), 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_ExtinctLW.data(), 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_VPD.data(), 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_T.data(), 80000 * sizeof(float));

    // 5) Crown site: 2601 ints
    out.write((char *)ctx.lookup.LookUp_Crown_site, 2601 * sizeof(int));

    return out.good();
}

// --- Load LUTs from file -----------------------------------------------------

//! Read the lookup tables from the cache file (CACHE_LUT option); false if missing, outdated or unreadable
//! @param ctx simulation context (all model state)
bool LoadLookUpTablesFromCache(Context &ctx)
{
    std::ifstream in(LutCacheFile(ctx), std::ios::binary);
    if (!in)
        return false; // cache does not exist

    // 0) Format header
    char magic[sizeof(LUT_CACHE_MAGIC)];
    in.read(magic, sizeof(magic));
    if (!in.good() || std::memcmp(magic, LUT_CACHE_MAGIC, sizeof(magic)) != 0)
        return false; // missing, old or foreign cache file: recompute

    // 1) Scalar sizes and bin accuracies
    in.read((char *)&ctx.lookup.nbTbins, sizeof(ctx.lookup.nbTbins));
    in.read((char *)&ctx.lookup.nbVPDbins, sizeof(ctx.lookup.nbVPDbins));
    in.read((char *)&ctx.lookup.nbHbins, sizeof(ctx.lookup.nbHbins));
    in.read((char *)&ctx.lookup.iTaccuracy, sizeof(ctx.lookup.iTaccuracy));
    in.read((char *)&ctx.lookup.iVPDaccuracy, sizeof(ctx.lookup.iVPDaccuracy));
    in.read((char *)&ctx.lookup.iHaccuracy, sizeof(ctx.lookup.iHaccuracy));
    if (!in.good())
        return false;

    // 2) Allocate raw arrays
    ctx.lookup.LookUp_KmT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_GammaT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_VcmaxT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_JmaxT.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_Rleaf.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_Rstem.assign(ctx.lookup.nbTbins, float());

#ifdef WATER
    ctx.lookup.LookUp_SLOPE.assign(ctx.lookup.nbTbins, float());
    ctx.lookup.LookUp_GRADN.assign(ctx.lookup.nbTbins, float());

    ctx.lookup.LookUp_INLR.resize(ctx.lookup.nbTbins);
    for (int i = 0; i < ctx.lookup.nbTbins; i++)
        ctx.lookup.LookUp_INLR[i].assign(ctx.lookup.nbVPDbins, float());

    ctx.lookup.LookUp_Wind.assign(ctx.lookup.nbHbins, float());
#endif

    ctx.lookup.LookUp_flux_absorption.assign(80000, float());
    ctx.lookup.LookUp_flux.assign(80000, float());
    ctx.lookup.LookUp_ExtinctLW.assign(80000, float());
    ctx.lookup.LookUp_VPD.assign(80000, float());
    ctx.lookup.LookUp_T.assign(80000, float());
    // ctx.lookup.LookUp_Crown_site is a static array; no allocation needed.

    // 3) Read back LUT data
    in.read((char *)ctx.lookup.LookUp_KmT.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_GammaT.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_VcmaxT.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_JmaxT.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_Rleaf.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_Rstem.data(), ctx.lookup.nbTbins * sizeof(float));

#ifdef WATER
    in.read((char *)ctx.lookup.LookUp_SLOPE.data(), ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_GRADN.data(), ctx.lookup.nbTbins * sizeof(float));

    for (int i = 0; i < ctx.lookup.nbTbins; ++i)
        in.read((char *)ctx.lookup.LookUp_INLR[i].data(), ctx.lookup.nbVPDbins * sizeof(float));

    in.read((char *)ctx.lookup.LookUp_Wind.data(), ctx.lookup.nbHbins * sizeof(float));
#endif

    in.read((char *)ctx.lookup.LookUp_flux_absorption.data(), 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_flux.data(), 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_ExtinctLW.data(), 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_VPD.data(), 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_T.data(), 80000 * sizeof(float));

    in.read((char *)ctx.lookup.LookUp_Crown_site, 2601 * sizeof(int));

    return in.good();
}