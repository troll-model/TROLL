#include "lookup_tables.hpp"
#include "internal/lookup_tables_cache.hpp"

#include <fstream>
#include <iostream>
#include <filesystem> // C++17 feature

static const char *LUT_CACHE_FILE = "./.troll_cache/lut.bin";

// --- Save LUTs to file -------------------------------------------------------

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

    std::ofstream out(LUT_CACHE_FILE, std::ios::binary);
    if (!out)
    {
        std::cerr << "Cannot open LUT cache for writing: " << LUT_CACHE_FILE << "\n";
        return false;
    }

    // 1) Scalar sizes
    out.write((char *)&ctx.lookup.nbTbins, sizeof(ctx.lookup.nbTbins));
    out.write((char *)&ctx.lookup.nbVPDbins, sizeof(ctx.lookup.nbVPDbins));
    out.write((char *)&ctx.lookup.nbHbins, sizeof(ctx.lookup.nbHbins));

    // 2) Arrays of size ctx.lookup.nbTbins
    out.write((char *)ctx.lookup.LookUp_KmT, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_GammaT, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_VcmaxT, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_JmaxT, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_Rleaf, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_Rstem, ctx.lookup.nbTbins * sizeof(float));

#ifdef WATER
    // 3) WATER-dependent LUTs
    out.write((char *)ctx.lookup.LookUp_SLOPE, ctx.lookup.nbTbins * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_GRADN, ctx.lookup.nbTbins * sizeof(float));

    // INLR: ctx.lookup.nbTbins × ctx.lookup.nbVPDbins
    for (int i = 0; i < ctx.lookup.nbTbins; ++i)
        out.write((char *)ctx.lookup.LookUp_INLR[i], ctx.lookup.nbVPDbins * sizeof(float));

    // Wind LUT
    out.write((char *)ctx.lookup.LookUp_Wind, ctx.lookup.nbHbins * sizeof(float));
#endif

    // 4) Flux LUTs: all size 80000
    out.write((char *)ctx.lookup.LookUp_flux_absorption, 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_flux, 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_ExtinctLW, 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_VPD, 80000 * sizeof(float));
    out.write((char *)ctx.lookup.LookUp_T, 80000 * sizeof(float));

    // 5) Crown site: 2601 ints
    out.write((char *)ctx.lookup.LookUp_Crown_site, 2601 * sizeof(int));

    return out.good();
}

// --- Load LUTs from file -----------------------------------------------------

bool LoadLookUpTablesFromCache(Context &ctx)
{
    std::ifstream in(LUT_CACHE_FILE, std::ios::binary);
    if (!in)
        return false; // cache does not exist

    // 1) Scalar sizes
    in.read((char *)&ctx.lookup.nbTbins, sizeof(ctx.lookup.nbTbins));
    in.read((char *)&ctx.lookup.nbVPDbins, sizeof(ctx.lookup.nbVPDbins));
    in.read((char *)&ctx.lookup.nbHbins, sizeof(ctx.lookup.nbHbins));
    if (!in.good())
        return false;

    // 2) Allocate raw arrays
    ctx.lookup.LookUp_KmT = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_GammaT = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_VcmaxT = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_JmaxT = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_Rleaf = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_Rstem = new float[ctx.lookup.nbTbins];

#ifdef WATER
    ctx.lookup.LookUp_SLOPE = new float[ctx.lookup.nbTbins];
    ctx.lookup.LookUp_GRADN = new float[ctx.lookup.nbTbins];

    ctx.lookup.LookUp_INLR = new float *[ctx.lookup.nbTbins];
    for (int i = 0; i < ctx.lookup.nbTbins; i++)
        ctx.lookup.LookUp_INLR[i] = new float[ctx.lookup.nbVPDbins];

    ctx.lookup.LookUp_Wind = new float[ctx.lookup.nbHbins];
#endif

    ctx.lookup.LookUp_flux_absorption = new float[80000];
    ctx.lookup.LookUp_flux = new float[80000];
    ctx.lookup.LookUp_ExtinctLW = new float[80000];
    ctx.lookup.LookUp_VPD = new float[80000];
    ctx.lookup.LookUp_T = new float[80000];
    // ctx.lookup.LookUp_Crown_site is a static array; no allocation needed.

    // 3) Read back LUT data
    in.read((char *)ctx.lookup.LookUp_KmT, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_GammaT, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_VcmaxT, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_JmaxT, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_Rleaf, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_Rstem, ctx.lookup.nbTbins * sizeof(float));

#ifdef WATER
    in.read((char *)ctx.lookup.LookUp_SLOPE, ctx.lookup.nbTbins * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_GRADN, ctx.lookup.nbTbins * sizeof(float));

    for (int i = 0; i < ctx.lookup.nbTbins; ++i)
        in.read((char *)ctx.lookup.LookUp_INLR[i], ctx.lookup.nbVPDbins * sizeof(float));

    in.read((char *)ctx.lookup.LookUp_Wind, ctx.lookup.nbHbins * sizeof(float));
#endif

    in.read((char *)ctx.lookup.LookUp_flux_absorption, 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_flux, 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_ExtinctLW, 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_VPD, 80000 * sizeof(float));
    in.read((char *)ctx.lookup.LookUp_T, 80000 * sizeof(float));

    in.read((char *)ctx.lookup.LookUp_Crown_site, 2601 * sizeof(int));

    return in.good();
}
