#// ============================================================================
#// Parameter Registry Implementation
#// ----------------------------------------------------------------------------
#// This file implements the parameter registry introduced to replace the old
#// AssignValueGlobal() system previously embedded inside troll.cpp.
#
#// GOALS
#// -----
#// - Centralize all parameter bindings and metadata (min/max/default/type)
#// - Remove long if/else chains
#// - Allow automated validation and conversion
#// - Prepare the codebase for migration away from global variables
#
#// CURRENT STATE (TEMPORARY)
#// -------------------------
#// The registry still points to GLOBAL variables declared in troll.cpp/
#// troll.hpp. This is why we temporarily forward-declare them below.
#
#// After the upcoming refactor:
#//   - All globals will live inside Context or Config structs
#//   - The registry will store pointers to those struct fields
#//   - All the forward declarations here will be removed
#
#// HOW IT WORKS
#// ------------
#// 1. RegisterParameters() builds a map<string, ParamSpec>.
#// 2. Each ParamSpec contains:
#//       * parameter type (INT/FLOAT/BOOL)
#//       * raw pointer to target variable
#//       * min value
#//       * max value
#//       * default value
#// 3. AssignParamFromRegistry() looks up an input parameter and uses the
#//    appropriate SetParameter() overload to assign validated values.
#
#// NOTE:
#// ----------
#// - This module does NOT read files.
#// - It only applies validated parameter values.
#// - ReadInputGeneral() in troll.cpp performs the actual scanning.
#
#// ============================================================================

#include "param_registry.hpp"
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <climits>
#include "../../include/troll_defines.hpp"

// NOTE (Temporary Forward Declarations):
// -------------------------------------
// These extern declarations exist ONLY because the current codebase still
// relies on a large number of global variables defined in troll.hpp.
//
// Once TROLL is properly refactored to:
//   (1) move all global parameters into a Context or Config struct,
//   (2) provide proper header separation,
//   (3) eliminate direct global access from unrelated modules,
//
// THIS ENTIRE FORWARD-DECLARATION SECTION WILL BE REMOVED.
//
// The parameter registry will then receive a reference to a parameter
// container (e.g., ctx.params), instead of touching globals directly.
//
// For now, these declarations prevent circular include dependencies.

// All globals formerly declared here have been migrated to Context sub-structs:
//   cols, rows, HEIGHT, length_dcell, NV, NH, leafdem_resolution → ctx.grid
//   p_nonvert, SWtoPPFD, klight, absorptance_leaves, theta, phi, g1, g0,
//   pheno_a0, pheno_b0, pheno_delta, vC, DBH0, H0, CR_min, CR_a, CR_b,
//   CD_a, CD_b, CD0, dens, fallocwood, falloccanopy, Cseedrain, nbs0,
//   p_tfsecondary, hurt_decay, m, m1, Cair, PRESS → ctx.params
//   shape_crown, crown_gap_fraction, extent_visual → ctx.crown
//   sigma/corr/cov intraspecific variation → ctx.intra
//   ModelOptions → ctx.opt

// ============================================================================
//   SetParameter numeric template
// ============================================================================
template <typename N>
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  N &parameter,
                  N parameter_min,
                  N parameter_max,
                  N parameter_default,
                  bool quiet)
{
    std::istringstream iss(parameter_value);
    N numeric;
    iss >> numeric;

    bool isnumeric = iss.eof() && !iss.fail();

    if (isnumeric)
    {
        if (numeric >= parameter_min * 0.99 && numeric <= parameter_max * 1.01)
        {
            if (numeric < parameter_min)
                parameter = parameter_min;
            else if (numeric > parameter_max)
                parameter = parameter_max;
            else
                parameter = numeric;

            if (!quiet)
                std::cout << parameter_name << ": " << parameter << std::endl;
        }
        else
        {
            parameter = parameter_default;
            if (!quiet)
                std::cout << "Warning: '" << parameter_name << "' value "
                          << numeric << " outside range (" << parameter_min
                          << ", " << parameter_max << "). Default: "
                          << parameter_default << std::endl;
        }
    }
    else
    {
        parameter = parameter_default;
        if (!quiet)
            std::cout << "Warning: '" << parameter_name
                      << "' invalid value '" << parameter_value
                      << "'. Default: " << parameter_default << std::endl;
    }
}

// explicit instantiations
//! Explicit instantiation of SetParameter for int parameters
template void SetParameter<int>(std::string &, std::string &, int &, int, int, int, bool);
//! Explicit instantiation of SetParameter for float parameters
template void SetParameter<float>(std::string &, std::string &, float &, float, float, float, bool);
//! Explicit instantiation of SetParameter for unsigned short parameters
template void SetParameter<unsigned short>(std::string &, std::string &, unsigned short &, unsigned short, unsigned short, unsigned short, bool);

// ============================================================================
//   SetParameter string overload
// ============================================================================
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  std::string &parameter,
                  std::string parameter_default,
                  bool quiet)
{
    if (!parameter_value.empty())
    {
        parameter = parameter_value;
        if (!quiet)
            std::cout << parameter_name << ": " << parameter << std::endl;
    }
    else
    {
        parameter = parameter_default;
        std::cout << "Warning: empty string for '" << parameter_name
                  << "'. Default: '" << parameter_default << "'" << std::endl;
    }
}

// ============================================================================
//   SetParameter boolean overload
// ============================================================================
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  bool &parameter,
                  bool parameter_min,
                  bool parameter_max,
                  bool parameter_default,
                  bool quiet)
{
    if (parameter_value == "1" || parameter_value == "true" || parameter_value == "True")
        parameter = true;
    else if (parameter_value == "0" || parameter_value == "false" || parameter_value == "False")
        parameter = false;
    else
    {
        parameter = parameter_default;
        if (!quiet)
            std::cout << "Warning: '" << parameter_name
                      << "' invalid bool '" << parameter_value
                      << "'. Default: " << parameter_default << std::endl;
        return;
    }

    if (parameter < parameter_min)
        parameter = parameter_min;
    if (parameter > parameter_max)
        parameter = parameter_max;

    if (!quiet)
        std::cout << parameter_name << ": " << parameter << std::endl;
}

// ============================================================================
//   Global registry map
// ============================================================================
std::unordered_map<std::string, ParamSpec> parameter_registry;

// ============================================================================
//   RegisterParameters()
// ============================================================================
void RegisterParameters(Context &ctx)
{
    parameter_registry.clear();

    // helper lambdas
    auto add_int = [&](const std::string &name, int &var,
                       int minv, int maxv, int defv)
    {
        parameter_registry.emplace(
            name, ParamSpec(ParamSpec::INT, &var, minv, maxv, defv));
    };

    auto add_float = [&](const std::string &name, float &var,
                         float minv, float maxv, float defv)
    {
        parameter_registry.emplace(
            name, ParamSpec(ParamSpec::FLOAT, &var, minv, maxv, defv));
    };

    auto add_bool = [&](const std::string &name, bool &var,
                        bool minv, bool maxv, bool defv)
    {
        parameter_registry.emplace(
            name, ParamSpec(ParamSpec::BOOL, &var, minv, maxv, defv));
    };

    // ============================================
    //  FULL registry table from your AssignValueGlobal()
    // ============================================

    add_int("cols", ctx.grid.cols, 0, INT_MAX, 400);
    add_int("rows", ctx.grid.rows, 0, INT_MAX, 400);
    add_int("HEIGHT", ctx.grid.HEIGHT, 0, 150, 70);
    add_int("length_dcell", ctx.grid.length_dcell, 0, INT_MAX, 25);
    add_int("nbiter", ctx.time.nbiter, 0, INT_MAX, 6000);

    add_float("NV", ctx.grid.NV, 0.0f, float(INT_MAX), 1.0f);
    add_float("NH", ctx.grid.NH, 0.0f, float(INT_MAX), 1.0f);

    add_int("nbout", ctx.time.nbout, 0, INT_MAX, 4);
    add_float("p_nonvert", ctx.params.p_nonvert, 0.0f, 1.0f, 0.05f);
    add_float("SWtoPPFD", ctx.params.SWtoPPFD, 0.0f, 5.0f, 2.29f);
    add_float("klight", ctx.params.klight, 0.0f, 1.0f, 0.5f);

    add_float("absorptance_leaves", ctx.params.absorptance_leaves, 0.0f, 1.0f, 0.9f);
    add_float("theta", ctx.params.theta, 0.0f, 10.0f, 0.7f);
    add_float("phi", ctx.params.phi, 0.0f, 1.0f, 0.06f);
    add_float("g1", ctx.params.g1, 0.0f, 1000.0f, 3.77f);

#ifdef G0
    add_float("g0", ctx.params.g0, 0.0f, 30.0f, 5.0f);
#endif

#ifdef PHENO_DROUGHT
    add_float("pheno_a0", ctx.params.pheno_a0, 0.0f, 1.0f, 0.5f);
    add_float("pheno_b0", ctx.params.pheno_b0, 0.0f, 1.0f, 0.5f);
    add_float("pheno_delta", ctx.params.pheno_delta, 0.0f, 1.0f, 0.1f);
#endif

    add_float("vC", ctx.params.vC, 0.0f, 1.0f, 0.05f);
    add_float("DBH0", ctx.params.DBH0, 0.0f, 2.5f, 0.005f);
    add_float("H0", ctx.params.H0, 0.0f, 100.0f, 0.95f);
    add_float("CR_min", ctx.params.CR_min, 0.0f, 50.0f, 0.2f);
    add_float("CR_a", ctx.params.CR_a, 0.0f, 5.0f, 2.13f);
    add_float("CR_b", ctx.params.CR_b, 0.0f, 50.0f, 0.63f);
    add_float("CD_a", ctx.params.CD_a, 0.0f, 0.5f, 0.0f);
    add_float("CD_b", ctx.params.CD_b, 0.0f, 1.0f, 0.2f);
    add_float("CD0", ctx.params.CD0, 0.0f, 50.0f, 0.1f);
    add_float("shape_crown", ctx.crown.shape_crown, 0.0f, 1.0f, 1.0f);
    add_float("dens", ctx.params.dens, 0.0f, 10.0f, 1.0f);
    add_float("fallocwood", ctx.params.fallocwood, 0.0f, 1.0f, 0.35f);
    add_float("falloccanopy", ctx.params.falloccanopy, 0.0f, 1.0f, 0.25f);
    add_float("Cseedrain", ctx.params.Cseedrain, 0.0f, 1000000.0f, 50000.0f);

    add_float("nbs0", ctx.params.nbs0, 0.0f, 10000.0f, 10.0f);

    add_float("sigma_height", ctx.intra.sigma_height, 0.0f, 1.0f, 0.19f);
    add_float("sigma_CR", ctx.intra.sigma_CR, 0.0f, 1.0f, 0.29f);
    add_float("sigma_CD", ctx.intra.sigma_CD, 0.0f, 1.0f, 0.0f);
    add_float("sigma_P", ctx.intra.sigma_P, 0.0f, 1.0f, 0.24f);
    add_float("sigma_N", ctx.intra.sigma_N, 0.0f, 1.0f, 0.12f);
    add_float("sigma_LMA", ctx.intra.sigma_LMA, 0.0f, 1.0f, 0.24f);
    add_float("sigma_wsg", ctx.intra.sigma_wsg, 0.0f, 0.5f, 0.06f);
    add_float("sigma_dbhmax", ctx.intra.sigma_dbhmax, 0.0f, 1.0f, 0.05f);

    add_float("sigma_leafarea", ctx.intra.sigma_leafarea, 0.0f, 1.0f, 0.05f);
    add_float("sigma_tlp", ctx.intra.sigma_tlp, 0.0f, 1.0f, 0.05f);

    add_float("corr_CR_height", ctx.intra.corr_CR_height, -1.0f, 1.0f, 0.0f);
    add_float("corr_N_P", ctx.intra.corr_N_P, -1.0f, 1.0f, 0.65f);
    add_float("corr_N_LMA", ctx.intra.corr_N_LMA, -1.0f, 1.0f, -0.43f);
    add_float("corr_P_LMA", ctx.intra.corr_P_LMA, -1.0f, 1.0f, -0.39f);

    add_int("leafdem_resolution", ctx.grid.leafdem_resolution, 0, INT_MAX, 30);

    add_float("p_tfsecondary", ctx.params.p_tfsecondary, 0.0f, 1.0f, 1.0f);
    add_float("hurt_decay", ctx.params.hurt_decay, 0.0f, 1.0f, 0.0f);
    add_float("crown_gap_fraction", ctx.crown.crown_gap_fraction, 0.0f, 1.0f, 0.0f);
    add_float("m", ctx.params.m, 0.0f, 1.0f, 0.013f);
    add_float("m1", ctx.params.m1, 0.0f, 1.0f, 0.013f);

    add_float("Cair", ctx.params.Cair, 0.0f, 1000000.0f, 400.0f);

#ifdef WATER
    add_float("PRESS", ctx.params.PRESS, 10.0f, 110.0f, 101.0f);
#endif

    add_bool("_LL_parameterization", ctx.opt._LL_parameterization, false, true, true);
    add_int("_LA_regulation", ctx.opt._LA_regulation, 0, 2, 2);
    add_bool("_sapwood", ctx.opt._sapwood, false, true, true);
    add_bool("_seedsadditional", ctx.opt._seedsadditional, false, true, false);

    add_int("_SOIL_LAYER_WEIGHT", ctx.opt._SOIL_LAYER_WEIGHT, 0, 2, 2);
    add_int("_WATER_RETENTION_CURVE", ctx.opt._WATER_RETENTION_CURVE, 0, 1, 0);

    add_bool("_NONRANDOM", ctx.opt._NONRANDOM, false, true, true);

    add_bool("_GPPcrown", ctx.opt._GPPcrown, false, true, false);
    add_bool("_BASICTREEFALL", ctx.opt._BASICTREEFALL, false, true, true);
    add_bool("_SEEDTRADEOFF", ctx.opt._SEEDTRADEOFF, false, true, false);
    add_bool("_NDD", ctx.opt._NDD, false, true, false);
    add_bool("_CROWN_MM", ctx.opt._CROWN_MM, false, true, false);
    add_bool("_OUTPUT_extended", ctx.opt._OUTPUT_extended, false, true, false);
    add_bool("_OUTPUT_inventory", ctx.opt._OUTPUT_inventory, false, true, false);

    add_int("extent_visual", ctx.crown.extent_visual, 0, INT_MAX, 0);
}

// ============================================================================
//   AssignParamFromRegistry()
// ============================================================================
void AssignParamFromRegistry(const std::string &name,
                             const std::string &value)
{
    auto it = parameter_registry.find(name);
    if (it == parameter_registry.end())
    {
        std::cerr << "Warning: Unknown parameter '" << name << "'\n";
        return;
    }

    std::string pname = name;
    std::string pvalue = value;

    ParamSpec &p = it->second;
    bool quiet = true;

    switch (p.type)
    {
    case ParamSpec::INT:
    {
        int &ref = *static_cast<int *>(p.target);
        SetParameter(pname, pvalue, ref,
                     (int)p.minv,
                     (int)p.maxv,
                     (int)p.def,
                     quiet);
        break;
    }
    case ParamSpec::FLOAT:
    {
        float &ref = *static_cast<float *>(p.target);
        SetParameter(pname, pvalue, ref,
                     (float)p.minv,
                     (float)p.maxv,
                     (float)p.def,
                     quiet);
        break;
    }
    case ParamSpec::BOOL:
    {
        bool &ref = *static_cast<bool *>(p.target);
        SetParameter(pname, pvalue, ref,
                     (bool)p.minv,
                     (bool)p.maxv,
                     (bool)p.def,
                     quiet);
        break;
    }
    }
}