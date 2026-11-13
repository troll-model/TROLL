// ============================================================================
// Parameter Registry System (Temporary System During Refactor)
// ----------------------------------------------------------------------------
// This module replaces the old AssignValueGlobal(). Instead of handling parameters through a long chain of
// if/else conditions, all parameters are now defined, validated, and applied
// through a unified registry table.
//
// WHAT THIS SYSTEM DOES
// ---------------------
// 1. Defines a ParamSpec structure describing:
//       - parameter type (INT, FLOAT, BOOL)
//       - target variable (currently a pointer to a GLOBAL variable)
//       - min, max, default allowed values for validation
//
// 2. Registers all parameters in one place through RegisterParameters().
//    This function fills an unordered_map<string, ParamSpec>, mapping each
//    parameter name (as appears in the input file) to its specification.
//
// 3. AssignParamFromRegistry(name, value) applies a string value to the correct
//    target variable and:
//       - Performs type conversion
//       - Checks ranges
//       - Applies defaults if invalid or missing
//
//
// HOW THIS HELPS THE GLOBALS → CONTEXT REFACTOR
// ----------------------------------------------
// This system has been designed to *simplify and accelerate*
// the migration away from global variables. Specifically:
//
// 1. Centralized parameter definitions
//    - All parameter metadata now lives in one place.
//    - When globals are moved into a Context/Config struct, only the pointer
//      locations in RegisterParameters() need to change.
//
// 2. Parsing is decoupled from storage
//    - AssignParamFromRegistry() does not know (or care) whether parameters
//      live in globals or inside ctx.params.<field>.
//    - Only the `target` pointer inside ParamSpec changes.
//
// 3. Enables removal of globals
//    - After refactoring, forward declarations of globals (currently in
//      param_registry.cpp) will be removed entirely.
//    - The registry will take a Context object and bind ParamSpecs to its
//      fields directly.
//
// 4. Safer refactor through strict validation
//    - All type checking, numeric bound checking, and default fallbacks happen
//      here, ensuring that incorrect refactor steps fail early and loudly.
//
//
// USAGE EXAMPLE
// -------------
// In troll.cpp (or any file reading parameters), instead of:
//
//      AssignValueGlobal(name, value);
//
// you now call:
//
//      AssignParamFromRegistry(name, value);
//
// RegisterParameters() must be invoked *once* at startup (e.g. in main.cpp):
//
//      int main() {
//          RegisterParameters();     // build registry table
//          ReadInputGeneral(ctx);    // parse and assign parameters
//          Initialise(ctx);
//          ...
//      }
//
//
// HOW TO ADD A NEW PARAMETER
// ---------------------------
// Adding a new parameter simple. Only ONE file must be
// modified: param_registry.cpp.
//
// 1. Make sure the underlying variable exists.
//    (currently as a global in troll.hpp; later it will live in Context/Config)
//
// 2. In RegisterParameters(), add exactly one line using the appropriate helper:
//
//        add_int   ("my_param", my_param, min, max, default);
//        add_float ("my_param", my_param, min, max, default);
//        add_bool  ("my_flag",  my_flag,  default);
//
//    Example:
//        add_float("new_param", new_param, 0.0, 5.0, 1.2);
//
// 3. The registry will automatically:
//       - link the name to the variable
//       - convert the string from the input file
//       - enforce valid ranges
//       - assign defaults when missing/invalid
//    No other files need to be touched.
// ============================================================================

#ifndef PARAM_REGISTRY_HPP
#define PARAM_REGISTRY_HPP

#include <string>
#include <unordered_map>
#include "../context.hpp"

// ==========================
// SetParameter declarations
// ==========================

// Numeric template (int, float)
//! Parse a numeric parameter value. Values within [parameter_min, parameter_max], or up to 1%
//! outside it (then clamped), are assigned; other or non-numeric values give parameter_default
//! with a warning. The assigned value is printed unless quiet.
//! @param parameter_name parameter name
//! @param parameter_value value read from the file (text)
//! @param parameter [out] variable receiving the value
//! @param parameter_min minimum allowed value
//! @param parameter_max maximum allowed value
//! @param parameter_default default value used when the value is invalid
//! @param quiet suppress console messages if true
template <typename N>
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  N &parameter,
                  N parameter_min,
                  N parameter_max,
                  N parameter_default,
                  bool quiet);

// String overload
//! Assign a string parameter (parameter_default if empty), printed unless quiet.
//! @param parameter_name parameter name
//! @param parameter_value value read from the file (text)
//! @param parameter [out] variable receiving the value
//! @param parameter_default default value used when the value is invalid
//! @param quiet suppress console messages if true
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  std::string &parameter,
                  std::string parameter_default,
                  bool quiet);

// Boolean overload
//! Assign a boolean parameter given as 0 or 1 (parameter_default otherwise), printed unless quiet.
//! @param parameter_name parameter name
//! @param parameter_value value read from the file (text)
//! @param parameter [out] variable receiving the value
//! @param parameter_min minimum allowed value
//! @param parameter_max maximum allowed value
//! @param parameter_default default value used when the value is invalid
//! @param quiet suppress console messages if true
void SetParameter(std::string &parameter_name,
                  std::string &parameter_value,
                  bool &parameter,
                  bool parameter_min,
                  bool parameter_max,
                  bool parameter_default,
                  bool quiet);

// ==========================
// Registry system
// ==========================

//! Entry of the global-parameter registry: target variable, type, allowed range and default
struct ParamSpec
{
    //! Type of the variable a registry entry assigns
    enum Type
    {
        INT,
        FLOAT,
        BOOL
    } type;          //!< Type of the target variable
    void *target;    //!< Address of the variable in Context that receives the value
    double minv;     //!< Minimum allowed value
    double maxv;     //!< Maximum allowed value
    double def;      //!< Default value

    //! Registry entry for variable tgt of type t, allowed range [mn, mx], default df
    //! @param t    type of the target variable
    //! @param tgt  address of the target variable
    //! @param mn   minimum allowed value
    //! @param mx   maximum allowed value
    //! @param df   default value
    ParamSpec(Type t, void *tgt, double mn, double mx, double df)
        : type(t), target(tgt), minv(mn), maxv(mx), def(df) {}
};

extern std::unordered_map<std::string, ParamSpec> parameter_registry; //!< Global parameters by name, filled by RegisterParameters()

void RegisterParameters(Context &ctx); //!< Register every global parameter of the parameter file, with its range and default

//! Assign the value read for parameter name, validated against its registry entry
//! @param name parameter name
//! @param value value read from the global parameter file (text)
void AssignParamFromRegistry(const std::string &name,
                             const std::string &value);

#endif