#ifndef QUAD_HPP
#define QUAD_HPP

#include <cmath>
#include <iostream>
#include "kernel_precision.hpp"

TROLL_FLOAT_KERNEL_BEGIN // part of the single-precision leaf kernel

//! Root of A x^2 + B x + C = 0: (-B + ROOT * sqrt(B^2 - 4AC)) / (2A), ROOT = 1 or -1
//! (for A > 0, ROOT = 1 is the larger root, ROOT = -1 the smaller one).
//! Used by Photosyn for the intercellular CO2 of the coupled photosynthesis-conductance
//! solution (larger root), where this textbook form is accurate in single precision
//! (max relative error 1.3e-6 measured over 2.3e7 calls). For the electron-transport rate,
//! use ElectronTransportRate() below.
//! Inline (4.1.2) so it is inlined into Photosyn even without link-time optimisation.
//! @param A     coefficient of x^2
//! @param B     coefficient of x
//! @param C     constant term
//! @param ROOT     1 or -1, sign in front of the square root
//! @return the selected root; NaN if B^2 - 4AC < 0 (a warning is printed)
inline float QUAD(float A, float B, float C, int ROOT)
{
    float quad = 0.0f;
    float discriminant = B * B - 4.0f * A * C;
    if (discriminant < 0.0f)
    {
        std::cerr << "IMAGINARY ROOTS IN QUADRATIC" << std::endl; // result is NaN below; callers check isnan
        quad = 0.0f;
    }
    if (A == 0.0f)
    {
        if (B == 0.0f)
            quad = 0.0f;
        else
            quad = -C / B;
    }
    else
        quad = (-B + float(ROOT) * std::sqrt(discriminant)) / (2.0f * A);
    return quad;
}

//! Electron-transport rate J (umol m-2 s-1): smaller root of the non-rectangular hyperbola
//!   theta J^2 - (I + Jmax) J + I Jmax = 0,
//! written as J = 2 I Jmax / ((I + Jmax) + sqrt((I + Jmax)^2 - 4 theta I Jmax)).
//! Algebraically identical to (I + Jmax - sqrt(...)) / (2 theta), but without its
//! cancellation at low light (I << Jmax), which in single precision lost up to three
//! significant digits (max relative error 1.9e-3; now ~2e-7). Same cost: one sqrt, one
//! division, no branches. Requires 0 < theta <= 1, I >= 0, Jmax >= 0 (discriminant >= 0).
//! GMD Part 1, Eq. 25.
//! @param theta  curvature of the light response (dimensionless)
//! @param I      light usable for electron transport, alpha x PPFD (umol e- m-2 s-1)
//! @param Jmax   maximum electron transport rate at leaf temperature (umol e- m-2 s-1)
//! @return J (umol e- m-2 s-1)
inline float ElectronTransportRate(float theta, float I, float Jmax)
{
    const float s = I + Jmax;
    const float denominator = s + std::sqrt(s * s - 4.0f * theta * I * Jmax);
    return denominator > 0.0f ? 2.0f * I * Jmax / denominator : 0.0f;
}

TROLL_FLOAT_KERNEL_END

#endif // QUAD_HPP
