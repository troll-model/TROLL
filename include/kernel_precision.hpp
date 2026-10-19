#ifndef KERNEL_PRECISION_HPP
#define KERNEL_PRECISION_HPP

//! Guards for the single-precision leaf-flux kernel (Tree::FluxesLeaf, Tree::Photosyn, QUAD).
//!
//! Since 4.1.1 the kernel computes in float throughout, which made the simulation loop
//! ~1.7x faster. One stray double literal (1.0 instead of 1.0f) or double constant silently
//! promotes a whole expression back to double. Code between TROLL_FLOAT_KERNEL_BEGIN and
//! TROLL_FLOAT_KERNEL_END is therefore compiled with -Wdouble-promotion as an error, in
//! every build, whatever the compiler flags. Use the *_F constants from constants.hpp and
//! f-suffixed literals there.

#if defined(__FAST_MATH__)
#error "Do not build TROLL with -ffast-math: the leaf kernel relies on isnan() checks, which fast-math removes."
#endif

#if defined(__GNUC__) // GCC and Clang
#define TROLL_FLOAT_KERNEL_BEGIN \
    _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic error \"-Wdouble-promotion\"")
#define TROLL_FLOAT_KERNEL_END _Pragma("GCC diagnostic pop")
#else
#define TROLL_FLOAT_KERNEL_BEGIN
#define TROLL_FLOAT_KERNEL_END
#endif

#endif // KERNEL_PRECISION_HPP
