// Unit tests for the quadratic solvers of the leaf kernel (include/quad.hpp).
// Run with `make unit-tests` (also part of `make check`).
#include <cmath>
#include <cstdio>
#include <random>
#include "quad.hpp"

static int failures = 0;
#define CHECK(cond, ...)                                     \
    do                                                       \
    {                                                        \
        if (!(cond))                                         \
        {                                                    \
            failures++;                                      \
            std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
            std::printf(__VA_ARGS__);                        \
            std::printf("\n");                               \
        }                                                    \
    } while (0)

// (-B + ROOT * sqrt(B^2 - 4AC)) / 2A in extended precision, without cancellation
static long double quad_ref(long double A, long double B, long double C, int ROOT)
{
    long double s = std::sqrt(B * B - 4.0L * A * C);
    long double q = -0.5L * (B + std::copysign(s, B));
    bool same_sign = (ROOT > 0) == !std::signbit(B);
    return same_sign ? C / q : q / A;
}

int main()
{
    // ---- QUAD: root selection (-B + ROOT sqrt(d)) / 2A, both signs of A and B ----
    CHECK(QUAD(1.0f, -3.0f, 2.0f, 1) == 2.0f, "x^2-3x+2, ROOT=1");
    CHECK(QUAD(1.0f, -3.0f, 2.0f, -1) == 1.0f, "x^2-3x+2, ROOT=-1");
    CHECK(QUAD(1.0f, 3.0f, 2.0f, 1) == -1.0f, "x^2+3x+2, ROOT=1");
    CHECK(QUAD(1.0f, 3.0f, 2.0f, -1) == -2.0f, "x^2+3x+2, ROOT=-1");
    CHECK(QUAD(-1.0f, 3.0f, -2.0f, 1) == 1.0f, "-x^2+3x-2, ROOT=1");
    CHECK(QUAD(-1.0f, 3.0f, -2.0f, -1) == 2.0f, "-x^2+3x-2, ROOT=-1");
    // ---- QUAD: edge cases ----
    CHECK(QUAD(0.0f, 2.0f, -4.0f, 1) == 2.0f, "linear 2x-4 -> 2");
    CHECK(QUAD(0.0f, 0.0f, 5.0f, 1) == 0.0f, "A = B = 0 -> 0");
    CHECK(QUAD(1.0f, 0.0f, 0.0f, 1) == 0.0f, "x^2 = 0 -> 0");
    std::fprintf(stderr, "(expected message below: negative discriminant)\n");
    CHECK(std::isnan(QUAD(1.0f, 0.0f, 1.0f, 1)), "negative discriminant -> NaN (callers test isnan)");

    // ---- QUAD: accuracy of the larger root when no cancellation occurs (Ci in Photosyn) ----
    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> mag(-4.0f, 4.0f), u(0.0f, 1.0f);
    double worst_quad = 0.0;
    long n_quad = 0;
    for (int i = 0; i < 1000000; i++)
    {
        float A = std::pow(10.0f, mag(rng)) * (u(rng) < 0.5f ? -1.0f : 1.0f);
        float B = std::pow(10.0f, mag(rng)) * (u(rng) < 0.5f ? -1.0f : 1.0f);
        float C = std::pow(10.0f, mag(rng)) * (u(rng) < 0.5f ? -1.0f : 1.0f);
        int ROOT = (B < 0.0f) ? 1 : -1; // -B and ROOT * sqrt(d) have the same sign: no cancellation
        long double d = (long double)B * B - 4.0L * A * C;
        if (d <= 0.0L || d < 0.25L * (long double)B * B) continue; // well-conditioned discriminant (B^2 - 4AC not itself cancelling)
        long double ref = quad_ref(A, B, C, ROOT);
        double err = std::fabs((double)((QUAD(A, B, C, ROOT) - ref) / ref));
        worst_quad = std::fmax(worst_quad, err);
        n_quad++;
    }
    std::printf("QUAD: %ld well-conditioned cases, max rel error %.3g\n", n_quad, worst_quad);
    CHECK(worst_quad < 2e-6, "QUAD max rel error %.3g >= 2e-6", worst_quad);

    // ---- ElectronTransportRate: accuracy over light, Jmax and theta ----
    double worst_j = 0.0;
    for (float theta : {0.05f, 0.3f, 0.7f, 0.9f, 0.99f})
        for (float Jmax = 1.0f; Jmax < 500.0f; Jmax *= 1.3f)
            for (float I = 1e-3f; I < 3000.0f; I *= 1.07f)
            {
                long double ref = quad_ref(theta, -((long double)I + Jmax), (long double)I * Jmax, -1);
                double err = std::fabs((double)((ElectronTransportRate(theta, I, Jmax) - ref) / ref));
                worst_j = std::fmax(worst_j, err);
            }
    std::printf("ElectronTransportRate: max rel error %.3g (theta 0.05-0.99, Jmax 1-500, I 1e-3-3000)\n", worst_j);
    CHECK(worst_j < 1e-6, "ElectronTransportRate max rel error %.3g >= 1e-6", worst_j);

    // The textbook form it replaces loses accuracy at low light: document the difference
    {
        float theta = 0.7f, I = 0.05f, Jmax = 150.0f;
        long double ref = quad_ref(theta, -((long double)I + Jmax), (long double)I * Jmax, -1);
        double err_old = std::fabs((double)((QUAD(theta, -(I + Jmax), I * Jmax, -1) - ref) / ref));
        double err_new = std::fabs((double)((ElectronTransportRate(theta, I, Jmax) - ref) / ref));
        std::printf("low light (I=%.2f, Jmax=%.0f): textbook form rel error %.3g, ElectronTransportRate %.3g\n", I, Jmax, err_old, err_new);
        CHECK(err_new < err_old || err_old < 1e-7, "stable form should not be less accurate");
    }

    // ---- ElectronTransportRate: edge cases ----
    CHECK(ElectronTransportRate(0.7f, 0.0f, 150.0f) == 0.0f, "no light -> J = 0");
    CHECK(ElectronTransportRate(0.7f, 500.0f, 0.0f) == 0.0f, "no Jmax -> J = 0");
    CHECK(ElectronTransportRate(0.7f, 0.0f, 0.0f) == 0.0f, "I = Jmax = 0 -> 0 (no 0/0)");
    CHECK(std::fabs(ElectronTransportRate(1.0f, 100.0f, 150.0f) - 100.0f) < 1e-4f, "theta = 1 -> min(I, Jmax)");
    {
        float J = ElectronTransportRate(1e-6f, 100.0f, 150.0f), lim = 100.0f * 150.0f / 250.0f;
        CHECK(std::fabs(J - lim) / lim < 1e-4f, "theta -> 0 -> I Jmax / (I + Jmax), got %g vs %g", J, lim);
    }
    {
        float theta = 0.7f, I = 800.0f, Jmax = 120.0f, J = ElectronTransportRate(theta, I, Jmax);
        float residual = theta * J * J - (I + Jmax) * J + I * Jmax;
        CHECK(std::fabs(residual) / (I * Jmax) < 1e-5f, "J solves theta J^2 - (I+Jmax) J + I Jmax = 0 (residual %g)", residual);
        CHECK(J > 0.0f && J <= std::fmin(I, Jmax), "J is the smaller root, within [0, min(I, Jmax)]");
    }

    std::printf(failures ? "%d FAILED\n" : "All quadratic-solver tests passed\n", failures);
    return failures ? 1 : 0;
}
