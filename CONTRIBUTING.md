# Contributing to TROLL

## Build and test

```bash
cd build
cmake .. -DTROLL_WARNINGS_AS_ERRORS=ON
make            # TROLL executable
make check      # unit tests + regression checks (about 2 minutes)
make threads-check   # results must be identical with 1 and 4 threads
make docs-check      # documentation: coverage report, fails on documentation errors
```

Every pull request runs these in CI (`.github/workflows/`).

## Changes that must not change results

Refactoring, performance work and documentation must leave outputs **bit-identical**. On
Linux x86-64 with GCC 13, check against the stored reference outputs with no tolerance:

```bash
sanity/check.sh --troll build/TROLL --case stand --nbiter 365 --golden golden/stand_n365 --rtol 0 --atol 0
```

## Changes that are meant to change results

Bug fixes and new science change outputs. Then:

1. verify the change on its own (what changes, by how much, why);
2. regenerate the stored references with `make golden`, review `git diff sanity/golden`;
3. commit the new references in their own commit, and describe the change in `NEWS.md`;
4. for changes to model behaviour, validate long-run statistics over several seeds.

Never regenerate the references to make a failing refactoring pass.

## Leaf kernel precision rules

`FluxesLeaf`, `Photosyn`, `QUAD` and `ElectronTransportRate` compute in single precision and are
compiled with `-Wdouble-promotion` as an error (`include/kernel_precision.hpp`). Use the `*_F`
constants of `constants.hpp` and `f`-suffixed literals there. Do not build with `-ffast-math`.

## Documentation rules

- Document every new function (`//!` before it, `@param` for each parameter) and data member (`//!<`).
- Give the **unit** of every physical quantity, and the **source** of every equation or constant
  (GMD Part 1 equation number, or the literature reference).
- Input and output files: update `docs/inputs.md` and `docs/outputs.md`.
- Discrepancies you cannot resolve go to `docs/open-questions.md`, not into silent fixes.
