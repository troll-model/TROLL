# TROLL regression ("sanity") tests

For developers. `check.sh` runs a **pinned reference version** of TROLL and the
current build on identical inputs, then compares every output file. The default
comparison is bit-exact, which is the right test for refactoring: any change in
any output, however small, fails.

## Usage

```bash
cd build
cmake ..
make sanity          # example inputs, 4 ha from bare ground, 20 days (~1 min)
make sanity1000      # same, 1000 days
make sanity-stand    # established 1 ha stand, 1 year (~2 min)
make check           # sanity + sanity-stand
```

Or call the script directly for other settings:

```bash
sanity/check.sh --troll build/TROLL --case stand --nbiter 730
sanity/check.sh --troll build/TROLL --case example --nbiter 365 --rtol 1e-6
sanity/check.sh --troll build/TROLL --ref <commit>          # other reference
sanity/check.sh --help
```

## Test cases

- **example**: the `./example` configuration (200 × 200 m) starting from bare ground.
  Fast, but during the first years all trees are seedlings, so adult growth,
  crown competition, mortality and treefall are barely exercised.
- **stand**: a 1 ha stand grown for 15 years from the example inputs
  (`fixtures/stand_1ha`, ~730 trees > 10 cm DBH/ha, AGB ~90 t/ha), restarted
  with `-f`/`-w`. This covers the adult-tree code paths and the inventory
  restart. Regenerate or lengthen it with `make_fixture.sh`.

## The reference

Since 4.1.2 the CMake targets compare against **stored reference outputs** in
`golden/` (all outputs under 256 KB of each case, with a README recording the
TROLL version and platform that produced them):

| Directory | Case | Target |
| --- | --- | --- |
| `golden/example_n20` | example, 20 days | `sanity` |
| `golden/stand_n365` | established stand, 365 days | `sanity-stand` |
| `golden/stand_n1000` | established stand, 1000 days | `sanity1000` |

4.1.2 changed results on purpose (Photosyn fix, accurate electron transport), so
an upstream commit can no longer serve as reference. `check.sh --ref <commit>`
still builds and runs any commit, e.g. to compare two branches.

**When results are meant to change** (bug fix, new science): verify the change
separately, run `make golden`, review the differences (`git diff sanity/golden`)
and commit the new references on their own, so the history records when and why
results changed. Never regenerate references to make a failing refactoring pass.

## Tolerances

On the platform that generated the references (Linux x86-64, GCC 13, glibc),
results are bit-identical to them, whatever the number of threads. For
refactoring work, check with no tolerance:

```bash
cmake -DTROLL_SANITY_RTOL_EXAMPLE=0 -DTROLL_SANITY_RTOL_LONG=0 \
      -DTROLL_SANITY_ATOL=0 -DTROLL_SANITY_FILES_LONG=. ..
```

The defaults allow for other compilers and math libraries: lookup tables are
computed with `exp`/`log`, whose last bits differ between libraries, and after a
few months such a difference tips a mortality draw, after which per-tree outputs
follow a different (statistically equivalent) trajectory.

| Variable | Default | Used by |
| --- | --- | --- |
| `TROLL_SANITY_RTOL_EXAMPLE` | `1e-2` | `sanity` (all stored files) |
| `TROLL_SANITY_RTOL_LONG` | `1e-2` | `sanity-stand`, `sanity1000` |
| `TROLL_SANITY_FILES_LONG` | `sumstats\|water_balance\|LAIprofile` | long runs: stand-level aggregates |
| `TROLL_SANITY_ATOL` | `1e-6` | all (ignores negligible values) |

Thread-count independence stays bit-exact: `make threads-check` compares 1 and
4 threads with no tolerance.

## Exact vs tolerance comparison

`compare.py` compares bit-exactly by default. After deliberate numerical
changes (compiler upgrade, float → double, reordering of sums), use
`--rtol`/`--atol`: it reports the largest relative difference per file.
Lines with CPU time and the random seed are ignored.

## Seeds

`--seed N` passes `-rN` to both executables. The reference must support
`-r` (it was introduced together with this test harness), so seed tests need a
`--ref` at or after that commit.

## Benchmarking

`sanity/bench.sh --troll build/TROLL --threads "1 2 4 8"` times the established
stand for each thread count and checks that outputs are identical to the
single-thread run (they must be: see the OpenMP notes in `src/simulation.cpp`).
