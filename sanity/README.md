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

The reference commit is pinned in `check.sh` (`PINNED_REF`), so the comparison
target never moves silently (the previous script cloned whatever `main` was at
the time). Reference builds and outputs are cached in `sanity/.ref/<commit>/`
and in CI.

When a change is *meant* to alter results (a bug fix, new science), verify the
change separately, then update `PINNED_REF` to the new commit in a dedicated
commit, so the history records when and why reference outputs changed.

## Tolerances since 4.1.1

4.1.1 computes the leaf kernel in single precision, a deliberate numerical change. The
pinned reference (4.0 / dev@4de3273) is therefore matched within tolerances measured for
that change, set as CMake cache variables:

| Variable | Default | Used by |
| --- | --- | --- |
| `TROLL_SANITY_RTOL_EXAMPLE` | `1e-2` | `sanity`, `sanity1000` |
| `TROLL_SANITY_RTOL_STAND` | `1e-3` | `sanity-stand` |
| `TROLL_SANITY_ATOL` | `1e-6` | all (ignores negligible values, e.g. deep-soil transpiration ~1e-12) |

Once a 4.1.1 commit is published, move `PINNED_REF` to it and set the tolerances to 0 to
restore bit-exact checking. Thread-count independence stays bit-exact: `make threads-check`
compares 1 and 4 threads with no tolerance.

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
