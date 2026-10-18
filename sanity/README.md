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

## Exact vs tolerance comparison

`compare.py` compares bit-exactly by default. After deliberate numerical
changes (compiler upgrade, float → double, reordering of sums), use
`--rtol`/`--atol`: it reports the largest relative difference per file.
Lines with CPU time and the random seed are ignored.

## Seeds

`--seed N` passes `-rN` to both executables. The reference must support
`-r` (it was introduced together with this test harness), so seed tests need a
`--ref` at or after that commit.
