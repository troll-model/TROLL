# TROLL sanity tests

For developpers. The test runs a few time steps from TROLL reference version (`main` branch, tag `v4.0`) and compare the outputs bit by bit with current version (`dev` branch), based on the `./example/` configuration. It is insufficient to catch long term trends break but it garantees that initialization and basics functions are not broken.

## Requirements

Same requirements for compiling TROLL, see build [README](../README.md)

## Usage

```bash
cd build
cmake ..
make test # similar to make sanity100
# make sanity1000 for longer run
# make sanity10000 for validating long term trends
```

## Expected outputs

```bash
Comparing test_runs ↔ base_results...
test_0_CHM.txt [OK]
test_0_deathrate.txt [OK]
test_0_death_snapshots.txt [OK]
test_0_death.txt [OK]
test_0_final_pattern.txt [OK]
test_0_final_SWC3D.txt [OK]
test_0_info.txt [OK]
test_0_initial_pattern.txt [OK]
test_0_LAIdynamics.txt [OK]
test_0_LAImature.txt [OK]
test_0_LAIold.txt [OK]
test_0_LAIprofile.txt [OK]
test_0_LAI.txt [OK]
test_0_LAIyoung.txt [OK]
test_0_litterfall.txt [OK]
test_0_phi_root.txt [OK]
test_0_ppfd0.txt [OK]
test_0_sdd.txt [OK]
test_0_site1.txt [OK]
test_0_site2.txt [OK]
test_0_site3.txt [OK]
test_0_site4.txt [OK]
test_0_site5.txt [OK]
test_0_site6.txt [OK]
test_0_soilproperties.txt [OK]
test_0_state_begin.txt [OK]
test_0_state_end.txt [OK]
test_0_state_mid.txt [OK]
test_0_sumstats_species.txt [OK]
test_0_sumstats.txt [OK]
test_0_SWC_begin.txt [OK]
test_0_SWC_end.txt [OK]
test_0_SWC_mid.txt [OK]
test_0_SWP_begin.txt [OK]
test_0_SWP_end.txt [OK]
test_0_SWP_mid.txt [OK]
test_0_vertd.txt [OK]
test_0_water_balance.txt [OK]
test_0_waterfluxes_begin.txt [OK]
test_0_waterfluxes_end.txt [OK]
test_0_waterfluxes_mid.txt [OK]
-----------------------------------
All 41 files match.
[100%] Built target sanity
```

