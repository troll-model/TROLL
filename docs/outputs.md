# TROLL output files

Every output file, its columns and their units, checked against the code that writes it
(TROLL 4.1.3, rechecked 4.1.6). Where the code does not establish a unit unambiguously, the
entry says so.

## Conventions

| Item | Rule |
| --- | --- |
| File names | `PREFIX_N_<name>.txt`, with `PREFIX` from option `-o` and run index `N` from option `-n` (default 0) |
| Format | tab-separated text with a header line; exceptions: `info.txt` (free text) and `final_SWC3D.txt` (no header) |
| `iter` | iteration (timestep) number, from 0; with the default daily timestep, one iteration is one day |
| Plot-level quantities | per hectare, unless stated otherwise |
| Fluxes | **per timestep**; the console summary prints GPP, NPP and litterfall annualised |
| Annual values | per-timestep flux × iterations per year (`iterperyear`, 365 by default); e.g. GPP 0.121 MgC ha⁻¹ per day × 365 = 44.2 MgC ha⁻¹ yr⁻¹ |
| Which files are written | `info.txt` and `sumstats.txt` always; `_OUTPUT_extended` writes most others; `_OUTPUT_inventory` writes the final inventory and soil water content |
| Column names | kept as written by the code, including historical misspellings (`transpitation`, `yougLA`), so existing scripts keep working |

## Overview

| File | Content | One row per |
| --- | --- | --- |
| `info.txt` | run summary (version, seed, grid, timing) | — (free text) |
| `sumstats.txt` | stand summary: stems, basal area, biomass, carbon fluxes | iteration |
| `sumstats_species.txt` | same, per species | species and iteration |
| `death.txt` | mortality and treefall counts | iteration |
| `LAIdynamics.txt` | mean ground-level LAI | iteration |
| `LAIyoung.txt`, `LAImature.txt`, `LAIold.txt` | leaf area density by leaf age and height | iteration |
| `LAIprofile.txt` | vertical leaf area profile | height and iteration |
| `water_balance.txt` | plot water balance, soil water per layer | iteration |
| `phi_root.txt` | mean root-zone water potential | iteration |
| `ppfd0.txt` | light diagnostics at ground level | iteration |
| `litterfall.txt` | litterfall diagnostics | iteration (when enabled) |
| `SWC_*`, `SWP_*`, `waterfluxes_*`, `soilproperties.txt`, `final_SWC3D.txt` | soil per layer and soil cell | soil layer or soil cell |
| `final_pattern.txt`, `initial_pattern.txt` | tree inventories (readable with `-f`) | tree |
| `state_begin/mid/end.txt` | tree snapshots | tree |
| `site1.txt` ... `site6.txt` | time series of six fixed sites | iteration |
| `death_snapshots.txt`, `deathrate.txt` | dead trees and death-rate diagnostics | dead tree |
| `CHM.txt`, `LAI.txt`, `vertd.txt`, `sdd.txt` | spatial fields and distributions | site, voxel or class |

## Run information

### `info.txt`

| Item | Meaning | Unit |
| --- | --- | --- |
| TROLL version | version that produced the run | — |
| Random seed | seed and its source (`command line`, `_NONRANDOM` or `time-based`) | — |
| Grid and options | plot dimensions and settings of the run | as in the global parameters |
| Average computation time | wall-clock time spent in the simulation loop | s |

## Stand-level time series

### `sumstats.txt` — stand summary

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `sum1` | living trees (all simulated trees, DBH ≥ ~1 cm) | stems ha⁻¹ |
| `sum10` | trees with DBH ≥ 10 cm | stems ha⁻¹ |
| `sum30` | trees with DBH ≥ 30 cm | stems ha⁻¹ |
| `ba` | basal area, all trees | m² ha⁻¹ |
| `ba10` | basal area, trees with DBH ≥ 10 cm | m² ha⁻¹ |
| `agb` | aboveground biomass | kg ha⁻¹ (the console prints t ha⁻¹) |
| `gpp` | gross primary production in the timestep | MgC ha⁻¹ per timestep |
| `npp` | net primary production in the timestep | MgC ha⁻¹ per timestep |
| `rday` | leaf respiration during daytime | MgC ha⁻¹ per timestep |
| `rnight` | leaf respiration at night | MgC ha⁻¹ per timestep |
| `rstem` | stem respiration | MgC ha⁻¹ per timestep |
| `litterfall` | leaf litter | Mg (dry mass) ha⁻¹ per timestep |

### `sumstats_species.txt` — summary per species

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `species` | species name | — |
| `sum1` ... `litterfall` | as in `sumstats.txt`, for the species | as in `sumstats.txt` |

### `death.txt` — mortality

| Column | Meaning | Unit |
| --- | --- | --- |
| `nbdead_n1` | trees that died in the timestep, all sizes | number per timestep (per-hectare scaling not verified) |
| `nbdead_n10` | same, DBH ≥ 10 cm | idem |
| `nbTreefall1` | trees that fell (treefall) | idem |
| `nbTreefall10` | same, DBH ≥ 10 cm | idem |

### `LAIdynamics.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `LAI` | leaf area index at ground level, mean over soil cells | m² m⁻² |

### `LAIyoung.txt`, `LAImature.txt`, `LAIold.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `h_0`, `h_1`, ... | leaf area density of young / mature / old leaves in each 1 m layer (`h_0` = ground to 1 m), plot mean | m² m⁻² per layer |

### `LAIprofile.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `height` | height | m |
| `LAI` | leaf area index at that height | m² m⁻² |

### `water_balance.txt` — plot water balance

| Column | Meaning | Unit |
| --- | --- | --- |
| `precipitation` | rainfall of the day (from the daily climate file) | mm per timestep |
| `interception` | rainfall intercepted by the canopy | **see known issue below** |
| `throughfall` | rainfall reaching the soil | see known issue below |
| `runoff` | surface runoff | m per timestep (plot mean) |
| `leak` | drainage below the deepest soil layer | m per timestep (plot mean) |
| `evaporation` | soil evaporation | m per timestep (plot mean) |
| `transpitation_<l>` | transpiration taken from soil layer *l* (column names keep the historical misspelling) | m per timestep (plot mean) |
| `transpiration1016` | transpiration of trees with DBH 10–16 cm (diagnostic) | m per timestep |
| `SWC_<l>` | soil water content of layer *l*, plot mean | m³ m⁻³ |
| `SWP_<l>` | soil water potential of layer *l*, plot mean | MPa |

Known issues (under review, unchanged since 4.1.3; see [open-questions.md](open-questions.md)):

| Column | Issue |
| --- | --- |
| `interception` | computed as the sum over soil cells of interception in mm, × 0.001, divided by the number of *sites* rather than of soil cells: smaller than the plot-mean interception (m) by (soil cells)/(sites), e.g. 16 / 10,000 for the 1 ha example |
| `throughfall` | computed in mm per soil cell in `simulation.cpp` but commented as m³ in `output.cpp`; check its unit before use |

### `phi_root.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `abundance_weighted_phi_root` | root-zone water potential, mean over trees | MPa |
| `abundance_weighted_phi_root_n10` | same, trees with DBH ≥ 10 cm | MPa |
| `agb_weighted_phi_root` | same, weighted by aboveground biomass | MPa |

### `ppfd0.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration | — |
| `variable` | name of the light diagnostic | — |
| `value` | value of the diagnostic | depends on `variable` |
| `residual` | residual of the diagnostic | depends on `variable` |

### `litterfall.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| — | litterfall diagnostics; the file is empty unless the corresponding diagnostic output is enabled | — |

## Soil (per layer and soil cell)

Soil cells (`dcell`) are squares of `length_dcell` metres (global parameters); in the per-cell
files, columns `dcell_0 ... dcell_K` are the cells and rows the soil layers.

| File | Content | Unit |
| --- | --- | --- |
| `SWC_begin/mid/end.txt` | soil water content at the beginning, middle and end of the run | m³ per layer and soil cell (as stored in the model) |
| `SWP_begin/mid/end.txt` | soil water potential | MPa |
| `waterfluxes_begin/mid/end.txt` | water fluxes per soil cell (`variable` names the flux) | as in `water_balance.txt` before averaging |
| `soilproperties.txt` | per layer: `depth` (m), saturated / maximum / field-capacity / residual / minimum water content (`sat` and `res` in m³ m⁻³; `max`, `fc`, `min` in m³ per soil cell), `Ksat` (mm s⁻¹), van Genuchten `a_vgm` (MPa) and `m_vgm` | see column |
| `final_SWC3D.txt` | soil water content at the end of the run, **no header**: one row per soil cell, first column the soil-cell index, then one column per layer; read back with option `-w` | m³ per layer and soil cell |

## Individual trees

### `final_pattern.txt`, `initial_pattern.txt` — inventories

One row per living tree at the end (or start) of the run; readable as an initial forest with
option `-f`.

| Column | Meaning | Unit |
| --- | --- | --- |
| `col`, `row`, `site` | position of the tree | — |
| `sp_lab` | species index | — |
| `from_Data` | 1 if the tree came from an input inventory | — |
| `LMA` | leaf mass per area | g m⁻² |
| `Nmass`, `Pmass` | leaf nitrogen and phosphorus per dry mass (as in the species file; corrected in 4.1.4) | g g⁻¹ |
| `wsg` | wood specific gravity | g cm⁻³ |
| `leafarea` | leaf area | cm² |
| `tlp` | turgor loss point | MPa |
| `Rdark`, `Vcmax`, `Jmax` | dark respiration, maximum carboxylation and electron transport rates | µmol m⁻² s⁻¹ |
| `dbhmax`, `hmax` | DBH threshold and asymptotic height | m |
| `LAImax` | maximum leaf area index | m² m⁻² |
| DBH, height | tree size | m |
| other columns | further state variables of the tree | as in `state_*` |

### `state_begin/mid/end.txt` — snapshots

One row per living tree at the beginning, middle and end of the run.

| Column | Meaning | Unit |
| --- | --- | --- |
| `site`, `sp_lab` | position and species index | — |
| `height`, `dbh` | tree height, diameter at breast height | m |
| `litter` | leaf litter in the timestep | g (dry mass) |
| `age` | age | years |
| `LA`, `yougLA`, `matureLA`, `oldLA` | total, young, mature, old leaf area (column name keeps the historical misspelling) | m² |
| `CR`, `CD` | crown radius, crown depth | m |
| `GPP`, `NPP`, `Rstem`, `Rnight` | carbon fluxes of the tree in the timestep | gC per timestep |
| `LAI3D_canopy`, `LAI3D_base_crown` | leaf area index above the crown top and above the crown base | m² m⁻² |
| `root_depth` | rooting depth | m |
| `phi_root` | water potential in the root zone | MPa |
| `WSFs`, `WSFns` | water stress factors on stomatal conductance and on non-stomatal processes | 0–1 |
| `transpiration` | transpiration of the tree in the timestep | unit not established from the code; check before use |
| `LAImax`, `LAmax` | maximum leaf area index and leaf area | m² m⁻², m² |
| `root_biomass_<l>` | root biomass in soil layer *l* | unit not established from the code |

### `site1.txt` ... `site6.txt`

| Column | Meaning | Unit |
| --- | --- | --- |
| `iter` | iteration (one row per iteration while a tree occupies the site) | — |
| other columns | as in `state_begin/mid/end.txt`, for the tree at one of six fixed sites | as in `state_*` |

### `death_snapshots.txt`, `deathrate.txt`

| File | Column | Meaning | Unit |
| --- | --- | --- | --- |
| `death_snapshots.txt` | `species` | species of the dead tree | — |
| `death_snapshots.txt` | `age` | age at death | years |
| `death_snapshots.txt` | `dbh`, `height` | size at death | m |
| `deathrate.txt` | `wsg` | wood specific gravity | g cm⁻³ |
| `deathrate.txt` | `dbh` | diameter at breast height | m |
| `deathrate.txt` | `basal` | basal death rate (as named in the code) | unit not established |
| `deathrate.txt` | `dr` | death rate (as named in the code) | unit not established |

## Spatial fields

| File | Content | Unit |
| --- | --- | --- |
| `CHM.txt` | canopy height model at the end of the run: per `site` (`row`, `col`), canopy `height` and spike-free height (m), and `LAI` (m² m⁻²) | m, m² m⁻² |
| `LAI.txt` | full 3D leaf area field: per site (`s`, `row`, `col`) and height `h` (m), cumulated `LAI3D` above that height | m² m⁻² |
| `vertd.txt` | vertical distribution of leaf area (`h`, `freq`) | — |
| `sdd.txt` | seed dispersal distances (`d` in m, `freq`) | — |

## Open points

| Entry | Meaning |
| --- | --- |
| "unit not established" | the code does not determine the unit unambiguously; check before use |
| "known issue" | listed in [open-questions.md](open-questions.md) for the TROLL authors; please report findings on the issue tracker |
