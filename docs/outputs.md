# TROLL output files

Each run with output prefix `PREFIX` (option `-o`) and run index `N` (option `-n`, default 0)
writes files named `PREFIX_N_<name>.txt`. All files are tab-separated text with a header line,
except `final_SWC3D.txt`. `iter` is the iteration (timestep) number, starting at 0; with the
default daily timestep, one iteration is one day.

Units below were checked against the code that writes each file (TROLL 4.1.3). Where the code
does not establish a unit unambiguously, the entry says so. Unless stated otherwise, plot-level
quantities are per hectare and fluxes are **per timestep** (the console summary prints some of
them annualised).

Which files are written depends on output options in the global parameter file
(`_OUTPUT_extended` writes most of the files below; `_OUTPUT_inventory` writes the final
inventory and soil water content). `info.txt` and `sumstats.txt` are always written.

## Run information

### `info.txt`
Free-text summary of the run: TROLL version, random seed and its source (`command line`,
`_NONRANDOM` or `time-based`), grid, timing and the "Average computation time" (wall-clock
seconds spent in the simulation loop). Not tab-separated.

## Stand-level time series (one row per iteration)

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

Annual values: multiply the per-timestep fluxes by the number of iterations per year
(`iterperyear` in the global parameters, 365 by default). Example: GPP 0.121 MgC ha⁻¹ per day
× 365 = 44.2 MgC ha⁻¹ yr⁻¹.

### `sumstats_species.txt` — summary per species
Same columns as `sumstats.txt`, with `species` (species name) after `iter`; one row per species
and iteration.

### `death.txt` — mortality
| Column | Meaning | Unit |
| --- | --- | --- |
| `nbdead_n1` | trees that died in the timestep, all sizes | number per timestep (per-hectare scaling not verified) |
| `nbdead_n10` | same, DBH ≥ 10 cm | idem |
| `nbTreefall1` | trees that fell (treefall) | idem |
| `nbTreefall10` | same, DBH ≥ 10 cm | idem |

### `LAIdynamics.txt`
Mean leaf area index over soil cells, at ground level (`LAI`, m² m⁻²), per iteration.

### `LAIyoung.txt`, `LAImature.txt`, `LAIold.txt`
Leaf area density of young, mature and old leaves by height layer: one column per 1 m layer
(`h_0` = ground to 1 m ...), plot mean (m² m⁻² per layer).

### `LAIprofile.txt`
Vertical leaf area profile: `height` (m) and `LAI` (m² m⁻²), several rows per iteration.

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

**Known issue (under review, unchanged in 4.1.3):** `interception` is computed as the sum over
soil cells of interception in mm, multiplied by 0.001 and divided by the number of *sites*,
not of soil cells. It is therefore smaller than the plot-mean interception (in m) by a factor
of (number of soil cells) / (number of sites), e.g. 16 / 10,000 for the 1 ha example.
`throughfall` is computed in mm per soil cell in `simulation.cpp` but commented as m³ in
`output.cpp`; its unit in this file needs checking before use.

### `phi_root.txt`
Water potential in the root zone (MPa), averaged over trees: weighted by abundance
(`abundance_weighted_phi_root`), by abundance of trees ≥ 10 cm DBH (`..._n10`) and by
aboveground biomass (`agb_weighted_phi_root`).

### `ppfd0.txt`
Diagnostics of light at ground level (`variable`, `value`, `residual`), per iteration.

### `litterfall.txt`
Litterfall diagnostics; empty unless the corresponding diagnostic output is enabled.

## Soil (per layer and soil cell)

Soil cells (`dcell`) are squares of `length_dcell` metres (global parameters); columns
`dcell_0 ... dcell_K` are the cells, rows the soil layers.

| File | Content | Unit |
| --- | --- | --- |
| `SWC_begin/mid/end.txt` | soil water content at the beginning, middle and end of the run | m³ per layer and soil cell (as stored in the model) |
| `SWP_begin/mid/end.txt` | soil water potential | MPa |
| `waterfluxes_begin/mid/end.txt` | water fluxes per soil cell (`variable` names the flux) | as in `water_balance.txt` before averaging |
| `soilproperties.txt` | per layer: `depth` (m), saturated / maximum / field-capacity / residual / minimum water content (`sat` and `res` in m³ m⁻³; `max`, `fc`, `min` in m³ per soil cell), `Ksat` (mm s⁻¹), van Genuchten `a_vgm` (MPa) and `m_vgm` | see column |
| `final_SWC3D.txt` | soil water content at the end of the run, **no header**: one row per soil cell, first column the soil-cell index, then one column per layer; read back with option `-w` | m³ per layer and soil cell |

## Individual trees

### `final_pattern.txt`, `initial_pattern.txt` — inventories
One row per living tree at the end (or start) of the run. Readable as an initial forest with
option `-f`. Columns: position (`col`, `row`, `site`), `sp_lab` (species index), `from_Data`
(1 if the tree came from an input inventory), and the tree's traits and state, in the units
of the species file and of the model: `LMA` (g m⁻²), `Nmass` and `Pmass` (mg g⁻¹), `wsg`
(g cm⁻³), `leafarea` (cm²), `tlp` (MPa), `Rdark`, `Vcmax`, `Jmax` (µmol m⁻² s⁻¹),
`dbhmax`, `hmax` (m), `LAImax` (m² m⁻²), DBH and height (m) and further state variables.

### `state_begin/mid/end.txt` — snapshots
One row per living tree at the beginning, middle and end of the run:

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
Time series of the same variables as the snapshots, for six fixed sites (one row per
iteration while a tree occupies the site).

### `death_snapshots.txt`, `deathrate.txt`
Trees that died (`species`, `age` in years, `dbh` and `height` in m), and death-rate
diagnostics (`wsg`, `dbh`, `basal`, `dr`).

## Spatial fields

| File | Content | Unit |
| --- | --- | --- |
| `CHM.txt` | canopy height model at the end of the run: per `site` (`row`, `col`), canopy `height` and spike-free height (m), and `LAI` (m² m⁻²) | m, m² m⁻² |
| `LAI.txt` | full 3D leaf area field: per site (`s`, `row`, `col`) and height `h` (m), cumulated `LAI3D` above that height | m² m⁻² |
| `vertd.txt` | vertical distribution of leaf area (`h`, `freq`) | — |
| `sdd.txt` | seed dispersal distances (`d` in m, `freq`) | — |

## Notes

- Fluxes in `sumstats` are per timestep; the console prints GPP, NPP and litterfall annualised.
- Column names are kept as written by the code, including historical misspellings
  (`transpitation`, `yougLA`), so that existing scripts keep working.
- Entries marked "unit not established" or "known issue" are open questions for the TROLL
  authors; please report findings on the issue tracker.
