# TROLL input files

TROLL reads five input files, given on the command line, plus two optional ones:

| Option | File | Content |
| --- | --- | --- |
| `-i` | global parameters | grid, simulation length, physiology and allometry parameters, options |
| `-s` | species parameters | one row per species: traits and allometric parameters |
| `-m` | daily climate | one row per day |
| `-d` | half-hourly climate | one row per daytime half-hour |
| `-p` | soil | one row per soil layer |
| `-f` (optional) | forest inventory | initial forest, e.g. the `final_pattern.txt` output of a previous run |
| `-w` (optional) | soil water content | initial soil water, e.g. the `final_SWC3D.txt` output of a previous run |

Units were established by cross-checking three sources: the code that reads and uses each
value, the GMD papers (Maréchaux et al. 2025, Part 1; Schmitt et al. 2025, Part 2), and the
magnitudes in `example/`. Where the sources disagree, the unit the **code** expects is given and
the disagreement is listed in [open-questions.md](open-questions.md).

## Species file (`-s`)

| Column | Meaning | Unit | Example range | Notes |
| --- | --- | --- | --- | --- |
| `s_name` | species name (`Genus_species`) | — | — | |
| `s_LMA` | leaf mass per area | g m⁻² | 55–177 | |
| `s_Nmass` | leaf nitrogen per dry mass | **g g⁻¹** | 0.014–0.027 | Part 1 states mg g⁻¹; the code multiplies by 1000 to get mg g⁻¹ for the Vcmax relationship |
| `s_Pmass` | leaf phosphorus per dry mass | **g g⁻¹** | 0.00038–0.00093 | as `s_Nmass` |
| `s_wsg` | wood specific gravity | g cm⁻³ | 0.43–0.83 | |
| `s_dbhmax` | DBH threshold (dbh_thres) | m | 0.17–1.03 | Part 2 text says cm, its Fig. 5 metres; the file uses metres |
| `s_hmax` | asymptotic height h_lim (Eq. 16) | m | 41–53 | |
| `s_ah` | height–DBH parameter a_h (Eq. 16) | m | 0.19–0.39 | |
| `s_seedmass` | seed mass | not established | 1 | used with the `_SEEDTRADEOFF` option |
| `s_regionalfreq` | relative frequency in the external seed rain | — (fraction) | 0.0092 | uniform in the example |
| `s_tlp` | leaf water potential at turgor loss point | MPa | −2.28 to −1.50 | |
| `s_leafarea` | leaf area | cm² | 3–494 | |

## Daily climate file (`-m`)

| Column | Meaning | Unit |
| --- | --- | --- |
| `NightTemperature` | mean night-time air temperature | °C |
| `Rainfall` | daily rainfall (treated as one event per day) | mm |

## Half-hourly climate file (`-d`)

One row per daytime half-hour (net radiation > 0), all days in sequence.

| Column | Meaning | Unit |
| --- | --- | --- |
| `DayJulian` | day of the year | — |
| `time_numeric` | time of day | h |
| `Temp` | air temperature above the canopy | °C |
| `Snet` | incoming shortwave radiation above the canopy | W m⁻² |
| `VPD` | vapour pressure deficit above the canopy | kPa |
| `WS` | wind speed at the meteorological station | m s⁻¹ |

## Soil file (`-p`)

One row per layer, from the surface down. Hydraulic properties are derived from these with
pedotransfer functions (Hodnett & Tomasella 2002; Cosby et al. 1984 for saturated conductivity).

| Column | Meaning | Unit |
| --- | --- | --- |
| `layer_thickness` | layer thickness | m |
| `proportion_Silt`, `proportion_Clay`, `proportion_Sand` | texture | % |
| `SOC` | soil organic carbon | % (converted to g kg⁻¹ in the code) |
| `DBD` | dry bulk density | g cm⁻³ |
| `pH` | pH | — |
| `CEC` | cation exchange capacity | not stated in the code (as expected by the Hodnett & Tomasella 2002 pedotransfer functions) |

## Global parameters (`-i`)

Generated from `example/global_inputs.txt` (TROLL 4.1.4): name, example value, unit, and the
description given in the file. Parameters beginning with `_` are options.

| Parameter | Example | Unit | Description |
| --- | --- | --- | --- |
| `cols` | 200 | cells (1 cell = 1/NH m) | nb of columns |
| `rows` | 200 | cells | nb of rows |
| `HEIGHT` | 70 | m | vertical extent of simulation |
| `length_dcell` | 25 | m | linear size of a dcell |
| `nbiter` | 365 | timesteps (days) | total nb of timesteps (in days so X ans = X*365) |
| `NV` | 1 | cells per m | vertical nb of cells (nb per m) |
| `NH` | 1 | cells per m | horizontal nb of cells (nb per m) |
| `nbout` | 4 | — | Number of outputs |
| `nbspp` | 109 | — (number of species; read from the species file in practice) | Number of species |
| `SWtoPPFD` | 2.27 | µmol photons J⁻¹ | convert short wave irradiance to PAR photons (cf. code) |
| `p_nonvert` | 0.05 | — | light incidence param (diff through turbid medium) |
| `klight` | 0.5 | — (k_geom, GMD Part 1 Eq. 1) | light attenuation in the canopy Beer-Lambert |
| `phi` | 0.10625 | µmol C µmol⁻¹ photons | quantum yield (in micromol C/micromol photon) |
| `absorptance_leaves` | 0.83 | — (fraction) | absorptance of individual leaves |
| `theta` | 0.7 | — (curvature, Eq. 25) | parameter of the Farquhar model |
| `g1` | 3.77 | kPa^0.5 | parameter g1 of Medlyn et al ’s stomatal conductance model |
| `g0` | 5.0 | **mmol H₂O m⁻² s⁻¹** (converted in the code; see note) | minimum leaf conductance to water vapour |
| `pheno_a0` | 0.2 | — (a_T,o, GMD Part 2 Eq. 1) | threshold for change in old leaf shedding rate, in proportion of TLP |
| `pheno_b0` | 0.02 | MPa (b_T,o in GMD Part 2 Eq. 1; the file describes it "in proportion of tree height", see open questions) | threshold for change in old leaf shedding rate, in proportion of tree height |
| `pheno_delta` | 0.2 | — (δ_o) | amplitude of change in old leaf shedding rate |
| `vC` | 0.021 | — | variance of the flexion moment |
| `DBH0` | 0.005 | m | initial dbh (m) |
| `H0` | 0.95 | m | initial height (m) |
| `CR_min` | 0.3 | m | minimum crown radius (in m) |
| `CR_a` | 1.85 | — (a_cr, Eq. 17, log scale) | CR log intercept or Michaelis Menten initial growth |
| `CR_b` | 0.4445 | — (b_cr, Eq. 17) | CR log slope or Michaelis Menten asymptotic CR |
| `CD_a` | 0.0 | m (a_cd, Eq. 18) | CD intercept (absolute value) |
| `CD_b` | 0.2 | — (b_cd, Eq. 18) | CD slope (as fraction of tree height) |
| `CD0` | 0.3 | m | initial crown depth(in m) |
| `shape_crown` | 0.72 | — (0 conical to 1 cylinder) | crown shape parameter |
| `dens` | 1.0 | m² m⁻² | initial leaf density (m^2/m^2) |
| `fallocwood` | 0.35 | — (fraction) | fraction of biomass allocated to above ground wood (branch turnover+stem) |
| `falloccanopy` | 0.25 | — (fraction) | fraction of biomass allocated to canopy (leaves + reproductive organs + twigs) |
| `Cseedrain` | 50000 | seeds ha⁻¹ per timestep (whole plot, all species) | constant used to scale total seed rain per hectare across species (in next computation) |
| `nbs0` | 10 | seeds per mature tree (per dispersal event) | nb of seeds produced and dispersed by each mature tree when SEEDTRADEOFF is not defined |
| `sigma_height` | 0.19 | — (SD, log scale) | intraspecific variation in tree height (lognormal) |
| `sigma_CR` | 0.29 | — (SD, log scale) | intraspecific variation in crown radius (lognormal) |
| `sigma_CD` | 0.0 | — (SD, log scale) | intraspecific variation in crown depth (lognormal) |
| `sigma_P` | 0.24 | — (SD, log scale) | intraspecific variation in leaf phosphorus (lognormal) |
| `sigma_N` | 0.12 | — (SD, log scale) | intraspecific variation in leaf nitrogen (lognormal) |
| `sigma_LMA` | 0.24 | — (SD, log scale) | intraspecific variation in LMA (lognormal) |
| `sigma_wsg` | 0.06 | g cm⁻³ (SD, normal) | intraspecific variation in wood specific gravity |
| `sigma_dbhmax` | 0.05 | — (SD) | intraspecific variation in maximum diameter |
| `sigma_leafarea` | 0.48 | — (SD, log scale) | intraspecific variation in leaf area (lognormal) |
| `sigma_tlp` | 0.10 | — (SD, log scale) | intraspecific variation in tlp (lognormal) |
| `corr_CR_height` | 0.0 | — | correlation coefficient between crown radius and tree height |
| `corr_N_P` | 0.65 | — | correlation coefficient between leaf nitrogen and leaf phosphorus |
| `corr_N_LMA` | -0.43 | — | correlation coefficient between leaf nitrogen and LMA |
| `corr_P_LMA` | -0.39 | — | correlation coefficient between leaf phosphorus and LMA |
| `leafdem_resolution` | 1 | — | resolution of leaf demography model |
| `p_tfsecondary` | 1.0 | — (probability) | probability of secondary treefall |
| `hurt_decay` | 0 | — | parameter determining how tree damages are repaired |
| `crown_gap_fraction` | 0.15 | — (fraction; f_gap) | fraction of gaps in the crown |
| `m` | 0.045 | not established (background mortality rate) | minimal death rate |
| `m1` | 0.045 | not established | m1 (slope of death rate) |
| `Cair` | 375 | µmol mol⁻¹ | atmospheric CO2 concentration in micromol/mol |
| `PRESS` | 101.0 | kPa | atmospheric pressure in kPa |
| `_SOIL_LAYER_WEIGHT` | 2 | option (see description) | soil layer weights: relative biomass, conductance, max transpiration (0,1, 2) |
| `_WATER_RETENTION_CURVE` | 1 | option (see description) | water retention curve (Brooks&Corey: 0; Van Genuchten Mualem: 1) |
| `_LL_parameterization` | 0 | option (see description) | LL parameterizations: empirical, Kikuzawa model, and Kikuzawa model with leaf plasticity (0,1,2) |
| `_LA_regulation` | 2 | option (see description) | dynamic LA regulation: off, 1.0, 0.75, or 0.5 (0,1,2,3) |
| `_sapwood` | 1 | option (see description) | sapwood parameterizations: constant thickness (0.04), Fyllas percentage, Fyllas lower limit (0,1,2) |
| `_seedsadditional` | 0 | option (see description) | excess biomass into seeds after maturation (0,1) |
| `_NONRANDOM` | 1 | option (see description) | If _NONRANDOM == 1, the seeds for the random number generators will be kept fixed at 1, default for bug fixing |
| `_GPPcrown` | 0 | option (see description) | This defines an option to compute only GPP from the topmost value of PPFD and GPP, instead of looping within the crown. |
| `_BASICTREEFALL` | 1 | option (see description) | if defined: treefall is a source of tree death |
| `_SEEDTRADEOFF` | 0 | option (see description) | if defined: the number of seeds produced is determined by NPP allocated to reproduction and seed mass, otherwise the number of seeds is fixed |
| `_CROWN_MM` | 0 | option (see description) | Michaelis Menten allometry for crowns instead of power law, parameters have to be changed in other input sheets accordingly |
| `_OUTPUT_extended` | 1 | option (see description) | extended set of output files |
| `_OUTPUT_inventory` | 1 | option (see description) | inventory set of ouput files |

**Note on `g0`:** the code reads `g0` in mmol H₂O m⁻² s⁻¹ and converts it to mol CO₂ m⁻² s⁻¹
(`g0 *= 0.001 * GSVGSC` in `input.cpp`). Up to 4.1.3 the description in the example file
said "mol m-2 s-1"; it was corrected in 4.1.4.
