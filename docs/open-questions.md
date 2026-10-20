# Open questions for the TROLL authors

Discrepancies found while documenting TROLL 4.1.x by cross-checking the code, the GMD papers
(Part 1: Maréchaux et al. 2025; Part 2: Schmitt et al. 2025) and the example files. None of
them was changed in the documentation releases: each needs a decision by the authors. Items
marked **results** would change model outputs if resolved in favour of the paper.

| # | Topic | Where | Discrepancy | Impact |
| --- | --- | --- | --- | --- |
| 1 | Leaf temperature update | `Tree::FluxesLeaf` | `TLEAF1 = Ta + 0.25 * TDIFF`: the comment says the factor slows convergence, but as written it divides the converged leaf-to-air temperature difference by 4 | **results** |
| 2 | Aerodynamic resistance of the soil surface | `simulation.cpp`, `r_aero` | code coefficient 43.17 = ln(1000)/κ²; Part 1 Eq. 15 gives (ln 1000)²/κ² = 298 (or 86 for ln(1000²)): soil evaporation ~7x higher in the code if the paper is right; Part 2 reports ET overestimated at one site | **results** |
| 3 | Canopy micro-climate parameters | `LookUp_T`, `LookUp_VPD` | ΔT = 3 °C, LAI_sat = 7 and C_VPD0 = 0.25 are hard-coded, while Part 1 presents them as set parameters; should they be input parameters? | design |
| 4 | Interception in `water_balance.txt` | `output.cpp` | summed in mm over soil cells, divided by the number of sites: too small by (soil cells)/(sites) | outputs only |
| 5 | Throughfall unit | `simulation.cpp` vs `output.cpp` | computed in mm per soil cell, commented as m³ when written; inside the model, interception/throughfall are in mm and evaporation/leakage in m³ | outputs; check fluxes |
| 6 | `g0` unit | `example/global_inputs.txt` | the code expects mmol H₂O m⁻² s⁻¹; the example file said mol m⁻² s⁻¹ (description corrected in 4.1.4); confirm and update user documentation (rcontroll) | users |
| 7 | Leaf N and P units | Part 1 Sect. 2.4.1 vs Part 2 and files | Part 1 gives mg g⁻¹, Part 2 and the species file g g⁻¹ (the code converts to mg g⁻¹); erratum? | documentation |
| 8 | DBH threshold and LMA units in Part 2 | Part 2 Sect. 2.1, Fig. 5 | dbh_thres in cm in the text, metres in Fig. 5 and the files; Fig. 5 gives LMA in g cm⁻³ instead of g m⁻² | documentation |
| 9 | `pheno_b0` | `global_inputs.txt` vs Part 2 Eq. 1 | described "in proportion of tree height" in the file; enters Eq. 1 as an absolute potential (MPa) | documentation |
| 10 | Units not established from the code | outputs `state_*`, `site*` | tree `transpiration`, `root_biomass_<l>` | outputs |
| 11 | Unused or broken code | `PROB_S`, `deltaR`; options `TROLL_NO_WATER`, `TROLL_NO_CROWN_UMBRELLA`, `TROLL_NO_FULL_CLIMATE`, `Output_ABC`, `TRACK_INDIVIDUALS` (uses `t_timeofyear_born`, declared `t_month_born`) | remove or repair? | maintenance |
| 12 | Open notes in the code | 47 comments such as "should be", "to discuss" | to be triaged into issues | maintenance |
