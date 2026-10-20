# Equations of the model description and their implementation

Where the equations of the TROLL 4.0 model description (Maréchaux et al. 2025, *Geosci. Model
Dev.* 18, 5143–5204, <https://doi.org/10.5194/gmd-18-5143-2025>, CC BY 4.0) are implemented in
TROLL 4.1.4, and whether the constants in the code match the paper. Covered so far: Eqs. 1–25
(Sects. 2.2–2.5.1). Constants were compared automatically against the code of TROLL 4.1.4.
"Matches" means the formula and its constants were checked; "mapped" means the implementing
function was identified but the formula was not checked line by line.

| Eq. | Process | Implementation | Constants in the code | Status |
| --- | --- | --- | --- | --- |
| 1 | PPFD above a layer, Beer–Lambert: PPFD_top exp(−k LAI) | `LookUp_flux` (`lookup_tables.cpp`) | k = `klight` × `absorptance_leaves` (`kpar`, `input.cpp`) | matches |
| 2–3 | PPFD absorbed in a layer | `LookUp_flux` / `LookUp_flux_absorption` | exp(−k LAI)(1 − exp(−k dens)) | matches |
| 4, 6 | Air temperature within the canopy: T_top − ΔT LAI / LAI_sat | `LookUp_T` | 0.4285714 = ΔT / LAI_sat with **ΔT = 3 °C, LAI_sat = 7** | matches; ΔT and LAI_sat are hard-coded, not input parameters (open question) |
| 5, 7 | VPD within the canopy | `LookUp_VPD` | 0.25 + √(0.0803571 (7 − LAI)): **C_VPD0 = 0.25**, (1 − 0.25)²/7 = 0.0803571; layer mean uses 0.188982 = ⅔ × 0.75/√7 | matches; C_VPD0 hard-coded |
| 8 | Wind above the canopy, logarithmic profile | `LookUp_Wind`, `TopWindSpeed_DCELL` | 16.67 = 1/0.06 (z₀ = 0.06 H), d = 0.8 H | matches |
| 9 | Wind within the canopy (Inoue) | `LookUp_Wind` | `alphaInoue` = 3 | matches |
| 10 | Soil water balance per layer | `UpdateField` (`simulation.cpp`) | — | mapped |
| 11 | Interception min(P, K LAI) | `Interception` (`simulation.cpp`) | K = 0.2 mm | matches (output scaling: open question) |
| 12–13 | Soil evaporation | `Evaporation` (`simulation.cpp`) | — | not checked in detail |
| 14 | Soil surface resistance | `r_soil` (`simulation.cpp`) | 8.206, 4.255 | matches |
| 15 | Aerodynamic resistance ln(Z/Z_m)² / (κ² u) | `r_aero` (`simulation.cpp`) | **43.17347 = ln(1000)/0.4²** | **differs**: the paper gives (ln 1000)²/0.16 = 298 (or ln(1000²)/0.16 = 86); open question |
| 16 | Height h_lim dbh / (a_h + dbh) | `UpdateHeight` | — | mapped |
| 17 | Crown radius exp(a_cr) dbh^b_cr | `UpdateCR` | `CR_a`, `CR_b` | mapped |
| 18 | Crown depth min(h/2, a_cd + b_cd h) | `UpdateCD` | `CD_a`, `CD_b` | mapped |
| 19 | Root depth 0.35 dbh^0.54, dbh in cm | `t_root_depth` (`tree.cpp`) | 0.35, 0.54, dbh × 100 | matches |
| 20 | Root profile exp(−3 z / RD) | root biomass per layer (`tree.cpp`) | 3.0 | matches |
| 21 | ψ_root weighting, ψ_R,min | soil layer weights (`tree.cpp`) | −3 MPa | matches |
| 22–23 | Soil-to-root conductance (Gardner) | soil layer weights (`tree.cpp`) | SRL = 10 m g⁻¹ | matches |
| 24 | FvCB net assimilation, A_j = J/4 (c_i − Γ*)/(c_i + 2Γ*) | `Tree::Photosyn` | same c_i in numerator and denominator | matches since 4.1.2 (was c_i of the Rubisco-limited solution) |
| 25 | Electron transport rate J | `ElectronTransportRate` (`quad.hpp`) | — | matches (numerically stable form since 4.1.2) |

Equations 26 onwards (stomatal conductance, energy balance, respiration, allocation, phenology,
recruitment, mortality, seed dispersal) remain to be mapped.
