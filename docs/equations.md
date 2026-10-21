# Equations of the model description and their implementation

Where the equations of the TROLL 4.0 model description (Maréchaux et al. 2025, *Geosci. Model
Dev.* 18, 5143–5204, <https://doi.org/10.5194/gmd-18-5143-2025>, CC BY 4.0) are implemented in
TROLL 4.1.5, and whether the constants in the code match the paper. Covered: Eqs. 1–68
(Sects. 2.2–2.7, the whole model description). Constants were compared automatically against the code of TROLL 4.1.5.

Status: **matches**: formula and constants checked against the paper; **mapped**: the implementing
code was identified, but the formula was not checked line by line; **differs**: the code differs
from the paper (see [open-questions.md](open-questions.md)).

## Summary

| Section | Equations | Entries | Matches | Mapped | Differs |
| --- | --- | --- | --- | --- | --- |
| 2.2 Light and above-ground micro-climate | 1–9 | 6 | 6 | 0 | 0 |
| 2.3 Soil water | 10–15 | 5 | 2 | 2 | 1 |
| 2.4 Tree structure (allometry and roots) | 16–23 | 7 | 4 | 3 | 0 |
| 2.5 Leaf physiology | 24–50 | 18 | 15 | 2 | 1 |
| 2.6 Carbon allocation and growth | 51–62 | 12 | 10 | 1 | 1 |
| 2.7 Demography (recruitment and mortality) | 63–68 | 7 | 5 | 2 | 0 |
| **All** | 1–68 | **55** | **42** | **10** | **3** |

## 2.2 Light and above-ground micro-climate (Eqs. 1–9)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 1 | PPFD above a layer, Beer–Lambert: PPFD_top exp(−k LAI) | `LookUp_flux` (`lookup_tables.cpp`) | k = `klight` × `absorptance_leaves` (`kpar`, `input.cpp`) | **matches** |
| 2–3 | PPFD absorbed in a layer | `LookUp_flux`, `LookUp_flux_absorption` (`lookup_tables.cpp`) | exp(−k LAI)(1 − exp(−k dens)) | **matches** |
| 4, 6 | Air temperature within the canopy: T_top − ΔT LAI / LAI_sat | `LookUp_T` (`lookup_tables.cpp`) | 0.4285714 = ΔT / LAI_sat with **ΔT = 3 °C, LAI_sat = 7** | **matches**; ΔT and LAI_sat are hard-coded, not input parameters (open question) |
| 5, 7 | VPD within the canopy | `LookUp_VPD` (`lookup_tables.cpp`) | 0.25 + √(0.0803571 (7 − LAI)): **C_VPD0 = 0.25**, (1 − 0.25)²/7 = 0.0803571; layer mean uses 0.188982 = ⅔ × 0.75/√7 | **matches**; C_VPD0 hard-coded |
| 8 | Wind above the canopy, logarithmic profile | `LookUp_Wind` (`lookup_tables.cpp`); `TopWindSpeed_DCELL` (`memory.cpp`) | 16.67 = 1/0.06 (z₀ = 0.06 H), d = 0.8 H | **matches** |
| 9 | Wind within the canopy (Inoue) | `LookUp_Wind` (`lookup_tables.cpp`) | `alphaInoue` = 3 | **matches** |

## 2.3 Soil water (Eqs. 10–15)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 10 | Soil water balance per layer | `UpdateField` (`simulation.cpp`) | — | **mapped** |
| 11 | Interception min(P, K LAI) | `Interception` (`memory.cpp`) | K = 0.2 mm | **matches** (output scaling: open question) |
| 12–13 | Soil evaporation | `Evaporation` (`memory.cpp`) | — | **mapped** |
| 14 | Soil surface resistance | `r_soil` (`simulation.cpp`) | 8.206, 4.255 | **matches** |
| 15 | Aerodynamic resistance ln(Z/Z_m)² / (κ² u) | `r_aero` (`simulation.cpp`) | **43.17347 = ln(1000)/0.4²** | **differs**: the paper gives (ln 1000)²/0.16 = 298 (or ln(1000²)/0.16 = 86); open question |

## 2.4 Tree structure (allometry and roots) (Eqs. 16–23)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 16 | Height h_lim dbh / (a_h + dbh) | `UpdateHeight` (`tree.cpp`) | — | **mapped** |
| 17 | Crown radius exp(a_cr) dbh^b_cr | `UpdateCR` (`tree.cpp`) | `CR_a`, `CR_b` | **mapped** |
| 18 | Crown depth min(h/2, a_cd + b_cd h) | `UpdateCD` (`tree.cpp`) | `CD_a`, `CD_b` | **mapped** |
| 19 | Root depth 0.35 dbh^0.54, dbh in cm | `t_root_depth` (`tree.cpp`) | 0.35, 0.54, dbh × 100 | **matches** |
| 20 | Root profile exp(−3 z / RD) |  | 3.0 | **matches** |
| 21 | ψ_root weighting, ψ_R,min |  | −3 MPa | **matches** |
| 22–23 | Soil-to-root conductance (Gardner) |  | SRL = 10 m g⁻¹ | **matches** |

## 2.5 Leaf physiology (Eqs. 24–50)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 24 | FvCB net assimilation, A_j = J/4 (c_i − Γ*)/(c_i + 2Γ*) | `Tree::Photosyn` (`tree.cpp`) | same c_i in numerator and denominator | **matches** since 4.1.2 (was c_i of the Rubisco-limited solution) |
| 25 | Electron transport rate J | `ElectronTransportRate` (`quad.hpp`) | — | **matches** (numerically stable form since 4.1.2) |
| 26 | Vcmax per leaf mass at 25 °C (Domingues et al. 2010) | `CalcVcmaxm` (`crown.cpp`) | −1.56, 0.43, 0.37 / −0.80, 0.45, 0.25, on log₁₀ values, N and P converted to mg g⁻¹ | **matches** the original log–log form; the paper omits the log₁₀ (erratum) |
| 27 | Jmax per leaf mass at 25 °C | `Tree::CalcJmaxm` (`tree.cpp`) | −1.50, 0.41, 0.45 / −0.74, 0.44, 0.32, log₁₀ form | **matches**, as Eq. 26 |
| 28–29 | Temperature dependence of Vcmax and Jmax (Bernacchi et al. 2003) | `LookUp_VcmaxT`, `LookUp_JmaxT` (`lookup_tables.cpp`) | 26.35, 65.33; 17.57, 43.54 | **matches**; R = 0.00831 in the code, 0.008314 in the paper (0.04%) |
| 30–31 | Γ* and Km (von Caemmerer 2000) | `LookUp_GammaT`, `LookUp_KmT` (`lookup_tables.cpp`) | 37, 23.4; 404, 59.36, 210, 248, 35.94 | **matches** (same R as above) |
| 32 | Leaf dark respiration at 25 °C (Atkin et al. 2015) | `CalcRdark` (`crown.cpp`) | code: area-based 1.3893 + 0.0728 N_area + 0.0015 P_area + 0.0095 Vcmax − 0.0358 × 26.2 | **differs**: the paper gives the mass-based model (8.5341 − 0.1306 N − …); open question |
| 33 | Temperature dependence of dark respiration | `LookUp_Rleaf` (`lookup_tables.cpp`) | (3.09 − 0.0215 (T + 25))^((T − 25)/10) | **matches** (0.043/2 = 0.0215) |
| 34–35 | Stomatal conductance (Medlyn et al. 2011) | `Tree::Photosyn` (`tree.cpp`) | g0, g1, 1.6 (`GSVGSC` 1.57 for CO2/H2O) | **mapped** |
| 36 | g1 = (6.53 − 3.97 wsg) × WSF_s (Lin et al. 2015) | `t_g1_0`, `Water_availability` (`tree.cpp`) | 6.53, 3.97; WSF_s applied | **matches** |
| 37 | Pre-dawn water potential ψ_root − 0.01 h | `Water_availability` (`tree.cpp`) | 0.01 | **matches** |
| 38–39 | Stomatal water stress factor exp(−2.23 ψ_pd / π_tlp) | `Water_availability`, `t_b` (`tree.cpp`); `t_WSF` (`lookup_tables.cpp`) | −2.23 | **matches** |
| 40 | Non-stomatal water stress factor (1 + (ψ_pd/π_tlp)⁶)⁻¹ | `Water_availability` (`tree.cpp`); `t_WSF_A` (`lookup_tables.cpp`) | exponent 6 | **matches** |
| 41 | Penman–Monteith leaf transpiration | `Tree::FluxesLeaf` (`tree.cpp`) | — | **mapped**; leaf temperature update: open question (factor 0.25) |
| 42 | Latent heat of vaporization | `FluxesLeaf`, `LHV` (`tree.cpp`) | 2.501 × 10⁶, 2365, 18 × 10⁻³ | **matches** |
| 43–44 | Isothermal net radiation, net longwave (Dilley & O'Brien 1998) | `LookUp_INLR`, `LookUp_ExtinctLW` (`lookup_tables.cpp`) | k_d = 0.8 | **matches** |
| 45–48 | Heat conductances: free and forced convection, radiation | `FluxesLeaf` (`tree.cpp`); `LookUp_GRADN` (`lookup_tables.cpp`) | 56.2341 = 0.5 × (1.6 × 10⁸)^0.25, 0.003, 4 σ T³/(Cp Ma) | **matches** |
| 49–50 | Conductance to water vapour | `FluxesLeaf` (`tree.cpp`) | 1.075 | **matches** |

## 2.6 Carbon allocation and growth (Eqs. 51–62)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 51 | NPP = GPP − maintenance − growth respiration (30%) | `Tree::CalcNPP` (`tree.cpp`) | 0.7 | **matches** |
| 52 | GPP weighting of young and old leaves, ρ = 0.5 | `ComputeRespGPP`, `effLA` (`tree.cpp`) | 0.5 (LA + LA_mature) | **matches** |
| 53 | Respiration and transpiration of young and old leaves, ρ′ = 0.75 | `ComputeRespGPP`, `effLA_night`, `effLAT` (`tree.cpp`) | 0.75 LA + 0.25 LA_mature | **matches** |
| 54 | Stem respiration 39.6 µmol m⁻³ s⁻¹, Q10 = 2 | `LookUp_Rstem` (`lookup_tables.cpp`) | 39.6, Q10 = 2 | **matches** |
| 55 | Sapwood area (Fyllas et al. 2014) | `Tree::UpdateSapwoodArea` (`tree.cpp`) | 0.066, 0.017, **−0.18**, 1.6; minimum 0.5 cm | **differs**: the paper gives δ₁ = −0.018; open question |
| 56 | Leaf age pools; f_leaves = 0.68 f_canopy, NPP_wood = 0.6 f_wood NPP | `Tree::UpdateLeafDynamics` (`tree.cpp`); `Growth` (`simulation.cpp`) | 0.68, 0.6 | **matches** |
| 57 | Leaf lifespan max(3, 12.755 exp(0.007 LMA − 0.565 N)) months | `Tree::CalcLeafLifespan` (`tree.cpp`) | 12.755, 0.007, 0.565, minimum 3 | **matches** (constants); N enters in g g⁻¹, open question on the unit |
| 58 | Old-leaf shedding threshold | `UpdateLeafDynamics`, `Thres` (`tree.cpp`) | `pheno_a0`, −0.01 h − `pheno_b0` | **matches** |
| 59 | Maximum carbon storage | `Tree::CalcCarbonStorageMax` (`tree.cpp`) | 1000 × 0.5 × 0.05 × 1.25 | **matches** |
| 60 | AGB = 0.0559 wsg dbh² h (Chave et al. 2014) | `Tree::CalcAGB` (`tree.cpp`) | 0.0559, dbh in cm | **matches** |
| 61 | Stem volume increment | `Tree::Growth` (`simulation.cpp`); `CalcIncrementDBH` (`tree.cpp`) | 0.5 (C fraction) | **mapped** |
| 62 | Size-related growth decline | `Tree::Growth` (`simulation.cpp`) | 3 − 2 dbh / dbh_thres | **matches** |

## 2.7 Demography (recruitment and mortality) (Eqs. 63–68)

| Eq. | Process | Implementation | Constants checked | Status |
| --- | --- | --- | --- | --- |
| 63 | Maturity threshold 0.5 dbh_thres | `Tree::Birth`, `t_dbhmature` (`tree.cpp`) | 0.5 | **matches** |
| 64 | External seed rain N_tot f_reg n_ha | `UpdateSeeds` (`simulation.cpp`) | — | **matches** (non-integer hectares fixed in 4.1.0) |
| 65 | Death rate as a sum of components | `Tree::DeathRate` (`tree.cpp`) | — | **mapped** |
| 66 | Background mortality m (1 − wsg / wsg_lim) | `Tree::DeathRate` (`tree.cpp`) | wsg_lim = 1 g cm⁻³ | **matches** |
| 67 | Treefall height threshold | `Tree::Birth` (`tree.cpp`); `TriggerTreefall` (`simulation.cpp`) | — | **mapped** |
| 68 | Lethal water potential −0.9842 + 3.1795 π_tlp (Bartlett et al. 2016) | `Tree::Birth`, `t_phi_lethal` (`tree.cpp`) | −0.9842, 3.1795 | **matches** |
| — | Recruitment conditions: top-soil ψ > π_tlp/2; initial LA = 0.25 LA_opt | `RecruitTree` (`simulation.cpp`); `Tree::Birth` (`tree.cpp`) | 0.5, 0.25 | **matches** |

Unit questions also remain on Eq. 57 and on Eqs. 26–27 as printed in the paper (see open questions 15 and 16).
