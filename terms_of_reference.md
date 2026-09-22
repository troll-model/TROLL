# TROLL individual based model - Terms of Reference

_Document version 1.0_

_Date: 7 July 2026_

## 1. Description

TROLL is a computer simulator written in C++ language. TROLL represents an idealized forest stand with a typical size of 1 to 100 ha within an explicit 3D grid of 1 m³ voxels aboveground, with at most one tree established per 1×1 m ground pixel. Climatic drivers include air temperature, vapour pressure deficit, wind speed, and light intensity above the canopy, as well as precipitation. Depending on the process, TROLL operates at half-hourly temporal resolution or at daily resolution. Forest trees are modelled individually and they are generated as instances of a species, with a set of plant functional traits controlling plant physiology, and especially the carbon and water uptake, as well as demographic processes. This trait-based parameterization, together with the spatial architecture are what lets TROLL jointly track detailed forest structure, carbon fluxes, and species-level biodiversity dynamics within the same simulation.

This code was initiated in 1999 as a research-oriented tropical forest simulator ([Chave 1999](https://doi.org/10.1016/S0304-3800(99)00171-4)), with major subsequent upgrades ([Maréchaux & Chave 2017](https://doi.org/10.1002/ecm.1271), [Maréchaux et al. 2025](https://doi.org/10.5194/gmd-18-5143-2025)). The current computer code is about 10,000 lines stored on https://github.com/troll-model/TROLL. The current version is 4.0, which is delivered in the form of a single file (mainTROLL4.0.cpp). 

Since 2023, a wrapper written in the R language, rcontroll, has been developed to manage the data inputs, run a version of TROLL called by R through the Rcpp language, and manage the data outputs ([Schmitt et al. 2023](https://doi.org/10.1111/2041-210X.14215)). The rcontroll code is available https://github.com/sylvainschmitt/rcontroll. 

The goal of this document is to explain the context and organisation of TROLL. 

## 2. Governance 

The TROLL project is led and managed by a group of five core developers, henceforth the core development team. This team is a self-assembled group of scientists, committed to advancing the scientific research in forest ecology, through the development, applications and delivery of software products. 

The core development team is responsible for reviewing and accepting major merges of the TROLL code. 

They may be contacted via the generic email address: dev@troll-model.org.

**Isabelle Maréchaux**, INRAE — AMAP, Univ Montpellier, CIRAD, CNRS, INRAE, IRD, 34000 Montpellier, France
Ecologist, forest modelling, plant physiology.
Role: project lead, lead of version 4.0 upgrade (including water module)

**Fabian Jörg Fischer** — School of Biological Sciences, University of Bristol, Bristol, BS8 1TQ, United Kingdom
Ecologist, forest modelling, statistical modeling, remote sensing.
Role: lead of spatial data integration (ALS), lead of version 3 upgrade (trait variability)

**Sylvain Schmitt**, CIRAD — CIRAD, UPR Forêts et Sociétés, 34398 Montpellier, France
Ecologist, forest modelling, statistical modeling.
Role: lead of rcontroll wrapper. 

**Philippe Verley**, IRD — AMAP, Univ Montpellier, CIRAD, CNRS, INRAE, IRD, 34000 Montpellier, France
Computer scientist 
Role: lead for code factorization and code acceleration; website management

Jérôme Chave, CNRS —  CRBE, Université de Toulouse, CNRS, IRD, Toulouse INP, 118 route de Narbonne, 31062 Toulouse, France
Ecologist, forest modelling, statistical modeling.
Role: lead of versions 1.0 and 2.0

The core development team agrees to meet twice a year to explore code updates, improvements, and discuss collaborative projects.

## 3. Licence and intellectual property

TROLL version 4.0 and further developments are publicly available on GitHub as a
standalone code. TROLL is licensed under the GNU General Public License v3.0 ([GPL3.0](https://github.com/sylvainschmitt/TROLL/blob/main/LICENSE)). “Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights”.

By inheritance, rcontroll is licensed under the GNU General Public License v3.0 (GPL3.0). 

Intellectual property resides with the code developers. 

## 4. Acknowledgements

Invitations to the core development team to collaborate on research projects are welcome.

Users of the TROLL code are encouraged to contact the core development team about any questions and/or feedback (for this, users are encouraged to make use of the “issues” option on github). 

In any case, users are encouraged to acknowledge the contribution 

## 5. References

Chave, J. (1999). Study of structural, successional and spatial patterns in tropical rain forests using TROLL, a spatially explicit forest model. *Ecological modelling*, 124(2-3), 233-254. https://doi.org/10.1016/S0304-3800(99)00171-4

Maréchaux, I., & Chave, J. (2017). An individual‐based forest model to jointly simulate carbon and tree diversity in Amazonia: description and applications. *Ecological Monographs*, 87(4), 632-664.  https://doi.org/10.1002/ecm.1271

Maréchaux, I., Fischer, F. J., Schmitt, S., & Chave, J. (2025). TROLL 4.0: representing water and carbon fluxes, leaf phenology, and intraspecific trait variation in a mixed-species individual-based forest dynamics model–Part 1: Model description. *Geoscientific Model Development*, 18(16), 5143-5204. https://doi.org/10.5194/gmd-18-5143-2025

Schmitt, S., Salzet, G., Fischer, F. J., Maréchaux, I., & Chave, J. (2023). rcontroll: An R interface for the individual‐based forest dynamics simulator TROLL. *Methods in Ecology and Evolution*, 14(11), 2749-2757. https://doi.org/10.1111/2041-210X.14215


