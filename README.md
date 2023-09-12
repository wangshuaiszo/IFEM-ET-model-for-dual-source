# IFEM remote sensing ET model for dual source

## 💡 Introduction

An Independent Framework-based Evapotranspiration Model (IFEM) for dual-source: from field to regional scale.

The IFEM model is a remote sensing-based land surface evapotranspiration retrieval model. It is a user-friendly model that allows for the easy estimation of surface evaporation and transpiration components with minimal required surface observation data.

![Guide_figures/Flow chart.jpg](Guide_figures/Flow%20chart.jpg)

The IFEM model is based on the Fractional vegetation cover (FVC) - land surface temperature (LST) feature space modeling approach. In contrast to the theoretical trapezoid framework proposed by Long et al. (2012), IFEM fixes five variables and parameters that influence the linear response relationship between evaporation efficiency and land surface temperature along the LST axis. These variables and parameters include surface albedo ($\alpha$), air temperature ($T_a$), saturation vapor pressure deficit (VPD), aerodynamic resistance ($r_a$), and the ratio of soil heat flux to net radiation ($\Gamma$). In other words, these five crucial variables and parameters do not change synchronously with land surface temperature (soil moisture availability), allowing the IFEM to construct a feature space domain with uniformly distributed evaporation efficiency contours for each pixel.

![Guide_figures/Flow chart.jpg](Guide_figures/Space.jpg)
