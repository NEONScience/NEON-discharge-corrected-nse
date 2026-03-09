# Calculate Nash-Sutcliffe Efficiencies (NSE) and Regression Coefficients for Corrected NEON Discharge Data

Respository for the R script used to analyze data and generate the table and figure published in the manuscript 'Evaluation of streamflow estimates generated at aquatic sites across the National Ecological Observatory Network, Water Years 2022-2024' by Harrison et al. (submitted for publication in 2026 to _Limnology and Oceanography Letters_).

## Data Processing Script

[/src/neon_discharge_nse.R](https://github.com/NEONScience/NEON-discharge-corrected-nse/blob/main/src/neon_discharge_nse.R)

**Description:** First, field and modeled discharge (15-min resolution) are downloaded from the NEON Data Portal. Each field discharge value is mapped to the nearest modeled discharge timestamp. Plots are generated for each site showing the linear fit of field vs. modeled discharge, then stacked into a single plot for output. Finally, an output table is generated that contains the following summary statistics for each site: Nash-Sutcliffe Efficiency (NSE), slope of linear equation, coefficient of determination (R^2) of linear equation, the percent of the period of record that contains NULL modeled discharge, the channel slope of the site, and the watershed area of the site.

## Citations: Nash-Sutcliffe Efficiency (NSE)

Nash, J.E. and J.V. Sutcliffe, 1970. “River flow forecasting through conceptual models part I – A discussion of principles”. Journal of Hydrology, no. 10: 282-290. https://doi.org/10.1016/0022-1694(70)90255-6.

Zambrano-Bigiarini, M. (2024) hydroGOF: Goodness-of-fit functions for comparison of simulated and observed hydrological time series R package version 0.6-0.1. https://cran.r-project.org/package=hydroGOF. doi:10.5281/zenodo.839854.

## Citations: NEON Data Products

NEON (National Ecological Observatory Network). Continuous discharge (DP4.00130.001), RELEASE-2026. https://doi.org/10.48443/4n6c-gc44. Dataset accessed from https://data.neonscience.org/data-products/DP4.00130.001/RELEASE-2026 on January 26, 2026.

NEON (National Ecological Observatory Network). Discharge field collection (DP1.20048.001), RELEASE-2026. https://doi.org/10.48443/qxb5-9d50. Dataset accessed from https://data.neonscience.org/data-products/DP1.20048.001/RELEASE-2026 on January 26, 2026.

NEON (National Ecological Observatory Network). Stream morphology map (DP4.00131.001), RELEASE-2026. https://doi.org/10.48443/57z7-s858. Dataset accessed from https://data.neonscience.org/data-products/DP4.00131.001/RELEASE-2026 on January 26, 2026.

NEON (National Ecological Observatory Network). Stage-discharge rating curves (DP4.00133.001), RELEASE-2026. https://doi.org/10.48443/pzm2-za41. Dataset accessed from https://data.neonscience.org/data-products/DP4.00133.001/RELEASE-2026 on January 26, 2026.








