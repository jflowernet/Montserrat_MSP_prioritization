# Marine spatial planning on the Caribbean island of Montserrat: Lessons for data-limited small islands
This repository contains code and data for:

Flower, J., Ramdeen, R., Estep, A., Thomas, L.R., Francis, S., Goldberg, G., Johnson, A.E., McClintock, W., Mendes, S.R., Mengerink, K., O’Garro, M., Rogers, L., Zischka, U. & Lester, S.E. (2020). Marine spatial planning on the Caribbean island of Montserrat: Lessons for data-limited small islands. Conservation Science and Practice, e158. https://doi.org/10.1111/csp2.158

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0002-6731-8182" href="https://orcid.org/0000-0002-6731-8182" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0002-6731-8182</a></div>

## Content:

* `data`: Contains all data used in the prioritization and spatial boundaries shown in the figures.
* `outputs`: Contains all figures and data outputs from scripts
* `scripts`: Contains all data processing and analysis scripts, specifically:
  + `marine_spatial_prioritization.R`: this is the main code file. Running this code file will repeat the spatial prioritization discussed in the manuscript (see Figure 4. (a) for the result), as well as reproduce supplementary Figures S2, S3 and S5.
  + `interp_sprichness.R`: a script called from the main code file (above) that interpolates the species richness data.

## Notes:

The systematic conservation prioritization package `prioritizr`used for the spatial prioritization uses the Gurobi library which requires a license. More information about installing Gurobi: https://prioritizr.net/articles/gurobi_installation.html

Other optimizers can be used to solve the `prioritizr` conservation problems: https://prioritizr.net/articles/saltspring.html#solving-the-problem
