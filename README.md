[![R build status](https://github.com/hadexversum/powerHaDeX/workflows/R-CMD-check/badge.svg)](https://github.com/hadexversum/powerHDX/actions)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# Power calculations for Hydrogen/Deuterium Exchange experiments 


This R package implements methods for for

  - simulation of isotope distribution of deuterated peptides (theoretical spectra) based on the [code provided by Zhong-Yuan Kan](https://github.com/kanzy/HX-MS-Simulations),
  - creating deuteration curves based on theoretical spectra with multiple technical replicates and a given error structure,
  - testing differential deuterium uptake,
  - power calculations for Hydrogen/Deuterium Exchange experiments.
  
  
Development version can be downloaded from Github:

```
if (!require(devtools)) {
  install.packages('devtools')
}
devtools::install_github("hadexversum/powerHaDeX")
```

