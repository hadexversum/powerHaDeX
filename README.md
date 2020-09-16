# Power calculations for Hydrogen/Deuterium Exchange experiments 


This R package implements methods for for

  - simulation of isotope distribution of deuterated peptides (theoretical spectra),
  - creating deuteration curves based on theoretical spectra with multiple technical replicates and a given error structure,
  - testing differential deuterium uptake,
  - power calculations for Hydrogen/Deuterium Exchange experiments.
  
  
Development version can be downloaded from Github:

```
if (!require(devtools)) {
  install.packages('devtools')
}
devtools::install_github("hadexversum/powerHDX")
```

