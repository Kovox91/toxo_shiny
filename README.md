# Dissertation Supplement: Poly(A) Shiny App
[![DOI](https://zenodo.org/badge/1198860727.svg)](https://doi.org/10.5281/zenodo.19898061)
Hosted app:
https://sascha-maschmann.shinyapps.io/shiny_app/

Local run:

```r
setwd("shiny_app")
renv::restore()
options(browser = "xdg-open")
shiny::runApp(launch.browser = FALSE)
```

Licensing:
- Code: MIT, see LICENSE.
- Data: CC BY 4.0, see LICENSE-data.

Citation metadata:
- CITATION.cff is included for GitHub and Zenodo DOI workflows.
