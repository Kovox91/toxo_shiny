# Dissertation Supplement: Poly(A) Shiny App

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
