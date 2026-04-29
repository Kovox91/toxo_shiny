# Poly(A) Shiny Supplement

Hosted app:
https://sascha-maschmann.shinyapps.io/shiny_app/

Run locally:

```r
setwd("shiny_app")
renv::restore()
options(browser = "xdg-open")
shiny::runApp(launch.browser = FALSE)
```

Required runtime files are included in this folder:
- app code in `app.R` and `R/`
- input data in `input/`
- reference tables in `reference/`
- reproducible package lockfile in `renv.lock`
