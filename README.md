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
