
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WolfPackR

<!-- badges: start -->
<!-- badges: end -->

The goal of WolfPackR is to identify and analyze wolf pack structures by
integrating genetic relatedness and spatial data. Using a relatedness
estimator, it constructs genetic networks to define potential packs,
then refines these groups by calculating Minimum Convex Polygons (MCPs)
to account for spatial territory overlap. The package identifies “lone
individuals” (genetically or spatially isolated) and “ugly ducklings”
(genetically linked but spatially disconnected), allowing for a nuanced
understanding of pack dynamics. WolfPackR leverages igraph for network
analysis, sf for spatial operations, and leaflet for
visualization, making it a comprehensive tool for ecologists studying
social structures in wolf populations.

## Installation

You can install the development version of WolfPackR from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("eboncourt/WolfPackR")
```
## Vignette
You can find a compiled vignette following [this link](https://github.com/eboncourt/WolfPackR/blob/master/inst/doc/WolfPackR_vignette.html) or from
``` r
vignette("WolfPackR_vignette")
```
