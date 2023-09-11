# _airportscreening_: Effectiveness of airport screening at detecting infected travellers

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

_airportscreening_ contains the code for a Shiny app whose purpose is to display the effectiveness of screening arrivals for infectious diseases.
We assume that the infection may have a latent stage during which travellers are infected but not symptomatic.
They may become symptomatic during their flight, which may be detected by entry screening, or after their arrival, and therefore no screening would have detected an infection.

This package has been created by converting the online Shiny app, which was previously found at https://cmmid-lshtm.shinyapps.io/traveller_screening/, and which was initially created as part of the CMMID response to the Covid-19 pandemic in early 2020.

## Quick start

The package can be installed from Github using:

```r
# install the devtools or remotes packages first
devtools::install_github("bquilty25/airport_screening")
```

Run the app using the `run_app` function.

```r
airportscreening::run_app()
```

| Authors |
| :-- |
| Mr. Billy Quilty |
| Assoc. Prof. Stefan Flasche |
| Dr. Sam Clifford |
| Assoc. Prof. Rosalind Eggo |
| Other members of CMMID at LSHTM |

![Screenshot of app](inst/info/figures/app_screenshot.png)

Users may adjust the following parameters:

* duration of flight
* mean and variance of the probability distribution describing time between infection and onset of symptoms
* mean and variance of the probability distribution describing time between onset of symptoms and severe symptoms
* sensitivity of exit screening
* sensitivity of entry screening

