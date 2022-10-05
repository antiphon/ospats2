---
bibliography: references.bib
---

# Optimal spatio-temporal stratification

## What

Reboot of the optimal spatial stratification and sample allocation package `Ospats`. The original package is not in CRAN, and the last update to the code, as far as I can tell found here <https://github.com/brendo1001/ospats/>, was some years ago. The package date is 2016-07-19.

## Why

The method described in @degruijter2015 has potential uses where we work (Luke, Finland), but the original package is no longer being developed. Specifically, we would like to have more freedom in inputting the spatial correlations, and incorporate temporal aspects to the stratification as well.

## Plan

-   [ ] The original `Ospats::ospatsF` will be refactored and included for legacy purposes

-   [ ] New user-facing function is constructed, superseding in input flexibility the original `ospatsF`
