# Optimal stratification under correlations

## What

The R-package `ospats2` is a reboot of the optimal spatial stratification and sample allocation package `Ospats`, implementing the eponymous algorithm describe by De Gruijter, Minasny, and Mcbratney (2015)[^readme-2], and written by Brendan Malone and Nicolas Saby.

The original package is not in CRAN (if ever was), and the closest thing to an official R-package as I can tell can be found here <https://github.com/brendo1001/ospats/>. Last relevant update in the repository dates back to 2016, the package date is 2016-07-19.

There is another repository from the authors containing Julia code for a 2017 by De Gruijter and others. Code can be found here <https://github.com/jjdegruijter/ospats>. The code is written for Julia versions before 1.0 so is very much deprecated (Updated code in this repo for comparisons).

## Why

The method described by De Gruijter, Minasny, and Mcbratney has potential uses where we work (Luke, Finland), but the original package is no longer being developed. We would like to use and develope the method further. For example, we would like more freedom in inputting the spatial correlation structure, and be able to work with spatio-temporal populations.


## Plan

-   [x] The original `Ospats::ospatsF` is refactored and included for comparison purposes (`ospatsF_ref`)
-   [x] The Julia version is translated to R (`ospatsF_julia`; Just the stratification parts)
-   [ ] New user-facing function is constructed, superseeding in input flexibility the original `ospatsF`
-   [ ] Cross-compatibility with other packages, both of inputs and outputs.
-   [ ] Optimize the code

Design principle is to be close to and/or compatible with the `spsurvey` package (which seems to be alive as of 2022-10-05). It implements the Generalized Random Tesselation Stratification (GRTS) algorithm which does a similar job to `Ospats`, but has also loads of functions for post-stratification tasks such as sampling and estimation.


[^readme-2]: De Gruijter, J. J., B. Minasny, and A. B. Mcbratney. 2015. "Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error." Journal of Survey Statistics and Methodology 3 (1): 19--42. <https://doi.org/10.1093/jssam/smu024>.
