#' @details
#' Reboot of `Ospats`.
#'
#' @section What:
#'The R-package `ospats2` is a reboot of the optimal spatial stratification and sample allocation package `Ospats`, implementing the eponymous algorithm describe by De Gruijter, Minasny, and Mcbratney (2015), and written by Brendan Malone and Nicolas Saby.
#'
#'The original package is not in CRAN (if ever was), and the closest thing to an official R-package as I can tell can be found here <https://github.com/brendo1001/ospats/>. Last relevant update in the repository dates back to 2016, the package date is 2016-07-19.
#'
#'There is another repository from the authors containing Julia code for a 2017 by De Gruijter and others. Code can be found here <https://github.com/jjdegruijter/ospats>. The code is written for Julia versions before 1.0 so is very much deprecated (Updated code in this repo for comparisons).
#' @section Why:
#'
#' The method described by De Gruijter, Minasny, and Mcbratney has potential uses where we work (Luke, Finland), but the original package is no longer being developed. We would like to use and develope the method further. For example, we would like more freedom in inputting the spatial correlation structure, and be able to work with spatio-temporal populations.
#'
#' See the vignettes for more details.
#'
#'
"_PACKAGE"
