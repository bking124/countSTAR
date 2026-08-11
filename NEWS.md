# countSTAR 1.2.1
* Added defensive checks and conditional evaluation across the package vignette and examples to gracefully handle failures in `Suggests` dependencies.
* Fixed package check errors on specific CRAN test environments.

# countSTAR 1.2.0
* Metadata-only release; no changes to package functionality

# countSTAR 1.1.0
* Reduced package dependencies
* New transformation `bnp` for `blm_star()` and `spline_star()`
* More efficient samplers for `blm_star()` and `spline_star()`
* New `HPDregion()` to compute highest posterior (predictive) density regions
* New `rdir()` for Dirichlet sampling
* Updated `g_inv()` for stable and efficient computing of inverse `bnp` transformation
* For `randomForest_star()` and `gbm_star()`, lowered the tolerance `tol` for convergence

# countSTAR 1.0.2
* Added `roaches` dataset to package
* Fixed bug in `blm_star()` function

# countSTAR 1.0.1
* Initial version
* Added a `NEWS.md` file to track changes to the package.
