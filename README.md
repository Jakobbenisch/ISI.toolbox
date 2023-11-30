# ISI.toolbox

## Development Tips
```R
## Installation

library(devtools)

devtools::install_git(
  "https://github.com/Jakobbenisch/ISI.toolbox.git"
)

## Packaging tools

library(devtools)

package.path = "/home/heinrich/gitrepos/DatapoolR/"     # path must point to the folder containing the package

## simulate a new package installation
load_all(package.path)

## run tests
test(package.path)            # this runs the tests in the `test` folder of the package
test_coverage(package.path)       # needs package 'covr' to be installed

## build documentation (uses Roxygen2)
document(package.path)

## run R CMD checks
check(package.path, cran=TRUE, manual=TRUE)

## build package for CRAN submission
build(pkg=package.path)
```
