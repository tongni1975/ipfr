## Test environments
* Windows 10, R 3.6.1 (local)
* Ubuntu 14.04.05 LTS, R 3.6.1 (Travis-CI)
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-Hub)
* Fedora Linux, R-devel, clang, gfortran (R-Hub)

## R CMD check results
There were no ERRORs or WARNINGs

check_rhub() created a note about possible word mispellings
(Reweighting, Xin, al, et, rebalancing). I checked each, and they are spelled
correctly.

## Downstream dependencies
There are currently no downstream dependencies for this package

## Changes from last submission
I removed the options(scipen=999) line from both vignettes. I also made a few
revisions to the vignette prose for clarity.