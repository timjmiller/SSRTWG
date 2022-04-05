# SSRTWG: Code Directory

All code from simulation studies

It looks like we will be primarily using WHAM for simulations and model fitting. To install you will need the devtools package:
```
install.packages("devtools")
```
The latest features are on devel branch. We may need some of these, but eventually they will be on master as a release version we can specify for reproducibility.
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
```
