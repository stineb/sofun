# Input processing for SOFUN

Beni Stocker, 20 Sep 2016

This repository holds scripts to process input data into a standard format, readable by SOFUN. Data is retrieved from original data files available on CX1 and is stored into a local standardised directory tree from where it's accessed by SOFUN. This is illustrated on https://docs.google.com/presentation/d/1XcZldCMMCRMk2pwlHinod5-8Q1PvCMH2ccZZO9h4U3A/edit?usp=sharing

## Get site-scale data
### Get climate data
To retrieve daily (and quasi-daily based on monthly) input data for variables temperature, precipitation, and vapour pressure, run the following R script:
```R
source( "get_climate.R" )
```
When working on Imperial's server CX1, make sure to load the following modules beforehand:
```sh
module load R
module load nco
module load netcdf
```
Use `qsub submit_get_climate.sh` to execute this on the cluster.

### Get fAPAR data
#### MODIS EVI cutouts at FLUXNET sites
Run in R:
```R
install.packages("devtools")
library(devtools)
install_github("seantuck12/MODISTools", build_vignettes=TRUE)
library( MODISTools )
source( "get_evi_modissubset.R" )
```

#### FAPAR3G data
Run in R:
```R
source( "get_fapar_fapar3g.R" )
```


### Get N deposition
Run in R:
```R
source( "get_ndep.R" )
```
