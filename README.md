# SOFUN 
Modelling framework for simulating terrestrial ecoystems and their biogeochemistry
SOFUN stands for Seasonal optimisation of fixation and uptake of nitrogen.

This is a modular framework that serves as a structure for simulating terrestrial ecosystem functioning (radiation, photosynthesis, allocation, soil organic dynamics, inorganic nitrogen dynamics). Code is programmed in Fortran 90 and can be compiled using PGI and Intel compilers.
Written, developed and maintained by Beni Stocker (b.stocker@imperial.ac.uk).

## Table of Contents
-------------------
[TOC]

## Repository Structure
----------------------
`src/`
This directory holds all source code.

`input/`
This directory holds all input files to drive the simulations.

  `input/sitedata/`
  Holds subdirectories for climate, co2, fapar, and ndep input data; structured for site-scale simulations.

`output/`
Holds model output files.


## Usage
---------
1. Clone git repository into a local directory called 'sofun' by `mkdir sofun; git clone https://<username>@bitbucket.org/labprentice/sofun.git sofun/`
2. Create input data (see repository 'getin')
3. Create simulation parameter files (see repository 'utils_sofun/prepare_input')
4. Compile
  This can be done at different levels of model integration. For example, to run the P-model setup only, compile only relevant modules for simulating GPP and the soil water balance by
```bash
  make pmodel
```
And execute the P-model by
```bash
  echo RUNNAME | ./runpmodel
```
To compile the full model, do 
```bash
  make cnmodel
```
Adjust Makefile and src/Makefile to use a different compiler.
Change 'RUNNAME' to any given simulation name. Parameter files containing 'RUNNAME' in their file name need to be adjusted accordingly (replacing the string 'RUNNAME' with your simulation name and adjusting parameter values.).


## Model structure
-------------------------

sofun is designed to be a modular framework that can adopt different formulation of processes by selecting modules. Irrespective of the module choice, a set of state variables must be calculated and updated by the respective subroutines which are contained in the chosen modules. These required state variables are declared in modules 'fluxes', 'pools', and 'treegeometry'. Other module-specific variables are declared within the module. Module-specific output variables are also declared only in the respective module. These "feature-"modules may also contain other subroutines to read/write necessary inputp/output that is not required by the model outside the module.

The separation between the main program 'sofun' and the subroutine 'biosphere' is somewhat arbitrary and I have not been able to follow a clear distinction between what is contained in modules and what is passed on as arguments.

Non-module-specific subroutines are defined in .F files while all module-specific subroutines are specified in their respective modlue (.mod.F). 

### Coding philosophy
- Functions are strictly self-contained: public variables are updated within functions
- Subroutines are used to update public variables

### State variables

In order to satisfy modularity, inter-changeable subroutines have to use the same set of global variables, but may treat them differently and use specific, locally defined variables. Still, the the core remains the same. That is, the vegetation is described by a set of state variables. These are:

*in module 'pools'
pleaf : leaf mass;
proot : root mass; 
plabl : labile pool; 
pexud : exudates in soil;
plitt_af : above-ground litter (fast turnover);
plitt_as : above-ground litter (slow turnover);
plitt_bg : below-ground litter;
psoil_sl : soil organic matter (slow turnover);
psoil_fs : soil organic matter (fast turnover);
wtot : soil water storage;
nh4 : ammonium;
no3 : nitrate;
no2 : nitrogen dioxide;
 
### Dimensionality of mass pools/fluxes

The major loops determine the dimensionality of variables. The highest-order subroutine is 'biosphere' which simulates C exchange over the selected domain. SR 'biosphere' is called each year. Within 'biosphere', the program runs over loops for each gridcell, month, and day (in this order). This may be changed to 'biosphere' being called each day, and only the gridcell loop being governed inside 'biosphere'. 

Model state variables generally have the lowest-possible number of dimensions. In most cases, this is 'pft' (for each plant functional type) or 'lu' (for each gridcell tile). An exception are pool variables which have a space dimension (representing gridcells 'jpngr'). Only input/output variables (to be kept separate from model state variables!) have an explicit time dimension (for month or day). 

Variables are generally named so that the first letter represents the time scale on which the variable is updated (d for daily, m for monthly, a for annually). Output variables are named, e.g., like 'outdgpp' (for daily GPP output), or 'luoutagpp' (for annual GPP, summed over all PFTs per LU). 

Fortran-derived types are used to define pools of organic matter (roots, sapwood, leaves, roots, litter, soil). This type consists of "dimensions" 'carbon' and 'nitrogen', which in turn consist of 'c12' and 'n14', respectively. There are also pools/fluxes that consist only of 'carbon' or 'nitrogen'. As a principle, the additional dimensionality of such derived-types should be introduced on the highest-possible level. E.g. it makes no sense to define GPP as 'carbon', as all variables NPP, Ra, ... carry the same signature of c12/c13/c14. The highest level is in general the step where a flux is added to a pool. In this case this is the biomass increment ('bm_inc'). I.e., only 'bm_inc' should be defined as 'carbon' and its c13/c14 signature should be explicitly defined.


### Units

All plant-related pools (pleaf, proot, psapw, pwood, plabl) are in units of gC/ind. and gN/ind. 'nind' is the number of individuals per m2, [nind]=ind./m2. Other pools (plitt, psoil, pexud, inorganic N pools) are in units of gC/m2 and/or gN/m2. Note that C and N transfers from plant to litter have to be multiplied by 'nind' (number of individuals). Fluxes (dgpp, dnpp, dexud, dnup, ...) are in units of gC/m2/day or gN/m2/day. 


### Year-to-year memory

In general, all pool variables carry on information from year to year, while all others don't. Exceptions are previous year's values to define "buffers", e.g., where moving average of previous N days is used to calculate stuff. This is the case for soil temperature (soiltemp_sitch), where the previous year's daily temperature field is used. This is defined within biosphere, using Fortran's 'SAVE' statement, and updated at the end of each year. For the first simulation year, the "previous" year's values are taken as the present year's values.


### Output

One output file is written for each output variable. Output variables have an according time dimension (daily and monthly output), and a space dimension (jpngr) and are kept separate from model state variables. To add a (module-independent) variable to be written to the output, add statements in the following files:
- outvars.mod.F: Declare non-module specific output variables.
- init.F: SR initoutput: Initialise output array, called daily.
- initio.F: Open file for output writing.
- getout.F: Copy daily updated state variables to according position in output array.
- writeout_ascii.F: Write output array into file opened by initio.F.
For module-specific output variables, add according statements contained in SR (named as 'initio_<modulename>') contained within respective module.


### Processes

*Phenology*
Temperature driven phenology ('dtphen') is determined by smoothed daily temperature (smoothing by interpolating from monthly mean air temperatures). For grasses, 'dtphen' is zero if smoothed temperature is <5 deg C and 1 if it is >5 deg C. When dtphen=1, then grass grows daily with LAI dynamically developping over the course of the season.




