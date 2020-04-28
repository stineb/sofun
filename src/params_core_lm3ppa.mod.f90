module md_params_core
  !////////////////////////////////////////////////////////////////
  ! This module contains parameters that are not modified, but needed
  ! to define variables, dimension lengths, etc.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  implicit none

  integer :: ntstepsyear                         ! 365 when daily
  integer, parameter :: ndayyear = 365           ! number of days in a year
  
  integer, parameter :: nmonth = 12              ! number of months in a year

  ! From LM3-PPA
  integer, parameter :: nlayers_soil = 3         ! number of soil layers
  integer, parameter :: n_dim_soil_types = 9     ! number of soil types
  integer, parameter :: MSPECIES = 15            ! number of species
  integer, parameter :: MAX_INIT_COHORTS = 10    ! Number of initial cohorts
  integer, parameter :: out_max_cohorts = 50     ! Try: Number of maximum cohorts

  integer, parameter :: nvars_hourly_tile = 15
  integer, parameter :: nvars_daily_tile = 35
  integer, parameter :: nvars_daily_cohorts = 27
  integer, parameter :: nvars_annual_tile = 44
  integer, parameter :: nvars_annual_cohorts = 23

  ! needed here
  real, parameter :: dummy = -9999.0             ! arbitrary dummy value

  integer, parameter, dimension(nmonth)   :: ndaymonth = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! number of days per month

end module md_params_core