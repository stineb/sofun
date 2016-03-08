module _params_core
  !////////////////////////////////////////////////////////////////
  !  This module contains parameters that are not modified, but needed
  !  to define variables, dimension lengths, etc.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  implicit none

  integer, parameter :: recycle = 5              ! length of standard recycling period
  integer, parameter :: ndayyear = 365           ! number of days in a year
  integer, parameter :: nmonth = 12              ! number of months in a year
  integer, parameter :: maxlen_runname = 30      ! maximum length of runname (arbitrary)
  integer, parameter :: maxgrid = 1              ! number of spatial gridcells (dummy dimension for later code extension)
  integer, parameter :: nbucket = 2              ! number of buckets for soil water model
  integer, parameter :: npft = 1                 ! number of PFTs
  integer, parameter :: nlu = 1                  ! number of land units (tiles)
  integer, parameter :: lunat = 1                ! ID of natural land unit
  integer, parameter :: lucrop = 2               ! ID of crop land unit

  integer, parameter, dimension(npft) :: pft_start = 1
  integer, parameter, dimension(npft) :: pft_end   = 1

  integer, parameter, dimension(nmonth)   :: ndaymonth = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! number of days per month
  integer, parameter, dimension(nmonth+1) :: middaymonth = (/16,44,75,105,136,166,197,228,258,289,319,350,381/) ! day of year of middle-month-day
  integer, parameter, dimension(nmonth)   :: cumdaymonth = (/31,59,90,120,151,181,212,243,273,304,334,365/)

  real, parameter    :: pi = 3.14159265359       ! pi - what else?
  real, parameter    :: dummy = -9999.0          ! arbitrary dummy value

  real, parameter :: c_molmass = 12.0107         ! g C / mol C
  real, parameter :: n_molmass = 14.0067         ! g N / mol N
  real, parameter :: h2o_molmass = 44.013        ! g H2O / mol H2O
  real, parameter :: c_content_of_biomass = 0.46 ! gC / g-dry mass


end module _params_core

