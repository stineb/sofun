module md_grid
  !////////////////////////////////////////////////////////////////
  ! Module for handling variables defining the model grid. For site-
  ! scale simulation this is trivial but still necesssary.
  !
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  implicit none

  private
  public gridtype !, getgrid, domaininfo_type, get_domaininfo, type_params_domain

  type gridtype
    integer :: ilon
    integer :: ilat
    real :: lon
    real :: lat
    real :: elv
    real :: landfrac
    real :: area
    logical :: dogridcell
    real :: nu               ! true anomaly (orbital parameter), recalculated each year for each gridcell in solar()
    real :: lambda           ! true longitude (orbital parameter), recalculated each year for each gridcell in solar()
  end type gridtype

end module md_grid

