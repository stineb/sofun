module md_interface

  use md_params_core, only: maxgrid, nmonth, nlu
  use md_grid, only: gridtype
  use md_forcing_siterun, only: landuse_type, climate_type, ninput_type
  use md_params_site, only: paramstype_site
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: outtype_steering, paramstype_siml

  implicit none

  private
  public interfacetype_biosphere, interface, initoutput_forcing, initio_forcing, &
    getout_daily_forcing, writeout_ascii_forcing

  type interfacetype_biosphere
    integer                                           :: year
    real                                              :: pco2
    type( gridtype )      , dimension(maxgrid)        :: grid
    type( paramtype_soil ), dimension(maxgrid)        :: soilparams
    type( landuse_type)   , dimension(maxgrid)        :: landuse
    type( climate_type )  , dimension(maxgrid)        :: climate
    type( ninput_type)    , dimension(maxgrid)        :: ninput_field
    real                  , dimension(nmonth,maxgrid) :: mfapar_field
    type( paramstype_site )                           :: params_site
    type( outtype_steering )                          :: steering
    type( paramstype_siml )                           :: params_siml
  end type interfacetype_biosphere

  type( interfacetype_biosphere ) :: interface

  !----------------------------------------------------------------
  ! Module-specific daily output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:) :: outdtemp

  !----------------------------------------------------------------
  ! Module-specific annual output variables
  !----------------------------------------------------------------
  real, dimension(maxgrid)     :: outatemp
  real, dimension(nlu,maxgrid) :: outanin


contains

  subroutine initoutput_forcing( interface )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, maxgrid

    ! argument
    type( interfacetype_biosphere ), intent(in) :: interface

    ! Allocate memory for daily output variables
    if (interface%steering%init .and. interface%params_siml%loutdtemp  ) allocate( outdtemp(ndayyear,maxgrid) )
    outdtemp(:,:) = 0.0

    ! annual output variables
    if (interface%params_siml%loutforcing) then
      outatemp(:) = 0.0
      outanin (:,:) = 0.0
    end if

  end subroutine initoutput_forcing


  subroutine initio_forcing( interface )
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    ! argument
    type( interfacetype_biosphere ), intent(in) :: interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------
    ! DAILY MEAN TEMPERATURE (DEG C)
    if (interface%params_siml%loutdtemp) then
      filnam=trim(prefix)//'.d.temp.out'
      open(950,file=filnam,err=999,status='unknown')
    end if 

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutforcing) then

      ! ANNUAL MEAN TEMPERATURE (DEG C) 
      filnam=trim(prefix)//'.a.temp.out'
      open(951,file=filnam,err=999,status='unknown')

      ! ANNUAL TOTAL N INPUT 
      filnam=trim(prefix)//'.a.nin.out'
      open(952,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_forcing


  subroutine getout_daily_forcing( interface, jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft

    ! arguments
    type( interfacetype_biosphere ), intent(in) :: interface
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutdtemp) outdtemp(doy,jpngr) = interface%climate(jpngr)%dtemp(doy)

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    if (interface%params_siml%loutforcing) then
      outatemp(jpngr)  = outatemp(jpngr)  + interface%climate(jpngr)%dtemp(doy) / ndayyear
      outanin(:,jpngr) = outanin(:,jpngr) + interface%ninput_field(jpngr)%dtot(doy)
    end if

  end subroutine getout_daily_forcing


  subroutine writeout_ascii_forcing( interface )
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    type( interfacetype_biosphere ), intent(in) :: interface

    ! local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    if (nlu>1) stop 'Output only for one LU category implemented.'

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if ( .not. interface%steering%spinup &
      .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
      .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      ! Write daily output only during transient simulation
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(interface%steering%outyear) + real(day-1)/real(ndayyear)
        
        if (interface%params_siml%loutdtemp ) write(950,999) itime, outdtemp(day,jpngr)

      end do
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutforcing) then

      itime = real(interface%steering%outyear)

      write(951,999) itime, outatemp(jpngr)
      write(952,999) itime, sum(outanin(:,jpngr))

    end if

    return

    999     format (F20.8,F20.8)

  end subroutine writeout_ascii_forcing


end module md_interface