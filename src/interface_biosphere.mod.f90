module md_interface

  use md_params_core, only: maxgrid, nlu, ndayyear, dummy
  use md_grid, only: gridtype, domaininfo_type
  use md_forcing, only: landuse_type, climate_type, ninput_type
  use md_params_domain, only: type_params_domain
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: outtype_steering, paramstype_siml

  implicit none

  private
  public interfacetype_biosphere, interface, initoutput_forcing, initio_forcing, &
    initio_nc_forcing, getout_daily_forcing, writeout_ascii_forcing, writeout_nc_forcing

  type interfacetype_biosphere
    integer                                             :: year
    real                                                :: pco2
    type( gridtype )      , dimension(:),   allocatable :: grid
    type( paramtype_soil ), dimension(:),   allocatable :: soilparams
    type( landuse_type)   , dimension(:),   allocatable :: landuse
    type( climate_type )  , dimension(:),   allocatable :: climate
    type( ninput_type)    , dimension(:),   allocatable :: ninput_field
    real                  , dimension(:,:), allocatable :: dfapar_field
    type( domaininfo_type )                             :: domaininfo
    type( outtype_steering )                            :: steering
    type( paramstype_siml )                             :: params_siml
  end type interfacetype_biosphere

  !----------------------------------------------------------------
  ! Interface instance is created here 
  ! (instead of locally defined and passed on as argument. Both are 
  ! ok but this has the advantage that unknown-size arguments are
  ! avoided).
  !----------------------------------------------------------------
  type( interfacetype_biosphere ) :: interface

  !----------------------------------------------------------------
  ! Module-specific daily output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:) :: outdtemp

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_temp
  character(len=*), parameter :: TEMP_NAME="temp"

  ! !----------------------------------------------------------------
  ! ! Module-specific annual output variables
  ! !----------------------------------------------------------------
  ! real, dimension(maxgrid)     :: outatemp
  ! real, dimension(nlu,maxgrid) :: outanin


contains

  subroutine initoutput_forcing( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    integer, intent(in) :: ngridcells

    ! Allocate memory for daily output variables
    if ( interface%steering%init .and. interface%params_siml%loutdtemp ) allocate( outdtemp(ndayyear,ngridcells) )
    if ( interface%params_siml%loutdtemp ) outdtemp(:,:) = 0.0

    ! ! annual output variables
    ! if (interface%params_siml%loutforcing) then
    !   outatemp(:) = 0.0
    !   outanin (:,:) = 0.0
    ! end if

  end subroutine initoutput_forcing


  subroutine initio_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens ascii output files.
    !----------------------------------------------------------------
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

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_forcing


  subroutine initio_nc_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D, check

    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: DOY_NAME  = "doy"
    character(len=*), parameter :: YEAR_NAME = "year"
    character(len=*), parameter :: filnamend = ".d.temp.nc"
    character(len=*), parameter :: varunits  = "degrees Celsius"
    character(len=*), parameter :: longnam   = "daily average 2 m temperature"
    character(len=*), parameter :: title     = "SOFUN GP-model output, module md_interface"
    character(len=4) :: year_char

    integer :: jpngr, doy
    integer, dimension(ndayyear) :: doy_vals

    write(year_char,999) interface%steering%outyear

    doy_vals = (/ (doy, doy = 1, ndayyear) /)

    if (interface%params_siml%lncoutdtemp) then

      prefix = "./output_nc/"//trim(interface%params_siml%runname)

      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      ncoutfilnam_temp = trim(prefix)//'.'//year_char//filnamend
      call init_nc_3D( filnam  = ncoutfilnam_temp, &
                      nlon     = interface%domaininfo%nlon, &
                      nlat     = interface%domaininfo%nlat, &
                      nz       = ndayyear, &
                      lon      = interface%domaininfo%lon, &
                      lat      = interface%domaininfo%lat, &
                      zvals    = doy_vals, &
                      recvals  = interface%steering%outyear, &
                      znam     = DOY_NAME, &
                      recnam   = YEAR_NAME, &
                      varnam   = TEMP_NAME, &
                      varunits = varunits, &
                      longnam  = longnam, &
                      title    = title &
                      )

    end if

    999  format (I4.4)

  end subroutine initio_nc_forcing


  subroutine getout_daily_forcing( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutdtemp) outdtemp(doy,jpngr) = interface%climate(jpngr)%dtemp(doy)

    ! !----------------------------------------------------------------
    ! ! ANNUAL SUM OVER DAILY VALUES
    ! ! Collect annual output variables as sum of daily values
    ! !----------------------------------------------------------------
    ! if (interface%params_siml%loutforcing) then
    !   outatemp(jpngr)  = outatemp(jpngr)  + interface%climate(jpngr)%dtemp(doy) / ndayyear
    !   outanin(:,jpngr) = outanin(:,jpngr) + interface%ninput_field(jpngr)%dtot(doy)
    ! end if

  end subroutine getout_daily_forcing


  subroutine writeout_ascii_forcing()
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    !-------------------------------------------------------------------------
    ! use md_params_siml, only: spinup, interface%params_siml%daily_out_startyr, &
    use md_params_core, only: ndayyear

    ! local variables
    real :: itime
    integer :: doy, moy, jpngr
    real, dimension(ndayyear) :: outdtemp_tot

    outdtemp_tot(:) = 0.0

    if (nlu>1) stop 'Output only for one LU category implemented.'

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutdtemp) then
      ! if ( .not. interface%steering%spinup &
      !   .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
      !   .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

        ! Write daily output only during transient simulation
        do doy=1,ndayyear

          ! Get weighted average
          do jpngr=1,size(interface%grid)
            outdtemp_tot(doy) = outdtemp_tot(doy) + outdtemp(doy,jpngr) * interface%grid(jpngr)%landfrac * interface%grid(jpngr)%area
          end do
          outdtemp_tot(doy) = outdtemp_tot(doy) / interface%domaininfo%landarea

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real( interface%steering%outyear ) + real( doy - 1 ) / real( ndayyear )
          
          write(950,999) itime, outdtemp_tot(doy)

        end do
      ! end if
    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_forcing


  subroutine writeout_nc_forcing()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_3D, check

    ! local variables
    integer :: doy, jpngr
    integer :: ncid
    integer :: varid_temp

    real, dimension(:,:,:,:), allocatable :: outarr

    ! if ( .not. interface%steering%spinup &
    !       .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
    !       .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      if (interface%params_siml%lncoutdtemp) then

        allocate( outarr(interface%domaininfo%nlon,interface%domaininfo%nlat,ndayyear,1) )
        outarr(:,:,:,:) = dummy        

        ! Populate output array
        do jpngr=1,size(interface%grid)
          if (interface%grid(jpngr)%dogridcell) then

            ! do doy=1,ndayyear
              ! option A
              ! call check( nf90_put_var( ncid, varid_temp, interface%climate(jpngr)%dtemp(1), start = (/ interface%grid(jpngr)%ilat, interface%grid(jpngr)%ilon /) ) )

              ! option B
              ! call check( nf90_put_var( ncid, varid_temp, outdtemp(doy,jpngr), start = (/ doy, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat /) ) )            
            ! end do
            ! call check( nf90_put_var( ncid, varid_temp, outdtemp(:,jpngr), start = (/ interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, 1 /), count = (/ 1, 1, ndayyear /) ) )

            ! populate array
            outarr(interface%grid(jpngr)%ilon,interface%grid(jpngr)%ilat,:,1) = outdtemp(:,jpngr)

          end if
        end do

        call write_nc_3D( ncoutfilnam_temp, TEMP_NAME, interface%domaininfo%nlon, interface%domaininfo%nlat, ndayyear, outarr(:,:,:,:)  )

        ! deallocate memory
        deallocate( outarr )

      end if

    ! end if

  end subroutine writeout_nc_forcing

end module md_interface
