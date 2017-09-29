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

    ! !////////////////////////////////////////////////////////////////
    ! ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    ! !----------------------------------------------------------------
    ! if (interface%params_siml%loutforcing) then

    !   ! ANNUAL MEAN TEMPERATURE (DEG C) 
    !   filnam=trim(prefix)//'.a.temp.out'
    !   open(951,file=filnam,err=999,status='unknown')

    !   ! ANNUAL TOTAL N INPUT 
    !   filnam=trim(prefix)//'.a.nin.out'
    !   open(952,file=filnam,err=999,status='unknown')

    ! end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_forcing


  subroutine initio_nc_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    integer, parameter :: ndims = 2

    integer :: ncid
    integer :: londimid, latdimid
    integer :: dimids(ndims)
    ! integer :: start(ndims), count(ndims)
    integer :: varid_temp

    integer :: jpngr

    ! real, dimension(:,:), allocatable :: dtemp_arr

    ! xxx test
    integer :: doy = 1

    print*,'in initio_nc_forcing...'

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    filnam=trim(prefix)//'.d.temp.nc'
    call check( nf90_create( trim(filnam), NF90_CLOBBER, ncid ) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim( ncid, "lon", size(interface%domaininfo%lon), londimid ) )
    call check( nf90_def_dim( ncid, "lat", size(interface%domaininfo%lat), latdimid ) )

    print*,'3'

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/ latdimid, londimid /)

    ! Define the variable. The type of the variable in this case is
    ! NF90_DOUBLE.
    call check( nf90_def_var( ncid, "dtemp", NF90_FLOAT, dimids, varid_temp ) )

    ! End define mode. This tells netCDF we are done defining metadata.
    call check( nf90_enddef( ncid ) )

    print*,'4'

    ! xxx test
    doy = 1

    ! Write the data, gridcell by gridcell
    do jpngr=1,size(interface%grid)
      if (interface%grid(jpngr)%dogridcell) then
        call check( nf90_put_var( ncid, varid_temp, interface%climate(jpngr)%dtemp(1), start = (/ interface%grid(jpngr)%ilat, interface%grid(jpngr)%ilon /) ) )
      end if
    end do

    print*,'5'

    ! ! Write the data to the file. Although netCDF supports
    ! ! reading and writing subsets of data, in this case we write all the
    ! ! data in one operation.
    ! call check( nf90_put_var( ncid, varid_temp, dtemp_arr ) )

    print*,'6'

    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close( ncid ) )
    print*,'6'

    stop 'check the nc file'

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

    ! !-------------------------------------------------------------------------
    ! ! ANNUAL OUTPUT
    ! ! Write annual value, summed over all PFTs / LUs
    ! ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    ! !-------------------------------------------------------------------------
    ! if (interface%params_siml%loutforcing) then

    !   itime = real(interface%steering%outyear)

    !   ! write(951,999) itime, outatemp(jpngr)
    !   ! write(952,999) itime, sum(outanin(:,jpngr))

    ! end if

    return

    999     format (F20.8,F20.8)

  end subroutine writeout_ascii_forcing


  subroutine writeout_nc_forcing
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    print*,'writeout_nc_forcing: doing nothing here'

  end subroutine writeout_nc_forcing


  subroutine check( status )
    !/////////////////////////////////////////////////////////////////////////
    ! Auxiliary subroutine handling NetCDF 
    !-------------------------------------------------------------------------
    use netcdf
    integer, intent (in) :: status
    if ( status /= nf90_noerr ) then 
      print *, trim( nf90_strerror(status) )
      stop "Stopped"
    end if
  end subroutine check  

end module md_interface
