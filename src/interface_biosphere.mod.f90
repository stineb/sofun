module md_interface

  use md_params_core, only: maxgrid, nlu, npft, ndayyear, dummy
  use md_grid, only: gridtype, domaininfo_type
  use md_forcing, only: landuse_type, climate_type, ninput_type
  use md_params_domain, only: type_params_domain
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: outtype_steering, paramstype_siml

  implicit none

  private
  public interfacetype_biosphere, interface, initoutput_forcing, &
    initio_nc_forcing, getout_daily_forcing, writeout_nc_forcing, &
    outtype_biosphere

  type paramstype_calib
    ! real :: k_decay_tissue
    real :: kphio
    real :: soilm_par_a
    real :: soilm_par_b
    real :: vpdstress_par_a
    real :: vpdstress_par_b
    real :: vpdstress_par_m
  end type paramstype_calib  

  type interfacetype_biosphere
    integer                                             :: year
    real                                                :: pco2
    type( gridtype )      , dimension(:),   allocatable :: grid
    type( paramtype_soil ), dimension(:),   allocatable :: soilparams
    type( landuse_type )  , dimension(:),   allocatable :: landuse
    type( climate_type )  , dimension(:),   allocatable :: climate
    type( ninput_type)    , dimension(:),   allocatable :: ninput_field
    real                  , dimension(:,:), allocatable :: dfapar_field
    type( domaininfo_type )                             :: domaininfo
    type( outtype_steering )                            :: steering
    type( paramstype_siml )                             :: params_siml
    real, dimension(:,:), allocatable                   :: fpc_grid
    type( paramstype_calib )                            :: params_calib    ! calibratable parameters
  end type interfacetype_biosphere

  !----------------------------------------------------------------
  ! Return variable of biosphere()
  !----------------------------------------------------------------
  ! This is the derived type-return variable of the function biosphere(),
  ! holding variables used for the cost function in sofun_calib.f90
  type outtype_biosphere
    real, dimension(ndayyear) :: gpp
    real, dimension(ndayyear) :: fapar
    real, dimension(ndayyear) :: transp
    real, dimension(ndayyear) :: latenth
  end type outtype_biosphere

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
  real, allocatable, dimension(:,:) :: outdfapar

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_fland
  character(len=256) :: ncoutfilnam_vegcover
  character(len=256) :: ncoutfilnam_temp
  character(len=256) :: ncoutfilnam_fapar

  character(len=*), parameter :: FLAND_NAME="fland"
  character(len=*), parameter :: VEGCOVER_NAME="vegcover"
  character(len=*), parameter :: TEMP_NAME="temp"
  character(len=*), parameter :: FAPAR_NAME="fapar"

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

    if (interface%params_siml%loutforcing) then

      ! Allocate memory for daily output variables
      if (interface%steering%init .and. interface%params_siml%loutdtemp ) allocate( outdtemp (interface%params_siml%outnt,ngridcells) )
      if (interface%steering%init .and. interface%params_siml%loutdfapar) allocate( outdfapar(interface%params_siml%outnt,ngridcells) )
      if (interface%params_siml%loutdtemp ) outdtemp (:,:) = 0.0
      if (interface%params_siml%loutdfapar) outdfapar(:,:) = 0.0

    end if

  end subroutine initoutput_forcing


  subroutine initio_nc_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_2D, init_nc_3D_time, init_nc_3D_pft, check

    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_interface"
    character(len=4) :: year_char

    integer :: jpngr, doy

    write(year_char,999) interface%steering%outyear

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    !----------------------------------------------------------------
    ! Land fraction
    !----------------------------------------------------------------
    if ( interface%steering%init ) then
      ncoutfilnam_fland = trim(prefix)//".fland.nc"
      print*,'initialising ', trim(ncoutfilnam_fland), '...'
      call init_nc_2D(  filnam   = trim(ncoutfilnam_fland), &
                        nlon     = interface%domaininfo%nlon, &
                        nlat     = interface%domaininfo%nlat, &
                        lon      = interface%domaininfo%lon, &
                        lat      = interface%domaininfo%lat, &
                        varnam   = FLAND_NAME, &
                        varunits = "unitless", &
                        longnam  = "gridcell fraction covered by land", &
                        title    = TITLE &
                        )
    end if

    !----------------------------------------------------------------
    ! Vegetation cover
    !----------------------------------------------------------------
    if ( interface%steering%init ) then
      ncoutfilnam_vegcover = trim(prefix)//".vegcover.nc"
      print*,'initialising ', trim(ncoutfilnam_vegcover), '...'
      call init_nc_3D_pft(  filnam   = trim(ncoutfilnam_vegcover), &
                            nlon     = interface%domaininfo%nlon, &
                            nlat     = interface%domaininfo%nlat, &
                            lon      = interface%domaininfo%lon, &
                            lat      = interface%domaininfo%lat, &
                            outnz    = npft, &
                            varnam   = VEGCOVER_NAME, &
                            varunits = "unitless", &
                            longnam  = "gridcell fraction covered by vegetation type", &
                            title    = TITLE &
                            )
    end if

    if ( .not. interface%steering%spinup &
         .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
         .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      !----------------------------------------------------------------
      ! Temperature output file 
      !----------------------------------------------------------------      
      if (interface%params_siml%loutdtemp) then
        ncoutfilnam_temp = trim(prefix)//'.'//year_char//".d.temp.nc"
        print*,'initialising ', trim(ncoutfilnam_temp), '...'
        call init_nc_3D_time( filnam   = trim(ncoutfilnam_temp), &
                              nlon     = interface%domaininfo%nlon, &
                              nlat     = interface%domaininfo%nlat, &
                              lon      = interface%domaininfo%lon, &
                              lat      = interface%domaininfo%lat, &
                              outyear  = interface%steering%outyear, &
                              outdt    = interface%params_siml%outdt, &
                              outnt    = interface%params_siml%outnt, &
                              varnam   = TEMP_NAME, &
                              varunits = "degrees Celsius", &
                              longnam  = "daily average 2 m temperature", &
                              title    = TITLE &
                              )
      end if

      !----------------------------------------------------------------
      ! fAPAR output file 
      !----------------------------------------------------------------
      if ( interface%params_siml%loutdfapar) then
        ncoutfilnam_fapar = trim(prefix)//'.'//year_char//".d.fapar.nc"
        print*,'initialising ', trim(ncoutfilnam_fapar), '...'
        call init_nc_3D_time( filnam   = trim(ncoutfilnam_fapar), &
                              nlon     = interface%domaininfo%nlon, &
                              nlat     = interface%domaininfo%nlat, &
                              lon      = interface%domaininfo%lon, &
                              lat      = interface%domaininfo%lat, &
                              outyear  = interface%steering%outyear, &
                              outdt    = interface%params_siml%outdt, &
                              outnt    = interface%params_siml%outnt, &
                              varnam   = FAPAR_NAME, &
                              varunits = "unitless", &
                              longnam  = "fraction of absorbed photosynthetically active radiation", &
                              title    = TITLE &
                              )
      end if

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

    ! local variables
    integer :: it

    if (interface%params_siml%loutforcing) then
      !----------------------------------------------------------------
      ! DAILY
      ! Collect daily output variables
      ! so far not implemented for isotopes
      !----------------------------------------------------------------
      it = floor( real( doy - 1 ) / real( interface%params_siml%outdt ) ) + 1
      if (interface%params_siml%loutdtemp)  outdtemp (it,jpngr) = outdtemp (it,jpngr) + interface%climate(jpngr)%dtemp(doy) / real( interface%params_siml%outdt )
      if (interface%params_siml%loutdfapar) outdfapar(it,jpngr) = outdfapar(it,jpngr) + interface%dfapar_field(doy,jpngr)   / real( interface%params_siml%outdt )

      ! !----------------------------------------------------------------
      ! ! ANNUAL SUM OVER DAILY VALUES
      ! ! Collect annual output variables as sum of daily values
      ! !----------------------------------------------------------------
      ! if (interface%params_siml%loutforcing) then
      !   outatemp(jpngr)  = outatemp(jpngr)  + interface%climate(jpngr)%dtemp(doy) / ndayyear
      !   outanin(:,jpngr) = outanin(:,jpngr) + interface%ninput_field(jpngr)%dtot(doy)
      ! end if
    end if

  end subroutine getout_daily_forcing


  subroutine writeout_nc_forcing()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_2D, write_nc_3D_time, write_nc_3D_pft, check

    !-------------------------------------------------------------------------
    ! land fraction
    !-------------------------------------------------------------------------
    if (interface%steering%init) then
      print*,'writing ', trim(ncoutfilnam_fland), '...'
      call write_nc_2D( trim(ncoutfilnam_fland), &
                        FLAND_NAME, &
                        interface%domaininfo%maxgrid, &
                        interface%domaininfo%nlon, &
                        interface%domaininfo%nlat, &
                        interface%grid(:)%ilon, &
                        interface%grid(:)%ilat, &
                        interface%grid(:)%dogridcell, &
                        interface%grid(:)%landfrac &
                        )
    end if

    !-------------------------------------------------------------------------
    ! vegetation cover
    !-------------------------------------------------------------------------
    if (interface%steering%init) then
      print*,'writing ', trim(ncoutfilnam_vegcover), '...'
      call write_nc_3D_pft( trim(ncoutfilnam_vegcover), &
                            VEGCOVER_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            npft, &
                            interface%grid(:)%dogridcell, &
                            interface%fpc_grid(:,:) &
                            )
    end if

    if ( .not. interface%steering%spinup &
         .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
         .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      if (interface%params_siml%loutforcing) then
        !-------------------------------------------------------------------------
        ! dtemp
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdtemp) print*,'writing ', trim(ncoutfilnam_temp), '...'
        if (interface%params_siml%loutdtemp) call write_nc_3D_time( trim(ncoutfilnam_temp), &
                                                                      TEMP_NAME, &
                                                                      interface%domaininfo%maxgrid, &
                                                                      interface%domaininfo%nlon, &
                                                                      interface%domaininfo%nlat, &
                                                                      interface%grid(:)%ilon, &
                                                                      interface%grid(:)%ilat, &
                                                                      interface%params_siml%outnt, &
                                                                      interface%grid(:)%dogridcell, &
                                                                      outdtemp(:,:) &
                                                                      )


        !-------------------------------------------------------------------------
        ! fapar
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdfapar) print*,'writing ', trim(ncoutfilnam_fapar), '...'
        if (interface%params_siml%loutdfapar) call write_nc_3D_time(trim(ncoutfilnam_fapar), &
                                                                      FAPAR_NAME, &
                                                                      interface%domaininfo%maxgrid, &
                                                                      interface%domaininfo%nlon, &
                                                                      interface%domaininfo%nlat, &
                                                                      interface%grid(:)%ilon, &
                                                                      interface%grid(:)%ilat, &
                                                                      interface%params_siml%outnt, &
                                                                      interface%grid(:)%dogridcell, &
                                                                      outdfapar(:,:) &
                                                                      )

      end if
    end if

  end subroutine writeout_nc_forcing

end module md_interface
