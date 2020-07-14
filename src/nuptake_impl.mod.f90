module md_nuptake_impl
  !////////////////////////////////////////////////////////////////
  ! Contains functions for using empirical relationships to model
  ! N uptake, driven by GPP 
  !--------------------------------------------------------------
  ! load modules
  use md_params_core
  use md_tile, only: tile_type, tile_fluxes_type
  use md_interface, only: interface
  use md_grid, only: gridtype, domaininfo_type

  implicit none
  private
  public nuptake_impl, get_preds_nimpl, getpar_nimpl, initio_nc_nimpl, initoutput_nimpl, &
    getout_annual_nimpl, writeout_nc_nimpl

  !----------------------------------------------------------------
  ! Variables without memory (not necessarily just fluxes; just define the type) 
  !----------------------------------------------------------------
  type plant_nimpl_fluxes_type
    real :: anpp        ! annual total net primary production (gC m-2 yr-1)
    real :: aanpp       ! annual aboveground net primary production (gC m-2 yr-1)
    real :: abnpp       ! annual belowground net primary production (gC m-2 yr-1)
    real :: alnpp       ! annual net primary production for leaf production (gC m-2 yr-1)
    real :: awnpp       ! annual net primary production for wood production (gC m-2 yr-1)
    real :: leafcn      ! the ratio of leaf carbon per mass to leaf nitrogen per mass (unitness)
    real :: alnf        ! annual laef nitrogen flux (gC m-2 yr-2)
    real :: awnf        ! annual wood nitrogen flux (gC m-2 yr-2)
    real :: abnf        ! annual belowground(root) nitrogen flux (gC m-2 yr-2)
  end type plant_nimpl_fluxes_type

  type tile_nimpl_fluxes_type
    type(plant_nimpl_fluxes_type), dimension(npft) :: plant
  end type tile_nimpl_fluxes_type

  !----------------------------------------------------------------
  ! object containing module-specific fluxes (create this new object)
  !----------------------------------------------------------------
  type(tile_nimpl_fluxes_type), dimension(nlu) :: tile_nimpl_fluxes

  !----------------------------------------------------------------
  ! Coefficients of statistical models
  !----------------------------------------------------------------
  type coef_nimpl_type

    ! BP:GPP model
    real :: cnsoil_bp
    real :: age_bp
    real :: fapar_bp
    real :: alpha_bp
    real :: intersect_bp

    ! ANPP:GPP model
    real :: cnsoil_anpp
    real :: age_anpp
    real :: fapar_anpp
    real :: alpha_anpp
    real :: intersect_anpp

    ! ALNPP:NPP model
    real :: ppfd_alnpp
    real :: tg_alnpp
    real :: vpd_alnpp
    real :: intersect_alnpp

    ! Leaf C:N model
    real :: vcmax25_leafcn
    real :: lma_leafcn
    real :: intersect_leafcn

    ! Constant ratio
    real :: root_cn
    real :: wood_cn

  end type coef_nimpl_type  

  type(coef_nimpl_type) :: coef_nimpl

  !----------------------------------------------------------------
  ! Predictor fields
  !----------------------------------------------------------------
  type preds_nimpl_type
    real :: cnsoil
    real :: age
    real :: fapar
    real :: alpha
    real :: ppfd
    real :: tg
    real :: vpd
    real :: lma
  end type preds_nimpl_type

  type(preds_nimpl_type), dimension(:), allocatable :: preds_nimpl

  !----------------------------------------------------------------
  ! Specify file and variable names for NetCDF reading (Vcmax25 may need additional input from SOFUN)
  !----------------------------------------------------------------
  character(len=100), parameter :: filnam_cnsoil = "./input/global/nimpl/CNrt.nc"
  character(len=100), parameter :: varnam_cnsoil = "CNrt"
  character(len=100), parameter :: filnam_age    = "./input/global/nimpl/age.nc"
  character(len=100), parameter :: varnam_age    = "age"
  character(len=100), parameter :: filnam_fapar  = "./input/global/nimpl/fAPAR.nc"
  character(len=100), parameter :: varnam_fapar  = "fAPAR"
  character(len=100), parameter :: filnam_alpha  = "./input/global/nimpl/alpha.nc"
  character(len=100), parameter :: varnam_alpha  = "alpha"
  character(len=100), parameter :: filnam_ppfd   = "./input/global/nimpl/PPFD.nc"
  character(len=100), parameter :: varnam_ppfd   = "PPFD"
  character(len=100), parameter :: filnam_Tg     = "./input/global/nimpl/Tg.nc"
  character(len=100), parameter :: varnam_Tg     = "Tg"
  character(len=100), parameter :: filnam_vpd    = "./input/global/nimpl/vpd.nc"
  character(len=100), parameter :: varnam_vpd    = "vpd"
  ! character(len=100), parameter :: filnam_vcmax25= "Vcmax25.nc"  ! needs further work
  ! character(len=100), parameter :: varnam_vcmax25= "Vcmax25"     ! needs further work
  character(len=100), parameter :: filnam_lma    = "./input/global/nimpl/LMA.nc"
  character(len=100), parameter :: varnam_lma    = "LMA"
  
  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! annual
  real, dimension(:), allocatable :: outanpp
  ! xxx complement

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_anpp
  ! xxx complement

  character(len=*), parameter :: NPP_NAME = "npp"
  ! xxx complement

contains

  subroutine nuptake_impl( jpngr, dogridcell, tile, tile_fluxes, init )
    !////////////////////////////////////////////////////////////////
    ! Determines all the downstream fluxes as a function of GPP
    ! using pre-defined statistical relationships
    !--------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    logical, intent(in) :: dogridcell
    type(tile_type), dimension(nlu), intent(in) :: tile
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    logical, intent(in) :: init

    ! local variable
    integer :: lu, pft

    !--------------------------------------------------------------
    ! Predict using statistical models
    !Note that some data was log-transfromed (cnsoil, age, PPFD, vpd, Vcmax25, lma) while some are not (alpha, fAPAR, Tg)
    !All ratios were using logit function (logit(y)=x), therefore we should convert them into y = 1/(1+ exp(-x))
    ! leaf c/n model were using log function, so it should be exp in advance
    !--------------------------------------------------------------
    ! Make tile_nimpl_fluxes a field
    if (dogridcell) then

      lu = 1
      do pft = 1,npft
        !print*,'1'
        tile_nimpl_fluxes(lu)%plant(pft)%anpp  = tile_fluxes(lu)%plant(pft)%agpp * (1/(1 + EXP(-(coef_nimpl%cnsoil_bp * LOG(preds_nimpl(jpngr)%cnsoil) + coef_nimpl%age_bp * LOG(preds_nimpl(jpngr)%age) + coef_nimpl%fapar_bp * preds_nimpl(jpngr)%fapar +coef_nimpl%alpha_bp * preds_nimpl(jpngr)%alpha + coef_nimpl%intersect_bp))))
        !print*,'2'
        tile_nimpl_fluxes(lu)%plant(pft)%aanpp = tile_fluxes(lu)%plant(pft)%agpp * (1/(1 + EXP(-(coef_nimpl%cnsoil_anpp * LOG(preds_nimpl(jpngr)%cnsoil) + coef_nimpl%age_anpp * LOG(preds_nimpl(jpngr)%age) + coef_nimpl%fapar_anpp * preds_nimpl(jpngr)%fapar +coef_nimpl%alpha_anpp * preds_nimpl(jpngr)%alpha + coef_nimpl%intersect_anpp))))
        !print*,'3'
        tile_nimpl_fluxes(lu)%plant(pft)%abnpp = tile_nimpl_fluxes(lu)%plant(pft)%anpp - tile_nimpl_fluxes(lu)%plant(pft)%aanpp
        !print*,'4'
        tile_nimpl_fluxes(lu)%plant(pft)%alnpp = tile_nimpl_fluxes(lu)%plant(pft)%aanpp * (1/(1+EXP(-(coef_nimpl%ppfd_alnpp * LOG(preds_nimpl(jpngr)%ppfd) + coef_nimpl%tg_alnpp * preds_nimpl(jpngr)%tg + coef_nimpl%vpd_alnpp * LOG(preds_nimpl(jpngr)%vpd) + coef_nimpl%intersect_alnpp))))
        !print*,'5'
        tile_nimpl_fluxes(lu)%plant(pft)%awnpp = tile_nimpl_fluxes(lu)%plant(pft)%aanpp - tile_nimpl_fluxes(lu)%plant(pft)%alnpp
        !print*,'tile_nimpl_fluxes(lu)%plant(pft)%awnpp', tile_nimpl_fluxes(lu)%plant(pft)%awnpp
        !print*,'6'
        if (tile(lu)%plant(pft)%vcmax25 > 0.0) then
          tile_nimpl_fluxes(lu)%plant(pft)%leafcn = EXP(coef_nimpl%vcmax25_leafcn * LOG(tile(lu)%plant(pft)%vcmax25) + coef_nimpl%lma_leafcn * LOG(preds_nimpl(jpngr)%lma) + coef_nimpl%intersect_leafcn)
        end if
        !print*,'7'
        if (tile_nimpl_fluxes(lu)%plant(pft)%leafcn > 0.0) then
          tile_nimpl_fluxes(lu)%plant(pft)%alnf = tile_nimpl_fluxes(lu)%plant(pft)%alnpp/tile_nimpl_fluxes(lu)%plant(pft)%leafcn 
        end if
        tile_nimpl_fluxes(lu)%plant(pft)%awnf = tile_nimpl_fluxes(lu)%plant(pft)%awnpp / coef_nimpl%wood_cn
        tile_nimpl_fluxes(lu)%plant(pft)%abnf = tile_nimpl_fluxes(lu)%plant(pft)%abnpp / coef_nimpl%root_cn
      end do

    else

      lu = 1
      !print*,'1'
      tile_nimpl_fluxes(lu)%plant(:)%anpp  = dummy
      !print*,'2'
      tile_nimpl_fluxes(lu)%plant(:)%aanpp = dummy
      !print*,'3'
      tile_nimpl_fluxes(lu)%plant(:)%abnpp = dummy
      !print*,'4'
      tile_nimpl_fluxes(lu)%plant(:)%alnpp = dummy
      !print*,'5'
      tile_nimpl_fluxes(lu)%plant(:)%awnpp = dummy
      !print*,'6'
      tile_nimpl_fluxes(lu)%plant(:)%leafcn = dummy
      !print*,'7'

      tile_nimpl_fluxes(lu)%plant(:)%alnf = dummy
      tile_nimpl_fluxes(lu)%plant(:)%awnf = dummy
      tile_nimpl_fluxes(lu)%plant(:)%abnf = dummy

    end if

  end subroutine nuptake_impl


  subroutine get_preds_nimpl( domaininfo, grid )
    !////////////////////////////////////////////////////////////////
    ! Some explanations XXX 
    !----------------------------------------------------------------
    ! arguments
    use md_params_core, only: npft

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid

    ! allocate memory
    allocate( preds_nimpl(domaininfo%maxgrid) )

    !--------------------------------------------------------------
    ! Read predictor fields from files, populates 'preds_nimpl'
    !--------------------------------------------------------------
    ! BP:GPP model
    call get_preds_nc_byvar( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
    call get_preds_nc_byvar( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
    call get_preds_nc_byvar( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
    call get_preds_nc_byvar( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

    ! ANPP:GPP model
    call get_preds_nc_byvar( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
    call get_preds_nc_byvar( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
    call get_preds_nc_byvar( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
    call get_preds_nc_byvar( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

    ! ALNPP:NPP model
    call get_preds_nc_byvar( trim(filnam_ppfd), trim(varnam_ppfd), domaininfo, grid, preds_nimpl(:)%ppfd )
    call get_preds_nc_byvar( trim(filnam_tg),   trim(varnam_tg),   domaininfo, grid, preds_nimpl(:)%tg )
    call get_preds_nc_byvar( trim(filnam_vpd),  trim(varnam_vpd),  domaininfo, grid, preds_nimpl(:)%vpd )

    ! Leaf C:N model
    ! call get_preds_nc_byvar( trim(filnam_vcmax25), trim(varnam_vcmax25), domaininfo, grid, preds_nimpl(:)%vcmax25 )
    call get_preds_nc_byvar( trim(filnam_lma),     trim(varnam_lma),     domaininfo, grid, preds_nimpl(:)%lma )

  end subroutine get_preds_nimpl


  subroutine getpar_nimpl()
    !////////////////////////////////////////////////////////////////
    !Extract all coeffients and constant in a dat file
    !--------------------------------------------------------------
    use md_sofunutils, only: getparreal
    ! BP/GPP model
    coef_nimpl%cnsoil_bp  = getparreal( 'params/params_nimpl.dat', 'cnsoil_bp' )
    coef_nimpl%age_bp     = getparreal( 'params/params_nimpl.dat', 'age_bp' )
    coef_nimpl%fapar_bp  = getparreal( 'params/params_nimpl.dat', 'fapar_bp' )
    coef_nimpl%alpha_bp     = getparreal( 'params/params_nimpl.dat', 'alpha_bp' )
    coef_nimpl%intersect_bp     = getparreal( 'params/params_nimpl.dat', 'intersect_bp' )

    ! ANPP/GPP model
    coef_nimpl%cnsoil_anpp  = getparreal( 'params/params_nimpl.dat', 'cnsoil_anpp' )
    coef_nimpl%age_anpp     = getparreal( 'params/params_nimpl.dat', 'age_anpp' )
    coef_nimpl%fapar_anpp  = getparreal( 'params/params_nimpl.dat', 'fapar_anpp' )
    coef_nimpl%alpha_anpp     = getparreal( 'params/params_nimpl.dat', 'alpha_anpp' )
    coef_nimpl%intersect_anpp     = getparreal( 'params/params_nimpl.dat', 'intersect_anpp' )

    ! ALNPP/NPP model
    coef_nimpl%ppfd_alnpp = getparreal( 'params/params_nimpl.dat', 'ppfd_alnpp' )
    coef_nimpl%tg_alnpp = getparreal( 'params/params_nimpl.dat', 'tg_alnpp' )
    coef_nimpl%vpd_alnpp = getparreal( 'params/params_nimpl.dat', 'vpd_alnpp' )
    coef_nimpl%intersect_alnpp     = getparreal( 'params/params_nimpl.dat', 'intersect_alnpp' )
    
    ! Leaf C:N model
    ! coef_nimpl%vcmax25_leafcn = getparreal( 'params/params_nimpl.dat', 'vcmax25_leafcn' )
    coef_nimpl%lma_leafcn = getparreal( 'params/params_nimpl.dat', 'lma_leafcn' )
    coef_nimpl%intersect_leafcn     = getparreal( 'params/params_nimpl.dat', 'intersect_leafcn' )

    ! Constant
    coef_nimpl%root_cn = getparreal( 'params/params_nimpl.dat', 'root_cn' )
    coef_nimpl%wood_cn = getparreal( 'params/params_nimpl.dat', 'wood_cn' )

  end subroutine getpar_nimpl


  subroutine get_preds_nc_byvar( filnam, varname, domaininfo, grid, pred )
    !////////////////////////////////////////////////////////////////
    ! xxx add explanation here
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: check

    ! arguments
    character(len=*) :: filnam
    character(len=*) :: varname
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    real, dimension(domaininfo%maxgrid), intent(out) :: pred

    ! local variables
    integer :: ncid, varid
    integer :: latdimid, londimid
    integer :: nlat_arr, nlon_arr
    real, allocatable, dimension(:)   :: lon_arr
    real, allocatable, dimension(:)   :: lat_arr
    real, allocatable, dimension(:,:) :: pred_arr

    integer :: i, pft, jpngr, ilon_arr, ilat_arr, n_noinfo
    integer, dimension(domaininfo%maxgrid) :: ilon
    integer, dimension(domaininfo%maxgrid) :: ilat
    real, allocatable, dimension(:) :: tmp
    real :: ncfillvalue
    real :: dlat, dlon
    character(len=3), parameter :: lonname = "lon"
    character(len=3), parameter :: latname = "lat"
    ! character(len=100), parameter :: dimname_pft = "z"
    ! character(len=100), parameter :: varname = "pftcover"
    ! character(len=100), parameter :: filnam = "./input/global/landcover/modis_landcover_halfdeg_2010_FILLED.nc"

    !----------------------------------------------------------------  
    ! Get vegetation cover information from file
    !----------------------------------------------------------------
    print*,'getting predictor from ', trim(filnam), ' ...'

    ! Read arrays of all months of current year from file  
    call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, trim(latname), latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, trim(lonname), londimid ) )

    ! Get latitude information: nlon
    call check( nf90_inquire_dimension( ncid, londimid, len = nlon_arr ) )

    ! for index association, get ilon and ilat vectors
    ! Allocate array sizes now knowing nlon and nlat 
    allocate( lon_arr(nlon_arr) )
    allocate( lat_arr(nlat_arr) )

    !! Get longitude and latitude values
    !call check( nf90_get_var( ncid, londimid, lon_arr ) )
    !call check( nf90_get_var( ncid, latdimid, lat_arr ) )

    ! xxx try:
    lon_arr = (/ (i, i = 1,nlon_arr) /)
    lon_arr = (lon_arr - 1) * 0.5 - 180.0 + 0.25 
    lat_arr = (/ (i, i = 1,nlat_arr) /)
    lat_arr = (lat_arr - 1) * 0.5 - 90.0 + 0.25 

    ! Check if the resolution of the climate input files is identical to the model grid resolution
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    !print*,'dlon', dlon
    !print*,'dlat', dlat

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of nimpl predictor file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of nimpl predictor file is not identical with model grid.'

    ! get index associations
    do jpngr=1,domaininfo%maxgrid
      ilon_arr = 1
      do while (grid(jpngr)%lon/=lon_arr(ilon_arr))
        ilon_arr = ilon_arr + 1
      end do
      ilon(jpngr) = ilon_arr

      ilat_arr = 1
      do while (grid(jpngr)%lat/=lat_arr(ilat_arr))
        ilat_arr = ilat_arr + 1
      end do
      ilat(jpngr) = ilat_arr
    end do

    ! allocate size of output array
    allocate( pred_arr(nlon_arr,nlat_arr) )
    ! allocate( tmp(npft_in))

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, trim(varname), varid ) )

    ! Read the array
    call check( nf90_get_var( ncid, varid, pred_arr, start=(/1, 1/), count=(/nlon_arr, nlat_arr/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define field    
    do jpngr=1,domaininfo%maxgrid
      
      pred(jpngr) = pred_arr(ilon(jpngr),ilat(jpngr))

      ! is data actually available?
      if (pred(jpngr)==ncfillvalue) then
        pred(jpngr) = dummy
        grid(jpngr)%dogridcell = .false.
      end if

    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( pred_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    return

  end subroutine get_preds_nc_byvar


  subroutine initio_nc_nimpl()
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific NetCDF output files.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D_time, check
    
    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_nimpl"
    character(len=12) :: coef_age_bp_char

    integer :: jpngr, doy

    character(len=4)  :: year_char
    write(year_char,999) interface%steering%outyear

    ! convert parameter values to charaters
    write(coef_age_bp_char,888) coef_nimpl%age_bp
    ! xxx complement. this is like as.character()

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    if ( .not. interface%steering%spinup ) then
      !----------------------------------------------------------------
      ! Annual NPP output file 
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_anpp = trim(prefix)//'.'//year_char//".a.npp.nc"
        print*,'initialising ', trim(ncoutfilnam_anpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_anpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = NPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual net primary productivivty", &
                          title    = TITLE, &
                          globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if

      ! xxx complement

    end if

    888  format (F12.6)  ! for numeric (decimals)
    999  format (I4.4)   ! for integers
    
  end subroutine initio_nc_nimpl


  subroutine initoutput_nimpl( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific output variables
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: ngridcells

    ! annual
    if (interface%params_siml%loutnimpl) then

      if (interface%steering%init) then
        allocate( outanpp(ngridcells) )
        ! xxx complement
      end if

      outanpp(:) = 0.0
      ! xxx complement
    
    end if

  end subroutine initoutput_nimpl


  subroutine getout_annual_nimpl( jpngr, tile )
    !////////////////////////////////////////////////////////////////
    ! Called once a year to gather annual output variables.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    type(tile_type), dimension(nlu), intent(inout) :: tile

    ! local variables
    integer :: pft, lu

    ! outanrlarea(jpngr) = anrlarea
    if (interface%params_siml%loutnimpl) then

      lu = 1
      if ( abs(sum(tile(lu)%plant(:)%fpc_grid) - 1.0) > eps ) stop 'getout_annual_nimpl(): fpc_grid does not sum up to 1.0'
      outanpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%anpp * tile(lu)%plant(:)%fpc_grid )
      
      ! xxx complement

    end if

  end subroutine getout_annual_nimpl


  subroutine writeout_nc_nimpl()
    !/////////////////////////////////////////////////////////////////////////
    ! Writes module-specific NetCDF output
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_2D, write_nc_3D_time, check
    
    if ( .not. interface%steering%spinup ) then
      !-------------------------------------------------------------------------
      ! Annual NPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_anpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_anpp), &
                                                              NPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outanpp(:) &
                                                              )

      ! xxx complement

    end if

  end subroutine writeout_nc_nimpl


end module md_nuptake_impl
