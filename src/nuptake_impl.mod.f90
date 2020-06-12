module md_nuptake_impl
  !////////////////////////////////////////////////////////////////
  ! Contains functions for using empirical relationships to model
  ! N uptake, driven by GPP 
  !--------------------------------------------------------------
  ! load modules
  use md_params_core
  use md_tile, only: tile_fluxes_type
  use md_interface, only: interface

  implicit none
  private
  public objectsxxx

  !----------------------------------------------------------------
  ! Variables without memory (not necessarily just fluxes; just define the type) 
  !----------------------------------------------------------------
  type canopy_nimpl_fluxes_type
    real :: anpp        ! annual total net primary production (gC m-2 yr-1)
    real :: aanpp       ! annual aboveground net primary production (gC m-2 yr-1)
    real :: abnpp       ! annual belowground net primary production (gC m-2 yr-1)
    real :: alnpp       ! annual net primary production for leaf production (gC m-2 yr-1)
    real :: awnpp       ! annual net primary production for wood production (gC m-2 yr-1)
    real :: leafcn      ! the ratio of leaf carbon per mass to leaf nitrogen per mass (unitness)
    real :: alnf        ! annual laef nitrogen flux (gC m-2 yr-2)
    real :: awnf        ! annual wood nitrogen flux (gC m-2 yr-2)
    real :: abnf        ! annual belowground(root) nitrogen flux (gC m-2 yr-2)
  end type canopy_nimpl_fluxes_type

  type tile_nimpl_fluxes_type
    type(canopy_nimpl_fluxes_type) :: canopy
  end type tile_nimpl_fluxes_type

  !----------------------------------------------------------------
  ! object containing module-specific fluxes (create this new object)
  !----------------------------------------------------------------
  type(tile_nimpl_fluxes_type), allocatable, dimension(:) :: tile_nimpl_fluxes

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
    real :: PPFD_alnpp
    real :: Tg_alnpp
    real :: vpd_alnpp
    real :: intersect_alnpp

    ! Leaf C:N model
    real :: Vcmax25_leafcn
    real :: LMA_leafcn
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
    ! BP:GPP model
    real :: cnsoil
    real :: age
    real :: fapar
    real :: alpha

    ! ANPP:GPP model
    real :: cnsoil
    real :: age
    real :: fapar
    real :: alpha

    ! ALNPP:NPP model
    real :: PPFD
    real :: Tg
    real :: vpd

    ! Leaf C:N model
    real :: Vcmax25
    real :: LMA

  end type preds_nimpl_type

  type(preds_nimpl_type), dimension(:), allocatable :: preds_nimpl

  !----------------------------------------------------------------
  ! Specify file and variable names for NetCDF reading xxx (variable name needs double check)
  !----------------------------------------------------------------
  character(len=100), parameter :: filnam_cnsoil = "soilcn.nc"
  character(len=100), parameter :: varnam_cnsoil = "actual_variable_name"
  character(len=100), parameter :: filnam_age    = "age.nc"
  character(len=100), parameter :: varnam_age    = "actual_variable_name"
  character(len=100), parameter :: filnam_fapar  = "fapar.nc"
  character(len=100), parameter :: varnam_fapar  = "actual_variable_name"
  character(len=100), parameter :: filnam_alpha  = "alpha.nc"
  character(len=100), parameter :: varnam_alpha  = "actual_variable_name"
  character(len=100), parameter :: filnam_PPFD   = "PPFD.nc"
  character(len=100), parameter :: varnam_PPFD   = "actual_variable_name"
  character(len=100), parameter :: filnam_Tg     = "Tg.nc"
  character(len=100), parameter :: varnam_Tg     = "actual_variable_name"
  character(len=100), parameter :: filnam_vpd    = "vpd.nc"
  character(len=100), parameter :: varnam_vpd    = "actual_variable_name"
  character(len=100), parameter :: filnam_Vcmax25= "Vcmax25.nc"
  character(len=100), parameter :: varnam_Vcmax25= "actual_variable_name"
  character(len=100), parameter :: filnam_LMA    = "LMA.nc"
  character(len=100), parameter :: varnam_LMA    = "actual_variable_name"
contains

  subroutine nuptake_impl( tile_fluxes, init )
    !////////////////////////////////////////////////////////////////
    ! Determines all the downstream fluxes as a function of GPP
    ! using pre-defined statistical relationships
    !--------------------------------------------------------------
    ! arguments
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    logical, intent(in) :: init

    if (init) then

      ! allocate memory
      allocate( preds_nimpl(size(interface%grid)) )
      allocate( tile_fluxes(nlu) )

      !--------------------------------------------------------------
      ! Read predictor fields from files, populates 'preds_nimpl'
      !--------------------------------------------------------------
      ! BP:GPP model
      call get_preds_nc( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
      call get_preds_nc( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
      call get_preds_nc( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
      call get_preds_nc( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

      ! ANPP:GPP model
      call get_preds_nc( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
      call get_preds_nc( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
      call get_preds_nc( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
      call get_preds_nc( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

      ! ALNPP:NPP model
      call get_preds_nc( trim(filnam_PPFD), trim(varnam_PPFD), domaininfo, grid, preds_nimpl(:)%PPFD )
      call get_preds_nc( trim(filnam_Tg),    trim(varnam_Tg),    domaininfo, grid, preds_nimpl(:)%Tg )
      call get_preds_nc( trim(filnam_vpd),  trim(varnam_vpd),  domaininfo, grid, preds_nimpl(:)%vpd )

      ! Leaf C:N model
      call get_preds_nc( trim(filnam_Vcmax25), trim(varnam_Vcmax25), domaininfo, grid, preds_nimpl(:)%Vcmax25 )
      call get_preds_nc( trim(filnam_LMA),    trim(varnam_LMA),    domaininfo, grid, preds_nimpl(:)%LMA )

      !--------------------------------------------------------------
      ! Read coefficients, populates 'coef_nimpl'
      !--------------------------------------------------------------
      ! BP:GPP model
      call getpar_nimpl()

    end if

    !--------------------------------------------------------------
    ! Predict using statistical models
    !Note that some data was log-transfromed (cnsoil, age, PPFD, vpd, Vcmax25, lma) while some are not (alpha, fAPAR, Tg)
    !All ratios were using logit function (logit(y)=x), therefore we should convert them into y = 1/(1+ exp(-x))
    ! leaf c/n model were using log function, so it should be exp in advance
    !--------------------------------------------------------------
    ! XXX loop or make tile_nimpl_fluxes a field?
    tile_nimpl_fluxes(lu)%canopy%anpp = tile_fluxes(lu)%canopy%agpp * (1/(1+EXP(-(coef_nimpl%cnsoil_bp * LOG(preds_nimpl(:)%cnsoil) + coef_nimpl%age_bp * LOG(preds_nimpl(:)%age) + coef_nimpl%fapar_bp * preds_nimpl(:)%fapar +coef_nimpl%alpha_bp * preds_nimpl(:)%alpha + coef_nimpl%intersect_bp))))
    tile_nimpl_fluxes(lu)%canopy%aanpp = tile_fluxes(lu)%canopy%agpp * (1/(1+EXP(-(coef_nimpl%cnsoil_anpp * LOG(preds_nimpl(:)%cnsoil) + coef_nimpl%age_anpp * LOG(preds_nimpl(:)%age) + coef_nimpl%fapar_anpp * preds_nimpl(:)%fapar +coef_nimpl%alpha_anpp * preds_nimpl(:)%alpha + coef_nimpl%intersect_anpp))))
    tile_nimpl_fluxes(lu)%canopy%abnpp = tile_nimpl_fluxes(lu)%canopy%anpp - tile_nimpl_fluxes(lu)%canopy%aanpp
    tile_nimpl_fluxes(lu)%canopy%alnpp = tile_nimpl_fluxes(lu)%canopy%aanpp * (1/(1+EXP(-(coef_nimpl%PPFD_alnpp * LOG(preds_nimpl(:)%PPFD) + coef_nimpl%Tg_alnpp * preds_nimpl(:)%Tg + coef_nimpl%vpd_alnpp * LOG(preds_nimpl(:)%vpd) + coef_nimpl%intersect_alnpp))))
    tile_nimpl_fluxes(lu)%canopy%awnpp = tile_nimpl_fluxes(lu)%canopy%aanpp - tile_nimpl_fluxes(lu)%canopy%alnpp
    tile_nimpl_fluxes(lu)%canopy%leafcn = EXP(coef_nimpl%Vcmax25_leafcn * LOG(preds_nimpl(:)%Vcmax25) + coef_nimpl%LMA_leafcn * LOG(preds_nimpl(:)%LMA) + coef_nimpl%intersect_leafcn)
    tile_nimpl_fluxes(lu)%canopy%alnf = tile_nimpl_fluxes(lu)%canopy%alnpp/tile_nimpl_fluxes(lu)%canopy%leafcn 
    tile_nimpl_fluxes(lu)%canopy%awnf = tile_nimpl_fluxes(lu)%canopy%awnpp/coef_nimpl%wood_cn
    tile_nimpl_fluxes(lu)%canopy%abnf = tile_nimpl_fluxes(lu)%canopy%abnpp/coef_nimpl%root_cn

  end subroutine nuptake_impl


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
    coef_nimpl%PPFD_alnpp = getparreal( 'params/params_nimpl.dat', 'PPFD_alnpp' )
    coef_nimpl%Tg_alnpp = getparreal( 'params/params_nimpl.dat', 'Tg_alnpp' )
    coef_nimpl%vpd_alnpp = getparreal( 'params/params_nimpl.dat', 'vpd_alnpp' )
    coef_nimpl%intersect_alnpp     = getparreal( 'params/params_nimpl.dat', 'intersect_alnpp' )
    
    ! Leaf C:N model
    coef_nimpl%Vcmax25_leafcn = getparreal( 'params/params_nimpl.dat', 'Vcmax25_leafcn' )
    coef_nimpl%LMA_leafcn = getparreal( 'params/params_nimpl.dat', 'LMA_leafcn' )
    coef_nimpl%intersect_leafcn     = getparreal( 'params/params_nimpl.dat', 'intersect_leafcn' )

    ! Constant
    coef_nimpl%root_cn = getparreal( 'params/params_nimpl.dat', 'root_cn' )
    coef_nimpl%wood_cn = getparreal( 'params/params_nimpl.dat', 'wood_cn' )

  end subroutine getpar_nimpl


  subroutine get_preds_nc( filnam, varname, domaininfo, grid, pred )
    !////////////////////////////////////////////////////////////////
    ! xxx add explanation here
    !----------------------------------------------------------------
    use netcdf

    ! arguments
    character(len=100) :: filnam
    character(len=100) :: varname
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
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

    ! Get longitude and latitude values
    call check( nf90_get_var( ncid, londimid, lon_arr ) )
    call check( nf90_get_var( ncid, latdimid, lat_arr ) )

    ! ! xxx try:
    ! lon_arr = (/ (i, i = 1,nlon_arr) /)
    ! lon_arr = (lon_arr - 1) * 0.5 - 180.0 + 0.25 
    ! lat_arr = (/ (i, i = 1,nlat_arr) /)
    ! lat_arr = (lat_arr - 1) * 0.5 - 90.0 + 0.25 

    ! Check if the resolution of the climate input files is identical to the model grid resolution
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of FPC input file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of FPC input file is not identical with model grid.'

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
    pred(:) = pred_arr(ilon(:),ilat(:))

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( pred_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    return

  end subroutine get_preds_nc


end module md_nuptake_impl
