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
  ! Variables without memory (not necessarily just fluxes)
  !----------------------------------------------------------------
  type canopy_nimpl_fluxes_type
    real :: anpp        ! annual total net primary production (gC m-2 yr-1)
    real :: aanpp       ! annual aboveground net primary production (gC m-2 yr-1)
    real :: abnpp       ! annual belowground net primary production (gC m-2 yr-1)
    real :: alnpp       ! annual net primary production for leaf production (gC m-2 yr-1)
    real :: awnpp       ! annual net primary production for wood production (gC m-2 yr-1)
  end type canopy_nimpl_fluxes_type

  type tile_nimpl_fluxes_type
    type(canopy_nimpl_fluxes_type) :: canopy
  end type tile_nimpl_fluxes_type

  !----------------------------------------------------------------
  ! object containing module-specific fluxes
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
    xxx

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
    xxx

  end type preds_nimpl_type

  type(preds_nimpl_type), dimension(:), allocatable :: preds_nimpl

  !----------------------------------------------------------------
  ! Specify file and variable names for NetCDF reading
  !----------------------------------------------------------------
  character(len=100), parameter :: filnam_cnsoil = "actual_file_name.nc"
  character(len=100), parameter :: varnam_cnsoil = "actual_variable_name"
  xxx

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
      xxx


      !--------------------------------------------------------------
      ! Read coefficients, populates 'coef_nimpl'
      !--------------------------------------------------------------
      ! BP:GPP model
      call getpar_nimpl()

    end if

    !--------------------------------------------------------------
    ! Predict using statistical models
    !--------------------------------------------------------------
    ! XXX loop or make tile_nimpl_fluxes a field?
    tile_nimpl_fluxes(lu)%canopy%anpp = coef_nimpl%cnsoil_bp * preds_nimpl(:)%cnsoil + coef_nimpl%age_bp * preds_nimpl(:)%age + xxx + coef_nimpl%intersect_bp

    ! xxx todo: other linear models for prediction

  end subroutine nuptake_impl


  subroutine getpar_nimpl()
    !////////////////////////////////////////////////////////////////
    ! xxx here comes a description
    !--------------------------------------------------------------
    use md_sofunutils, only: getparreal

    coef_nimpl%cnsoil_bp  = getparreal( 'params/params_nimpl.dat', 'cnsoil_bp' )
    coef_nimpl%age_bp     = getparreal( 'params/params_nimpl.dat', 'age_bp' )

    xxx read other coefficients and C:N ratios in different tissues from the same file

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
