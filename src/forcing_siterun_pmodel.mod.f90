module md_forcing
  !////////////////////////////////////////////////////////////////
  ! Module contains forcing variables (climate, co2, ...), and
  ! subroutines used to read forcing input files for a specific year
  ! ('forcingyear'), specifically for site scale simulations.
  ! This module is only used on the level of 'sofun', but not
  ! within 'biosphere', as these variables are passed on to 'biosphere'
  ! as arguments.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: nmonth, ndaymonth, lunat, ndayyear, maxgrid, nlu, dummy, npft
  use md_sofunutils, only: daily2monthly, read1year_daily, read1year_monthly, &
    getvalreal, monthly2daily_weather, monthly2daily, calc_patm

  use md_grid, only: gridtype, domaininfo_type  

  implicit none

  private
  public getco2, getfapar, getclimate, climate_type, get_fpc_grid, vegcover_type, &
    getninput, ninput_type, gettot_ninput, getlanduse, landuse_type

  type climate_type
    real :: dtemp  ! deg C
    real :: dprec  ! mm d-1
    real :: dsnow  ! mm d-1 water equivalents
    real :: dfsun  ! unitless
    real :: dvpd   ! Pa
    real :: dppfd  ! mol m-2 d-1
    real :: dpatm  ! Pa
    real :: dnetrad! W m-2
    real :: dvwind  ! m s-1
  end type climate_type

  type vegcover_type
    real :: dfapar ! fraction of absorbed photosynthetically active radiation
  end type vegcover_type

  type landuse_type
    real, dimension(nlu)         :: lu_area
    logical, dimension(ndayyear) :: do_grharvest
  end type landuse_type

  type ninput_type
    real, dimension(ndayyear) :: dnoy
    real, dimension(ndayyear) :: dnhx
    real, dimension(ndayyear) :: dtot
  end type ninput_type

contains

  function getclimate( domaininfo, grid, init, climateyear, in_ppfd, in_netrad ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! This function invokes file format specific "sub-functions/routines"
    ! to read from NetCDF. This nesting is necessary because this 
    ! cannot be done file-specific in SR sofun, but can be done here
    ! as this module is compilation-specific (only for global simulations)
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_ppfd
    logical, intent(in) :: in_netrad

    ! local variables
    character(len=4) :: climateyear_char
    character(len=256) :: snowf_name
    logical :: snowf_exists

    ! function return variable
    type( climate_type ), dimension(ndayyear,domaininfo%maxgrid) :: out_climate

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    out_climate(:,1)%dtemp = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dtemp_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    out_climate(:,1)%dprec = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dprec_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    out_climate(:,1)%dvpd  = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dvpd_' //trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    if (in_ppfd) then
      out_climate(:,1)%dppfd = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dppfd_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    else
      out_climate(:,1)%dppfd = dummy
    end if
    if (in_netrad) then
      out_climate(:,1)%dnetrad = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dnetrad_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    else
      out_climate(:,1)%dnetrad = dummy
    end if
    if ( in_netrad .and. in_ppfd ) then
      out_climate(:,1)%dfsun = dummy
    else
      out_climate(:,1)%dfsun = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dfsun_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    end if

    ! read snow file if it exists
    snowf_name = 'sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dsnow_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt'
    INQUIRE(FILE = trim(snowf_name), EXIST = snowf_exists)
    if (snowf_exists) then
      out_climate(:,1)%dsnow = read1year_daily(snowf_name)
    else
      out_climate(:,1)%dsnow = 0.0
    end if

    out_climate(:,1)%dpatm = calc_patm( grid(1)%elv)

    return
    999  format (I4.4)

  end function getclimate


  function getco2( runname, domaininfo, forcingyear, const_co2_year, firstyeartrend, co2_forcing_file ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in) :: forcingyear
    integer, intent(in) :: const_co2_year
    integer, intent(in) :: firstyeartrend
    character(len=*), intent(in) :: co2_forcing_file

    ! function return variable
    real :: pco2

    ! local variables 
    integer :: readyear

    if (const_co2_year/=int(dummy)) then
      readyear = const_co2_year
    else  
      readyear = forcingyear
    end if
    ! write(0,*) 'GETCO2: use CO2 data of year ', readyear
    pco2 = getvalreal( 'sitedata/co2/'//trim(domaininfo%domain_name)//'/'//trim(co2_forcing_file), readyear )

  end function getco2


  function getfapar( domaininfo, grid, year, fapar_forcing_source ) result( out_vegcover )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
    integer, intent(in) :: year
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    type( vegcover_type ), dimension(ndayyear,domaininfo%maxgrid) :: out_vegcover

    ! local variables 
    integer :: readyear
    character(len=4) :: faparyear_char

    if (trim(fapar_forcing_source)=='NA') then
      ! If in simulation parameter file 'NA' is specified for 'fapar_forcing_source'
      out_vegcover(:,1)%dfapar = dummy

    else
      ! Prescribed. Read monthly fAPAR value from file
      ! create 4-digit string for year  
      write(faparyear_char,999) min( max( 2000, year ), 2014 )
      out_vegcover(:,1)%dfapar = read1year_daily( 'sitedata/fapar/'//trim(domaininfo%domain_name)//'/'//faparyear_char//'/'//'dfapar_'//trim(domaininfo%domain_name)//'_'//faparyear_char//'.txt' )

      ! ! "Correct" fAPAR
      ! print*,"WARNING: normalising fAPAR to within 0.12 and 1.0."
      ! out_vegcover(:,1)%dfapar = max((out_vegcover(1)%dfapar(:) - 0.12), 0.0)/(1.0 - 0.12)

    end if

    return
    999  format (I4.4)

  end function getfapar


  function get_fpc_grid( domaininfo, grid, params_siml ) result( fpc_grid_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns the fractional land cover by vegetation types 
    ! based on the 10 IGBP types in the input file (MODIS Landcover)
    ! 1: ENF: type2 = "evergreen needleleaf forest" ;
    ! 2: EBF: type3 = "evergreen broadleaf forest" ;
    ! 3: DNF: type4 = "deciduous needleleaf forest" ;
    ! 4: DBF: type5 = "deciduous broadleaf forest" ;
    ! 5: MF:  type6 = "mixed forest" ;
    ! 6: SHR: type7+type8 = "closed shrublands" + "open shrublands";
    ! 7: SAV: type9+type10 = "savannas" plus "woody savannas"
    ! 8: GRA: type11 = "grasslands" ;
    ! 9: WET: type12 = "permanent wetlands" ;
    ! 10:CRO: type13 + type15 = "croplands" + "cropland (natural vegetation mosaic)";
    !----------------------------------------------------------------
    use md_params_siml, only: paramstype_siml
    use md_params_core, only: npft

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
    type( paramstype_siml ), intent(in) :: params_siml

    ! function return variable
    real, dimension(npft,domaininfo%maxgrid) :: fpc_grid_field

    ! local variables
    integer :: jpngr, pft

    jpngr = 1

    ! get binary information of PFT presence from simulation parameters
    fpc_grid_field(:,jpngr) = 0.0

    ! Code below must follow the same structure as in 'plant_pmodel.mod.f90'
    pft = 0
    if ( params_siml%lTrE ) then
      ! xxx dirty: call all non-grass vegetation types 'TrE', see indeces above
      pft = pft + 1
      fpc_grid_field(pft,jpngr) = 1.0
    end if 

    if ( params_siml%lGr3 ) then
      ! xxx dirty: call all grass vegetation types 'Gr3'
      pft = pft + 1
      fpc_grid_field(pft,jpngr) = 1.0
    end if

    if ( params_siml%lGr4 ) then
      ! xxx dirty: call all grass vegetation types 'Gr3'
      pft = pft + 1
      fpc_grid_field(pft,jpngr) = 1.0
    end if

    if (pft==0) stop 'get_fpc_grid: no PFT activated accoring to simulation parameter file.'

    if (pft/=npft) stop 'GET_FPC_GRID: Adjust npft manually in params_core.mod.f90'

    return
    999  format (I2.2)

  end function get_fpc_grid


  function getninput( ntype, runname, domaininfo, forcingyear, firstyeartrend, const_ninput_year, ninput_noy_forcing_file, ninput_nhx_forcing_file, climate ) result( out_getninput )
    !////////////////////////////////////////////////////////////////
    ! Dummy function
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: ntype   ! either "nfert" or "ndep"
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in)          :: forcingyear
    integer, intent(in) :: firstyeartrend
    integer, intent(in) :: const_ninput_year
    character(len=*), intent(in) :: ninput_noy_forcing_file
    character(len=*), intent(in) :: ninput_nhx_forcing_file
    type( climate_type ), dimension(maxgrid), intent(in) :: climate

    ! function return variable
    type( ninput_type ), dimension(maxgrid) :: out_getninput 

    ! local variables
    integer :: jpngr
    
    do jpngr=1,maxgrid
      out_getninput(jpngr)%dnoy(:) = dummy
      out_getninput(jpngr)%dnhx(:) = dummy
      out_getninput(jpngr)%dtot(:) = dummy
    end do

  end function getninput


  function gettot_ninput( ninput1, ninput2 ) result( out_gettot_ninput )
    !////////////////////////////////////////////////////////////////
    ! Function returns totals of two ninput type variables with 
    ! dimension maxgrid
    !----------------------------------------------------------------
    ! arguments
    type( ninput_type ), dimension(maxgrid), intent(in) :: ninput1, ninput2 

    ! local variables
    integer :: jpngr

    ! function return variable
    type( ninput_type ), dimension(maxgrid) :: out_gettot_ninput 

    do jpngr=1,maxgrid
      out_gettot_ninput(jpngr)%dnoy(:) = dummy
      out_gettot_ninput(jpngr)%dnhx(:) = dummy
      out_gettot_ninput(jpngr)%dtot(:) = dummy
    end do

  end function gettot_ninput


  function getlanduse( runname, domaininfo, forcingyear, do_grharvest_forcing_file, const_lu_year, firstyeartrend ) result( out_landuse )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's annual landuse state and harvesting regime (day of above-ground harvest)
    ! Grass harvest forcing file is read for specific year, if none is available,
    ! use earliest forcing file available. 
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in)          :: forcingyear
    character(len=*), intent(in), optional :: do_grharvest_forcing_file
    integer, intent(in) :: const_lu_year
    integer, intent(in) :: firstyeartrend

    ! local variables
    integer :: doy
    integer :: findyear
    real, dimension(ndayyear) :: tmp
    character(len=4) :: landuseyear_char
    character(len=245) :: filnam
    integer :: readyear
    logical :: file_exists

    ! function return variable
    type( landuse_type ) :: out_landuse

    ! xxx dummy
    out_landuse%lu_area(lunat)  = 1.0
    out_landuse%do_grharvest(:) = .false.

  end function getlanduse

end module md_forcing

