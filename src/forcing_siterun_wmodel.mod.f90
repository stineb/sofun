module md_forcing_siterun
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
  use md_params_core, only: nmonth, ndaymonth, lunat, ndayyear, maxgrid, nlu, dummy
  use md_sofunutils, only: daily2monthly, read1year_daily, read1year_monthly, &
    getvalreal, monthly2daily_weather, monthly2daily

  implicit none

  private
  public getco2, getninput, ninput_type, gettot_ninput, getfapar, getclimate_site, getlanduse, landuse_type, climate_type

  type climate_type
    real, dimension(ndayyear) :: dtemp  ! deg C
    real, dimension(ndayyear) :: dprec  ! mm d-1
    real, dimension(ndayyear) :: dfsun  ! unitless
    real, dimension(ndayyear) :: dvpd   ! Pa
    real, dimension(ndayyear) :: dnetrad! W m-2
  end type climate_type

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

  function getco2( runname, sitename, forcingyear, const_co2_year, firstyeartrend, co2_forcing_file ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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
    pco2 = dummy

  end function getco2


  function getninput( ntype, runname, sitename, forcingyear, firstyeartrend, const_ninput_year, ninput_noy_forcing_file, ninput_nhx_forcing_file, climate ) result( out_getninput )
    !////////////////////////////////////////////////////////////////
    ! Dummy function
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: ntype   ! either "nfert" or "ndep"
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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
      out_gettot_ninput(jpngr)%dnoy(:) = ninput1(jpngr)%dnoy(:) + ninput2(jpngr)%dnoy(:)
      out_gettot_ninput(jpngr)%dnhx(:) = ninput1(jpngr)%dnhx(:) + ninput2(jpngr)%dnhx(:)
      out_gettot_ninput(jpngr)%dtot(:) = ninput1(jpngr)%dtot(:) + ninput2(jpngr)%dtot(:)
    end do

  end function gettot_ninput


  function getfapar( runname, sitename, forcingyear, fapar_forcing_source ) result( fapar_field )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in)          :: forcingyear
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    real, dimension(nmonth,maxgrid) :: fapar_field

    ! local variables 
    integer :: jpngr
    integer :: readyear
    character(len=4) :: faparyear_char

    fapar_field(:,jpngr) = dummy

    return
    999  format (I4.4)

  end function getfapar


  function getclimate_site( sitename, climateyear, in_netrad ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    ! Read year-2013 data after 2013
    !----------------------------------------------------------------    
    ! arguments
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_netrad

    ! local variables
    integer :: day
    integer :: jpngr = 1
    character(len=4) :: climateyear_char

    ! function return variable
    type( climate_type ), dimension(maxgrid) :: out_climate

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    jpngr = 1

    out_climate(jpngr)%dtemp(:)   = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dtemp_'//trim(sitename)//'_'//climateyear_char//'.txt')
    out_climate(jpngr)%dprec(:)   = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dprec_'//trim(sitename)//'_'//climateyear_char//'.txt')
    out_climate(jpngr)%dfsun(:)   = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dfsun_'//trim(sitename)//'_'//climateyear_char//'.txt')
    out_climate(jpngr)%dvpd(:)    = dummy
    if (in_netrad) then
      out_climate(jpngr)%dnetrad(:) = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dnetrad_'//trim(sitename)//'_'//climateyear_char//'.txt')
    else
      out_climate(jpngr)%dnetrad(:) = dummy
    end if

    return
    999  format (I4.4)

  end function getclimate_site


  function getlanduse( runname, sitename, forcingyear, do_grharvest_forcing_file, const_lu_year, firstyeartrend ) result( out_landuse )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's annual landuse state and harvesting regime (day of above-ground harvest)
    ! Grass harvest forcing file is read for specific year, if none is available,
    ! use earliest forcing file available. 
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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

end module md_forcing_siterun

