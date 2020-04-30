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
  use, intrinsic :: iso_fortran_env, dp=>real64, sp=>real32, in=>int32
  use md_params_core, only: ntstepsyear, ndayyear, kTkelvin

  implicit none

  private
  public climate_type, getclimate, getco2, forcingData

  type :: climate_type
    integer :: year          ! Year
    integer :: doy           ! day of the year
    real    :: hod           ! hour of the day
    real    :: PAR           ! umol m-2 s-1
    real    :: radiation     ! W/m2
    real    :: Tair          ! air temperature,  K
    real    :: Tsoil         ! soil temperature, K
    real    :: RH            ! relative humidity
    real    :: rain          ! kgH2O m-2 s-1
    real    :: windU         ! wind velocity (m s-1)
    real    :: P_air         ! pa
    real    :: CO2           ! mol CO2/mol dry air
    real    :: soilwater     ! soil moisture, vol/vol

    ! new:
    real    :: vpd           ! vapour pressure deficit (Pa)

  end type climate_type

  ! Input forcing data
  type(climate_type), pointer, save :: forcingData(:)


contains

  function getclimate( nt, ntstepsyear, ntstepsyear_forcing, daily, forcing, climateyear_idx, climateyear, elv ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! This function invokes file format specific "sub-functions/routines"
    ! to read from NetCDF. This nesting is necessary because this 
    ! cannot be done file-specific in SR sofun, but can be done here
    ! as this module is compilation-specific (only for global simulations)
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: nt ! number of total time steps in forcing array 
    integer, intent(in) :: ntstepsyear   ! number of time steps per year of model
    integer, intent(in) :: ntstepsyear_forcing  ! number of time steps per year of forcing data
    logical, intent(in) :: daily
    type(climate_type), dimension(nt), intent(in) :: forcing
    integer, intent(in) :: climateyear_idx, climateyear
    real, intent(in)    :: elv

    ! local variables
    integer :: idx_start, idx_end, it

    ! function return variable
    type(climate_type), dimension(ntstepsyear) :: out_climate

    idx_start = (climateyear_idx - 1) * ntstepsyear_forcing + 1
    idx_end   = idx_start + ntstepsyear_forcing - 1

    out_climate(:) = forcing(idx_start:idx_end)      

    ! get additional variables
    do it=1,ntstepsyear
      out_climate(it)%vpd  = calc_vpd_rh( out_climate(it)%RH, (out_climate(it)%Tair - kTkelvin) )
    end do

    ! aggregate to daily
    if (daily) then
      out_climate(:) = aggregate_climate_byday( forcing(idx_start:idx_end) )
    end if

  end function getclimate


  function aggregate_climate_byday(forcing) result(forcing_agg)
    !////////////////////////////////////////////////////////////////
    ! Takes mean over fast time steps provided in input, and 
    ! returns aggregated values.
    !----------------------------------------------------------------
    type(climate_type), dimension(:), intent(in) :: forcing

    ! function return variable
    type(climate_type), dimension(ndayyear) :: forcing_agg

    ! local
    integer :: nt  ! number of time steps per year (may vary)
    integer :: doy, idx_start, idx_end
    real :: nt_day

    nt = size(forcing, 1)
    nt_day = nt / ndayyear

    do doy = 1, ndayyear

      idx_start = (doy - 1) * nt_day + 1
      idx_end = idx_start + nt_day - 1

      forcing_agg(doy)%year      = forcing(idx_start)%year
      forcing_agg(doy)%doy       = forcing(idx_start)%doy
      forcing_agg(doy)%hod       = 12.0
      forcing_agg(doy)%PAR       = sum( forcing(idx_start:idx_end)%PAR ) / nt_day         ! umol m-2 s-1
      forcing_agg(doy)%radiation = sum( forcing(idx_start:idx_end)%radiation ) / nt_day   ! W/m2
      forcing_agg(doy)%Tair      = sum( forcing(idx_start:idx_end)%Tair ) / nt_day        ! air temperature,  K
      forcing_agg(doy)%Tsoil     = sum( forcing(idx_start:idx_end)%Tsoil ) / nt_day       ! soil temperature, K
      forcing_agg(doy)%RH        = sum( forcing(idx_start:idx_end)%RH ) / nt_day          ! relative humidity
      forcing_agg(doy)%rain      = sum( forcing(idx_start:idx_end)%rain ) / nt_day        ! kgH2O m-2 s-1
      forcing_agg(doy)%windU     = sum( forcing(idx_start:idx_end)%windU ) / nt_day       ! wind velocity (m s-1)
      forcing_agg(doy)%P_air     = sum( forcing(idx_start:idx_end)%P_air ) / nt_day       ! pa
      forcing_agg(doy)%CO2       = sum( forcing(idx_start:idx_end)%CO2 ) / nt_day         ! ppm
      forcing_agg(doy)%soilwater = sum( forcing(idx_start:idx_end)%soilwater ) / nt_day   ! soil moisture, vol/vol

    end do

  end function aggregate_climate_byday


  function getco2( nt, ntstepsyear, ntstepsyear_forcing, daily, forcing, forcingyear_idx, forcingyear ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    integer,  intent(in) :: nt ! number of time steps
    integer, intent(in) :: ntstepsyear   ! number of time steps per year of model
    integer, intent(in) :: ntstepsyear_forcing  ! number of time steps per year of forcing data
    logical, intent(in) :: daily
    type(climate_type), dimension(nt), intent(in) :: forcing
    integer, intent(in) :: forcingyear_idx, forcingyear
    ! function return variable
    real, dimension(ntstepsyear) :: pco2
    ! local variables 
    integer :: idx_start, idx_end, year
    idx_start = (forcingyear_idx - 1) * ntstepsyear_forcing + 1
    idx_end   = idx_start + ntstepsyear_forcing - 1
    if (daily) then
      pco2(:) = aggregate_co2_byday( forcing(idx_start:idx_end)%CO2 )
    else
      pco2(:) = forcing(idx_start:idx_end)%CO2  
    end if
    if (forcing(idx_start)%year /= forcingyear) stop 'getco2(): forcingyear does not correspond to index read from forcing'
  
  end function getco2


  function aggregate_co2_byday(vec) result(vec_agg)
    !////////////////////////////////////////////////////////////////
    ! Takes mean over fast time steps provided in input, and 
    ! returns aggregated values.
    !----------------------------------------------------------------
    real, dimension(:), intent(in) :: vec
    ! function return variable
    real, dimension(ndayyear) :: vec_agg
    ! local
    integer :: nt  ! number of time steps per year (may vary)
    integer :: doy, idx_start, idx_end
    real :: nt_day

    nt = size(vec)
    nt_day = nt / ndayyear

    do doy = 1, ndayyear
      idx_start = (doy - 1) * nt_day + 1
      idx_end = idx_start + nt_day - 1
      vec_agg(doy) = sum(vec(idx_start:idx_end)) / nt_day
    end do
  
  end function aggregate_co2_byday


  function calc_vpd_rh( rh, tc ) result( vpd )
    !////////////////////////////////////////////////////////////////////////
    ! Calculates vapor pressure deficit
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: rh      ! relative humidity (fraction, <1)
    real, intent(in)    :: tc      ! daily mean air temperature (deg C), daily varying from WATCH-WFDEI (ACTUALLY NOT USED)

    ! function return variable
    real :: vpd         ! vapor pressure deficit (Pa)

    ! local variables
    real :: esat       ! saturation water vapor pressure (Pa) at given air temperature

    esat = 611.0 * exp( (17.27 * tc)/(tc + 237.3) )

    vpd = esat * (1.0 - rh)

  end function calc_vpd_rh

end module md_forcing