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
  use md_params_core, only: ntstepsyear

  implicit none

  private
  public climate_type, getclimate, getco2, forcingData, aggregate_climate

  ! type climate_type
  !   integer(kind=in), dimension(ntstepsyear) :: year          ! Year
  !   integer(kind=in), dimension(ntstepsyear) :: doy           ! day of the year
  !   real(kind=sp), dimension(ntstepsyear)    :: hod           ! hour of the day
  !   real(kind=sp), dimension(ntstepsyear)    :: PAR           ! umol m-2 s-1
  !   real(kind=sp), dimension(ntstepsyear)    :: radiation     ! W/m2
  !   real(kind=sp), dimension(ntstepsyear)    :: Tair          ! air temperature,  K
  !   real(kind=sp), dimension(ntstepsyear)    :: Tsoil         ! soil temperature, K
  !   real(kind=sp), dimension(ntstepsyear)    :: RH            ! relative humidity
  !   real(kind=sp), dimension(ntstepsyear)    :: rain          ! kgH2O m-2 s-1
  !   real(kind=sp), dimension(ntstepsyear)    :: windU         ! wind velocity (m s-1)
  !   real(kind=sp), dimension(ntstepsyear)    :: P_air         ! pa
  !   real(kind=sp), dimension(ntstepsyear)    :: CO2           ! ppm
  !   real(kind=sp), dimension(ntstepsyear)    :: soilwater     ! soil moisture, vol/vol
  ! end type climate_type

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
     real    :: CO2           ! ppm
     real    :: soilwater     ! soil moisture, vol/vol
  end type climate_type

  ! Input forcing data
  type(climate_type), pointer, save :: forcingData(:)


contains

  function aggregate_climate(forcing) result(forcing_agg)
    !////////////////////////////////////////////////////////////////
    ! Takes mean over fast time steps provided in input, and 
    ! returns aggregated values.
    !----------------------------------------------------------------
    type(climate_type), dimension(:), intent(in) :: forcing

    type(climate_type) :: forcing_agg
    integer :: nt

    nt = size(forcing)

    forcing_agg%year      = forcing(1)%year
    forcing_agg%doy       = forcing(1)%doy
    forcing_agg%hod       = 12.0
    forcing_agg%PAR       = sum( forcing(:)%PAR ) / nt         ! umol m-2 s-1
    forcing_agg%radiation = sum( forcing(:)%radiation ) / nt   ! W/m2
    forcing_agg%Tair      = sum( forcing(:)%Tair ) / nt        ! air temperature,  K
    forcing_agg%Tsoil     = sum( forcing(:)%Tsoil ) / nt       ! soil temperature, K
    forcing_agg%RH        = sum( forcing(:)%RH ) / nt          ! relative humidity
    forcing_agg%rain      = sum( forcing(:)%rain ) / nt        ! kgH2O m-2 s-1
    forcing_agg%windU     = sum( forcing(:)%windU ) / nt       ! wind velocity (m s-1)
    forcing_agg%P_air     = sum( forcing(:)%P_air ) / nt       ! pa
    forcing_agg%CO2       = sum( forcing(:)%CO2 ) / nt         ! ppm
    forcing_agg%soilwater = sum( forcing(:)%soilwater ) / nt   ! soil moisture, vol/vol

  end function aggregate_climate


  function getclimate( nt, forcing, climateyear_idx, climateyear ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! This function invokes file format specific "sub-functions/routines"
    ! to read from NetCDF. This nesting is necessary because this 
    ! cannot be done file-specific in SR sofun, but can be done here
    ! as this module is compilation-specific (only for global simulations)
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: nt ! number of time steps

    ! real(kind=dp),  dimension(nt,13), intent(in)  :: forcing  ! array containing all temporally varying forcing data (rows: time steps; columns: 1=air temperature, 2=rainfall, 3=vpd, 4=ppfd, 5=net radiation, 6=sunshine fraction, 7=snowfall, 8=co2, 9=N-deposition) 
    type(climate_type), dimension(nt), intent(in) :: forcing

    integer, intent(in) :: climateyear_idx, climateyear

    ! local variables
    integer :: idx_start, idx_end
    ! integer, dimension(2) :: shape_forcing
    real, parameter :: timestep = 1.0

    ! function return variable
    type(climate_type), dimension(ntstepsyear) :: out_climate

    idx_start = (climateyear_idx - 1) * ntstepsyear + 1
    idx_end   = idx_start + ntstepsyear - 1
    
    ! ! Test if forcing dimensions are correct
    ! shape_forcing = shape(forcing)
    ! if (idx_end>shape_forcing(1)) then
    !   stop 'forcing array size does not have enough rows.'
    ! end if

    ! if (int(forcing(idx_start, 1)) /= climateyear) stop 'getclimate(): climateyear does not correspond to index read from forcing'
    ! if (forcing(idx_start)%year /= climateyear) then
    !   print*,'forcing(idx_start)%year ', forcing(idx_start)%year
    !   print*,'climateyear ', climateyear
    !   print*,'climateyear_idx ', climateyear_idx
    !   stop 'getclimate(): climateyear does not correspond to index read from forcing'
    ! end if

    ! ! This is to read from ORNL file
    ! out_climate%year      = int(forcing(idx_start:idx_end, 1))           ! Year
    ! out_climate%doy       = int(forcing(idx_start:idx_end, 2))           ! day of the year
    ! out_climate%hod       = real(forcing(idx_start:idx_end, 3))           ! hour of the day
    ! out_climate%PAR       = real(forcing(idx_start:idx_end, 4)) * 2.0     ! umol/m2/s           ! umol m-2 s-1
    ! out_climate%radiation = real(forcing(idx_start:idx_end, 5))           ! W/m2
    ! out_climate%Tair      = real(forcing(idx_start:idx_end, 6)) + 273.16  ! air temperature, K
    ! out_climate%Tsoil     = real(forcing(idx_start:idx_end, 7)) + 273.16  ! soil temperature, K
    ! out_climate%RH        = real(forcing(idx_start:idx_end, 8)) * 0.01    ! relative humidity (0.xx)
    ! out_climate%rain      = real(forcing(idx_start:idx_end, 9))/(timestep * 3600) ! kgH2O m-2 s-1
    ! out_climate%windU     = real(forcing(idx_start:idx_end, 10))           ! wind velocity (m s-1)
    ! out_climate%P_air     = real(forcing(idx_start:idx_end, 11))           ! pa
    ! out_climate%soilwater = 0.8                                           ! soil moisture, vol/vol

    ! ! This is to read from ORNL file
    ! print*,'getclimate(): forcing'
    ! print*, forcing(idx_start)%year
    ! print*, forcing(idx_start)%doy
    ! print*, forcing(idx_start)%hod
    ! print*, forcing(idx_start)%PAR
    ! print*, forcing(idx_start)%radiation
    ! print*, forcing(idx_start)%Tair
    ! print*, forcing(idx_start)%Tsoil
    ! print*, forcing(idx_start)%rain
    ! print*, forcing(idx_start)%windU
    ! print*, forcing(idx_start)%P_air
    ! print*, forcing(idx_start)%RH
    ! print*, forcing(idx_start)%CO2
    ! print*, forcing(idx_start)%soilwater

    out_climate(:) = forcing(idx_start:idx_end)

    ! print*,'getclimate(): out_climate'
    ! print*, out_climate(1)%year
    ! print*, out_climate(1)%doy
    ! print*, out_climate(1)%hod
    ! print*, out_climate(1)%PAR
    ! print*, out_climate(1)%radiation
    ! print*, out_climate(1)%Tair
    ! print*, out_climate(1)%Tsoil
    ! print*, out_climate(1)%rain
    ! print*, out_climate(1)%windU
    ! print*, out_climate(1)%P_air
    ! print*, out_climate(1)%RH
    ! print*, out_climate(1)%CO2
    ! print*, out_climate(1)%soilwater

  end function getclimate


  function getco2( nt, forcing, forcingyear_idx, forcingyear ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    integer,  intent(in) :: nt ! number of time steps

    ! real(kind=dp),  dimension(nt,13), intent(in)  :: forcing  ! array containing all temporally varying forcing data (rows: time steps; columns: 1=air temperature, 2=rainfall, 3=vpd, 4=ppfd, 5=net radiation, 6=sunshine fraction, 7=snowfall, 8=co2, 9=N-deposition) 
    type(climate_type), dimension(nt), intent(in) :: forcing

    integer, intent(in) :: forcingyear_idx, forcingyear

    ! function return variable
    real, dimension(ntstepsyear) :: pco2

    ! local variables 
    integer :: idx_start, idx_end, year
    real, parameter :: timestep = 1.0

    idx_start = (forcingyear_idx - 1) * ntstepsyear + 1
    idx_end   = idx_start + ntstepsyear - 1

    ! if (int(forcing(idx_start, 1)) /= forcingyear) stop 'getco2(): forcingyear does not correspond to index read from forcing'
    if (forcing(idx_start)%year /= forcingyear) stop 'getco2(): forcingyear does not correspond to index read from forcing'

    ! pco2 = real(forcing(idx_start:idx_end, 12)) * 1.0e-6  ! mol/mol
    pco2 = forcing(idx_start:idx_end)%CO2

  end function getco2  


end module md_forcing

