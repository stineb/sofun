program main
  !////////////////////////////////////////////////////////////////
  ! Demo program for calculating GPP using the Fortran 90 
  ! implementation of the P-model.
  !
  ! Compile by:
  !   make demo_pmodel
  !
  ! Run by:
  !   echo <temp> <vpd> <co2> <ppfd> <fapar> <elv> | ./rundemo_pmodel
  !
  ! Author: Benjamin D. Stocker
  !----------------------------------------------------------------
  use md_gpp, only: getlue, calc_dgpp
  use md_plant, only: plant_type, params_pft_plant_type, init_plant, init_params_pft_plant

  implicit none

  !----------------------------------------------------------------
  ! Local variables
  !----------------------------------------------------------------
  real :: temp    ! relevant temperature (deg C)
  real :: vpd     ! vapour pressure deficit (Pa)
  real :: co2     ! ambient CO2 (ppm)
  real :: ppfd    ! daily total photosynthetic photon flux (mol m-2 d-1)
  real :: fapar   ! fraction of absorbed photosynthetically active radiation (unitless) 
  real :: elv     ! elevation above sea level (m)

  real :: lue     ! light use efficiency (mol CO2 (mol photon)-1)
  real :: gpp     ! gross primary productivity (g C m-2 d-1)

  type(plant_type) :: plant
  type(params_pft_plant_type) :: params_pft_plant

  !----------------------------------------------------------------
  ! Read arguments from standard input.
  ! These have to be specified in a fixed order.
  !----------------------------------------------------------------
  read (*,*) temp vpd co2 ppfd fapar elv

  !----------------------------------------------------------------
  ! Initialise stuff (needed unfortunately to make it compatible 
  ! with SOFUN structure)
  !----------------------------------------------------------------
  plant = init_plant()
  params_pft_plant = init_params_pft_plant()

  !----------------------------------------------------------------
  ! Calculate light use efficiency (LUE).
  !----------------------------------------------------------------
  lue = getlue( co2, temp, vpd, elv )

  !----------------------------------------------------------------
  ! Calculate GPP, given LUE
  !----------------------------------------------------------------
  gpp = calc_dgpp( fapar, 1.0, ppfd, lue, tempstress, soilmstress )

  !----------------------------------------------------------------
  ! Write gpp to standard output
  !----------------------------------------------------------------
  write(0,*) gpp

end program main
