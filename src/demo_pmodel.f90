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
  use md_interface, only: get_interface
  use md_plant, only: plant_type, plant_fluxes_type, initglobal_plant, params_pft_plant, getpar_modl_plant, initdaily_plant
  use md_gpp, only: pmodel, calc_dgpp, params_pft_gpp, outtype_pmodel, getpar_modl_gpp

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

  real :: gpp     ! gross primary productivity (g C m-2 d-1)
  real :: tempstress = 1.0
  real :: soilmstress = 1.0 

  integer :: pft = 1  ! index of plant functional type dimension

  type(plant_type),        allocatable, dimension(:,:) :: plant
  type(plant_fluxes_type), allocatable, dimension(:)   :: plant_fluxes
  type(outtype_pmodel) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

  !----------------------------------------------------------------
  ! Read arguments from standard input.
  ! These have to be specified in a fixed order.
  !----------------------------------------------------------------
  read (*,*) temp, vpd, co2, ppfd, fapar, elv

  !----------------------------------------------------------------
  ! Initialise the interface derived type. This is unfortunately unavoidable because 'interface' is referred to inside the md_gpp module
  !----------------------------------------------------------------
  call get_interface()

  !----------------------------------------------------------------
  ! Initialise stuff (needed unfortunately to make it compatible 
  ! with SOFUN structure)
  !----------------------------------------------------------------
  call getpar_modl_plant()
  call getpar_modl_gpp()
      
  allocate( plant( 1, 1 ) )
  allocate( plant_fluxes( 1 ) )

  call initglobal_plant( plant(:,:), 1 )
  call initdaily_plant( plant_fluxes(:) )

  !----------------------------------------------------------------
  ! Calculate GPP
  !----------------------------------------------------------------
  pft = 1
  out_pmodel = pmodel( params_pft_gpp(pft)%kphio, fapar = fapar, ppfd = ppfd, co2 = co2, tc = temp, vpd = vpd, elv = elv, c4 = params_pft_plant(pft)%c4, method_optci = "prentice14", method_jmaxlim = "wang17" )

  ! !----------------------------------------------------------------
  ! ! Calculate GPP, given LUE
  ! !----------------------------------------------------------------
  ! gpp = calc_dgpp( fapar = fapar, fpc_grid = 1.0, dppfd = ppfd, lue = out_pmodel%lue, tempstress = tempstress, soilmstress = soilmstress )

  !----------------------------------------------------------------
  ! Write gpp to standard output
  !----------------------------------------------------------------
  print*, out_pmodel
  ! write(0,*) gpp

end program main
