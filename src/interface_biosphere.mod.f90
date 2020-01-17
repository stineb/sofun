module md_interface

  use, intrinsic :: iso_fortran_env, dp=>real64

  use md_forcing, only: climate_type  
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: paramstype_siml, outtype_steering
  use md_params_core, only: MSPECIES, ntstepsyear, ndayyear, MAX_INIT_COHORTS
  use md_grid, only: gridtype !, domaininfo_type

  implicit none

  private
  public interfacetype_biosphere, outtype_biosphere, myinterface  

  ! type paramstype_calib
  !   real :: kphio
  ! end type paramstype_calib  

  type paramstype_tile
    integer :: soiltype
    real :: FLDCAP
    real :: WILTPT
    real :: K1
    real :: K2
    real :: K_nitrogen
    real :: etaN
    real :: MLmixRatio
    real :: l_fract
    real :: retransN
    real :: f_initialBSW
  end type paramstype_tile
  
  type paramstype_species
    real, dimension(0:MSPECIES) :: lifeform
    real, dimension(0:MSPECIES) :: phenotype
    real, dimension(0:MSPECIES) :: pt
    real, dimension(0:MSPECIES) :: seedlingsize
    real, dimension(0:MSPECIES) :: LMA
    real, dimension(0:MSPECIES) :: phiRL
    real, dimension(0:MSPECIES) :: LNbase
    real, dimension(0:MSPECIES) :: laimax
    real, dimension(0:MSPECIES) :: LAI_light
    real, dimension(0:MSPECIES) :: Nfixrate0
    real, dimension(0:MSPECIES) :: NfixCost0
    real, dimension(0:MSPECIES) :: phiCSA
    real, dimension(0:MSPECIES) :: mortrate_d_c
    real, dimension(0:MSPECIES) :: mortrate_d_u
    real, dimension(0:MSPECIES) :: maturalage
    real, dimension(0:MSPECIES) :: fNSNmax
    real                        :: f_N_add
  end type paramstype_species

  type inittype_cohort 
    real, dimension(MAX_INIT_COHORTS) :: init_cohort_species
    real, dimension(MAX_INIT_COHORTS) :: init_cohort_nindivs
    real, dimension(MAX_INIT_COHORTS) :: init_cohort_bsw
    real, dimension(MAX_INIT_COHORTS) :: init_cohort_bHW
    real, dimension(MAX_INIT_COHORTS) :: init_cohort_nsc
  end type inittype_cohort

  type inittype_soil 
    real :: init_fast_soil_C
    real :: init_slow_soil_C
    real :: init_Nmineral
    real :: N_input
  end type inittype_soil

  type interfacetype_biosphere
    integer                                       :: year
    real, dimension(:), allocatable               :: pco2
    type(gridtype)                                :: grid
    type(climate_type), dimension(:), allocatable :: climate
    type(outtype_steering)                        :: steering
    type(paramstype_siml)                         :: params_siml
    real, dimension(:), allocatable               :: fpc_grid   ! allocatable because we don't know number of PFTs a priori
    ! type(paramstype_calib)                      :: params_calib    ! calibratable parameters
    type(paramstype_species)                      :: params_species
    type(paramtype_soil)                          :: params_soil
    type(paramstype_tile)                         :: params_tile
    type(inittype_cohort)                         :: init_cohort
    type(inittype_soil)                           :: init_soil
    integer                                       :: datalines
    integer                                       :: steps_per_day
    real                                          :: dt_fast_yr
    real                                          :: step_seconds
  end type interfacetype_biosphere

  type(interfacetype_biosphere) :: myinterface

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

contains


end module md_interface
