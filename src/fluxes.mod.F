  module fluxes
    !////////////////////////////////////////////////////////////////
    ! Module contains daily model state variables. These variables have no 
    ! spatial dimension. At the end of the daily loop, values may be 
    ! copied to spatial arrays.
    ! All these variables are required by the model and are independent
    ! of which modules are used.
    ! xxx This is to be made more modular.
    !----------------------------------------------------------------
    use params_core
    use classdefs

    implicit none

    type(carbon), dimension(nlu)  :: dnep     ! net ecosystem production [gC/m2/d]
    real, dimension(npft)         :: dgpp     ! gross primary production [gC/m2/d]
    type(carbon), dimension(npft) :: dnpp     ! net primary production [gC/m2/d]
     
    type(carbon), dimension(npft) :: drsoil   ! soil respiration (only from exudates decomp.) [gC/m2/d]
    type(carbon), dimension(nlu)  :: drhet    ! heterotrophic respiration [gC/m2/d]
    real, dimension(npft) :: drauto           ! autotrophic respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: drleaf           ! leaf maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: drroot           ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: drsapw           ! sapwood maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: drgrow           ! growth respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: dcex             ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft) :: netmin_litt      ! net N mineralisation in litter decomposition
    real, dimension(nlu)  :: netmin_soil      ! net N mineralisation in soil decomposition

    real, dimension(nlu)   :: doc             ! daily soil turnover, used as surrogate for labile C (former 'doc') [gC/m2/d]

    type(nitrogen), dimension(npft) :: dnup   ! daily N uptake [gN/m2/d]

    real, dimension(nlu)  :: dnfix_free       ! N fixation by free-living organisms [gN/m2/d]
    ! xxx add dnfix_symb

    real, dimension(nlu)  :: dn2o             ! soil N2O emissions [gN/m2/d]
    real, dimension(nlu)  :: dn2              ! soil N2 emissions [gN/m2/d]
    real, dimension(nlu)  :: dno              ! soil NO emissions [gN/m2/d]
    real, dimension(nlu)  :: dnloss           ! total N loss (gaseous+leaching) [gN/m2/d]
    real, dimension(nlu)  :: ddenitr          ! gross denitrification [gN/m2/d]
    real, dimension(nlu)  :: dnitr            ! gross nitrification [gN/m2/d]
    real, dimension(nlu)  :: dvol             ! N volatilisation[gN/m2/d]
    real, dimension(nlu)  :: dleach           ! N leaching [gN/m2/d]

    ! xxx this module is to be renamed, 'fluxes' is not appropriate
    real, dimension(nlu) :: dwn               ! ! soil water mass [mm = kg/m2]
    real, dimension(nlu) :: dro               ! daily runoff [mm]
    real, dimension(nlu) :: aet_d             ! daily AET (mm)
    real :: rn_d                              ! daytime net radiation [J/m2]


    ! variables without spatial dimension: have NO memory from one to next year
    type(orgpool), dimension(nlu) :: aestab   ! annual C and N fixation due to establishment (=acflux_estab,an_uptake_estab)

  end module fluxes

