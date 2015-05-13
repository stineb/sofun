  module _fluxes
    !////////////////////////////////////////////////////////////////
    ! Module contains daily model state variables. These variables have no 
    ! spatial dimension. At the end of the daily loop, values may be 
    ! copied to spatial arrays.
    ! All these variables are required by the model and are independent
    ! of which modules are used.
    ! xxx This is to be made more modular.
    !----------------------------------------------------------------
    use _params_core
    use _classdefs

    implicit none

    real, dimension(npft)           :: dgpp     ! gross primary production [gC/m2/d]
    type(carbon), dimension(nlu)    :: dnep     ! net ecosystem production [gC/m2/d]
    type(carbon), dimension(npft)   :: dnpp     ! net primary production [gC/m2/d]
    
    type(carbon), dimension(npft)   :: drsoil   ! soil respiration (only from exudates decomp.) [gC/m2/d]
    type(carbon), dimension(nlu)    :: drhet    ! heterotrophic respiration [gC/m2/d]
    real, dimension(npft)           :: drauto           ! autotrophic respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft)           :: drleaf           ! leaf maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft)           :: drroot           ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft)           :: drsapw           ! sapwood maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft)           :: drgrow           ! growth respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real, dimension(npft)           :: dcex             ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    
    real, dimension(nlu)            :: doc             ! daily soil turnover, used as surrogate for labile C (former 'doc') [gC/m2/d]
    
    type(nitrogen), dimension(npft) :: dnup   ! daily N uptake [gN/m2/d]
    
                                   ! variables without spatial dimension: have NO memory from one to next year
    type(orgpool), dimension(nlu)   :: aestab   ! annual C and N fixation due to establishment (=acflux_estab,an_uptake_estab)

  end module _fluxes

