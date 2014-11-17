      module statevars
      !////////////////////////////////////////////////////////////////
      !  Module contains all daily variables. These variables have no 
      !  spatial dimension. At the end of the daily loop, values may be 
      !  copied to spatial arrays.
      !----------------------------------------------------------------
        use params_core
        use classdefs

        implicit none
        real, dimension(npft) :: dtransp          ! daily transpiration [mm/d]

        type(carbon), dimension(nlu)  :: dnep     ! net ecosystem production [gC/m2/d]
        real, dimension(npft)          :: dgpp    ! gross primary production [gC/m2/d]
        type(carbon), dimension(npft) :: dnpp     ! net primary production [gC/m2/d]
         
        type(carbon), dimension(npft) :: drsoil   ! soil respiration (only from exudates decomp.) [gC/m2/d]
        type(carbon), dimension(nlu)  :: drhet    ! heterotrophic respiration [gC/m2/d]
        real, dimension(npft) :: drauto           ! autotrophic respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
        real, dimension(npft) :: drleaf           ! leaf maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
        real, dimension(npft) :: drroot           ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
        real, dimension(npft) :: drgrow           ! growth respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
        real, dimension(npft) :: dcex             ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
        real, dimension(npft) :: netmin_litt      ! net N mineralisation in litter decomposition
        real, dimension(nlu)  :: netmin_soil      ! net N mineralisation in soil decomposition

        real,dimension(nlu)   :: doc               ! daily soil turnover, used as surrogate for labile C (former 'doc') [gC/m2/d]

        type(nitrogen),dimension(npft) :: dnup    ! daily N uptake [gN/m2/d]
        real,dimension(npft) :: dnup_pas          ! daily passive N uptake [gN/m2/d]
        real,dimension(npft) :: dnup_act          ! daily active N uptake [gN/m2/d]  
        real,dimension(npft) :: dnup_fix          ! daily N uptake by plant symbiotic N fixation [gN/m2/d]
        real,dimension(npft) :: dnup_ret          ! daily N uptake [gN/m2/d]

        real,dimension(nlu)  :: dnfix_free        ! N fixation by free-living organisms [gN/m2/d]

        real,dimension(nlu)  :: dn2o              ! soil N2O emissions [gN/m2/d]
        real,dimension(nlu)  :: dn2               ! soil N2 emissions [gN/m2/d]
        real,dimension(nlu)  :: dno               ! soil NO emissions [gN/m2/d]
        real,dimension(nlu)  :: dnloss            ! total N loss (gaseous+leaching) [gN/m2/d]
        real,dimension(nlu)  :: ddenitr           ! gross denitrification [gN/m2/d]
        real,dimension(nlu)  :: dnitr             ! gross nitrification [gN/m2/d]
        real,dimension(nlu)  :: dvol              ! N volatilisation[gN/m2/d]
        real,dimension(nlu)  :: dleach            ! N leaching [gN/m2/d]


        !////////////////////////////////////////////////////////////////
        !  In LPX, these state variables were updated annually
        !----------------------------------------------------------------
        ! variables with 'jpngr' dimension: have memory from one to next year (veg. variables)
        ! therefore to be treated like pools.
        logical, dimension(npft,maxgrid)      :: ispresent      ! boolean whether PFT is present

        real, dimension(npft,maxgrid)         :: nind           ! number of individuals [1/m2]
        real, dimension(npft,maxgrid)         :: crownarea  
        real, dimension(npft,maxgrid)         :: lai_ind
        real, dimension(npft,maxgrid)         :: fpc_grid

        ! variables without spatial dimension: have NO memory from one to next year
        type(orgpool), dimension(nlu) :: aestab         ! annual C and N fixation due to establishment (=acflux_estab,an_uptake_estab)

      end module statevars

