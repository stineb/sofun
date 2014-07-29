       module params_modl
       !////////////////////////////////////////////////////////////////
       !  Module contains (constrainable) model parameters.
       !  Model parameters adopted here are from LPX C3 grass PFT
       !  Litter and soil turnover parameters are divided by 365 to 
       !  convert from [1/yr] to [1/d].
       !----------------------------------------------------------------
       use params_core
       use classdefs

       implicit none

       ! NON PFT-DEPENDENT PARAMETERS
       real :: klitt_af10        !  1 above-ground fast (leaf) litter decay rate [1/d]
       real :: klitt_as10        !  2 above-ground slow (woody) litter decay rate [1/d] 
       real :: klitt_bg10        !  3 below-ground (root) litter decay rate [1/d] 
       real :: kexu10            !  4 exudates decay rate [1/d]
       real :: ksoil_fs10        !  5 fast soil pool decay rate [1/d]
       real :: ksoil_sl10        !  6 slow soil pool decay rate [1/d]
       real :: ntoc_crit1        !  7 factor for "Manzoni Equation" (XPXXX) [1]
       real :: ntoc_crit2        !  8 exponent for "Manzoni Equation" (XPXXX) [1]
       real :: cton_microb       !  9 C:N ratio of microbial biomass [1]
       real :: fastfrac          ! 10 fraction of litter input to fast soil pool [1]
       real :: rcton_LS          ! 11 constant ratio between leaf and sapwood C/N ratios (Friend et al.,2008 in XP08)
       real :: rcton_LR          ! 12 constant ratio between leaf and root C/N ratios (Friend et al., 2008 in XP08)
       real :: MAXNITR           ! 13 maximum nitrification rate
       real :: NON               ! 14 maximum NO from nitrification (day-1)
       real :: N2ON              ! 15 maximum N2O from nitrification (day-1)
       real :: KN                ! 16 Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
       real :: KDOC              ! 17 Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
       real :: DOCMAX            ! 18 DOCMAX
       real :: DNITR2N2O         ! 19 Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
       real :: lmtorm0           ! 20 leaf to root ratio under non-water stressed conditionss
       real :: kbeer             ! 21 canopy light extinction coefficient
       real :: kresp             ! 22 additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")


       ! PFT-DEPENDENT PARAMETERS
       real, dimension(npft) :: cton_soil           ! 1 C:N ratio of soil organic matter [1]
       real, dimension(npft) :: lu_category         ! 2 land use category associated with PFT
       logical, dimension(npft) :: grass            ! boolean for growth form 'grass'
       logical, dimension(npft) :: tree             ! boolean for growth form 'tree'
       real, dimension(npft) :: lai0                ! sapling (or grass on initialisation) LAI (=lai_sapl)
       type(orgpool), dimension(npft) :: lm_sapl    ! sapling leaf mass
       type(orgpool), dimension(npft) :: rm_sapl    ! sapling root mass
       real, dimension(npft) :: long_leaf           ! leaf longeivity [yr]
       real, dimension(npft) :: long_sapw           ! sapwood longeivity [yr]
       real, dimension(npft) :: long_root           ! root longeivity [yr]
       real, dimension(npft) :: sla                 ! specific leaf area [m2/gC]
       real, dimension(npft) :: cton_pro            ! C:N ratio of new production [1]
       integer, dimension(npft) :: pftcode          ! code for identification of PFT (1=C3 grass)
       real, dimension(npft) :: respcoeff           ! maintenance respiration coefficient

       contains

       subroutine getpar_modl
       !////////////////////////////////////////////////////////////////
       !  Subroutine reads model parameters from input file
       !----------------------------------------------------------------

       implicit none

       ! local variables
       integer :: pft
       integer, parameter :: npar = 22
       real, dimension(npar) :: params_array

       integer, parameter :: npar_pft = 10
       real, dimension(npar_pft) :: pft_params_array

       real, dimension(npft) :: lu_category_prov    ! land use category associated with PFT (provisional)


       !----------------------------------------------------------------
       ! NON-PFT DEPENDENT PARAMETERS
       ! read parameter input file and store values in single array
       !----------------------------------------------------------------
       open(unit=01,file='params_modl.dat',status='OLD')      
       read (01,*) params_array
       close (01)

       ! above-ground fast (leaf) litter decay rate [1/d]
       klitt_af10 = params_array(1)/ndayyear

       ! above-ground slow (woody) litter decay rate [1/d] 
       klitt_as10 = params_array(2)/ndayyear

       ! below-ground (root) litter decay rate [1/d] 
       klitt_bg10 = params_array(3)/ndayyear

       ! exudates decay rate [1/d]
       kexu10 = params_array(4)/ndayyear

       ! fast soil pool decay rate [1/d]
       ksoil_fs10 = params_array(5)/ndayyear

       ! slow soil pool decay rate [1/d]
       ksoil_sl10 = params_array(6)/ndayyear

       ! factor for "Manzoni Equation" (XPXXX) [1]
       ntoc_crit1 = params_array(7)

       ! exponent for "Manzoni Equation" (XPXXX) [1]
       ntoc_crit2 = params_array(8)

       ! C:N ratio of microbial biomass [1]
       cton_microb = params_array(9)

       ! fraction of litter input to fast soil pool [1]
       fastfrac = params_array(10)

       ! constant ratio between leaf and sapwood C/N ratios (Friend et al.,2008 in XP08)
       rcton_LS = params_array(11)

       ! constant ratio between leaf and root C/N ratios (Friend et al., 2008 in XP08)
       rcton_LR = params_array(12)

       ! maximum nitrification rate
       MAXNITR = params_array(13)

       ! maximum NO from nitrification (day-1)
       NON = params_array(14)

       ! maximum N2O from nitrification (day-1)
       N2ON = params_array(15)

       ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
       KN = params_array(16)

       ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
       KDOC = params_array(17)

       ! DOCMAX
       DOCMAX = params_array(18)

       ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
       DNITR2N2O = params_array(19)

       ! leaf to root ratio under non-water stressed conditions
       lmtorm0 = params_array(20)

       ! canopy light extinction coefficient
       kbeer = params_array(21)

       ! additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")
       kresp = params_array(22)


       !----------------------------------------------------------------
       ! PFT DEPENDENT PARAMETERS
       ! read parameter input file and store values in single array
       !----------------------------------------------------------------
       open(unit=02,file='params_pft.dat',status='OLD')      
       read (02,*) pft_params_array
       close (02)

       ! initialise
       grass(:) = .false.
       tree(:)  = .false.

       do pft=1,npft

              ! C:N ratio of soil organic matter [1]
              cton_soil(pft) = pft_params_array(1)

              ! C:N ratio of new production
              cton_pro(pft) = pft_params_array(6)

              ! land use category associated with PFT (provisional)
              lu_category_prov(pft) = pft_params_array(2)
              if (lu_category_prov(pft)==1.) lu_category(pft) = lunat

              ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
              if (pft_params_array(3)==3.) then
                     grass(pft) = .true.
              else
                     tree(pft) = .true.
              end if

              ! sapling (or grass on initialisation) LAI
              lai0(pft) = pft_params_array(4)

              ! leaf mass of sapling
              if (grass(pft)) then
                     lm_sapl(pft)%c%c12 = lai0(pft)/sla(pft)
                     lm_sapl(pft)%n%n14 = lm_sapl(pft)%c%c12/cton_pro(pft)
              end if
              rm_sapl(pft)%c%c12 = (1.0/lmtorm0)*lm_sapl(pft)%c%c12
              rm_sapl(pft)%n%n14 = rm_sapl(pft)%c%c12/cton_pro(pft)

              ! leaf longevity [yr]
              long_leaf(pft) = pft_params_array(5)

              ! specific leaf area [m2/gC]
              ! defined as a function of 'long_leaf'
              ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
              ! Equation based on Reich et al 1997, Fig 1f:
              sla(pft)=2.0d-4*exp(6.15d0-0.46d0*log(long_leaf(pft)*12.0d0))

              ! PFT identification code (1=C3 grass)
              pftcode(pft) = int(pft_params_array(7))

              ! sapwood turnover period
              long_sapw(pft) = pft_params_array(8)

              ! root turnover period
              long_root(pft) = pft_params_array(9)

              ! maintenance respiration coefficient
              respcoeff(pft) = pft_params_array(10)

       end do

       end subroutine getpar_modl


       end module params_modl
