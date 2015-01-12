      module params_modl
      !////////////////////////////////////////////////////////////////
      !  Module contains (constrainable) model parameters.
      !  Model parameters adopted here are from LPX C3 grass PFT
      !  Litter and soil turnover parameters are divided by 365 to 
      !  convert from [1/yr] to [1/d].
      !----------------------------------------------------------------
#include "sofun_module_control.inc"
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
      real :: kbeer             ! 20 canopy light extinction coefficient
      real :: kresp             ! 21 additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")
      real :: MINIMUMCOSTFIX    ! 22 minimum cost of N fixation (at optimal temperature)
      real :: FIXOPTIMUM        ! 23 optimum temperature for N fixation
      real :: FIXWIDTH          ! 24 shape parameter for width of N fixation cost function
      real :: KN_ACTIVE_NO3     ! 25 N-availability constant in cost function for active NO3 uptake
      real :: KN_ACTIVE_NH4     ! 26 N-availability constant in cost function for active NH4 uptake
      real :: KC_ACTIVE_NO3     ! 27 root exploration constant in cost function for active NO3 uptake
      real :: KC_ACTIVE_NH4     ! 28 root exploration constant in cost function for active NH4 uptake
      real :: F_NRETAIN         ! 29 fraction of N retained at leaf abscission 
      real :: fpc_tree_max      ! 30 maximum fractional plant coverage of trees


      ! PFT-DEPENDENT PARAMETERS
      real, dimension(npft) :: cton_soil           ! C:N ratio of soil organic matter [1]
      integer, dimension(npft) :: lu_category      ! land use category associated with PFT
      logical, dimension(npft,nlu) :: islu         ! islu(ipft,ilu) is true if ipft belongs to ilu
      logical, dimension(npft) :: grass            ! boolean for growth form 'grass'
      logical, dimension(npft) :: tree             ! boolean for growth form 'tree'
      real, dimension(npft) :: lai0                ! sapling (or grass on initialisation) LAI (=lai_sapl)
      type(orgpool), dimension(npft) :: lm_sapl    ! sapling leaf mass
      type(orgpool), dimension(npft) :: sm_sapl    ! sapling leaf mass
      type(orgpool), dimension(npft) :: hm_sapl    ! sapling leaf mass
      type(orgpool), dimension(npft) :: rm_sapl    ! sapling root mass
      real, dimension(npft) :: long_leaf           ! leaf longeivity [days]
      real, dimension(npft) :: long_sapw           ! sapwood longeivity [days]
      real, dimension(npft) :: long_root           ! root longeivity [days]
      real, dimension(npft) :: sla                 ! specific leaf area [m2/gC]
      real, dimension(npft) :: cton_pro            ! C:N ratio of new production [1]
      integer, dimension(npft) :: pftcode          ! code for identification of PFT (1=C3 grass)
      real, dimension(npft) :: rresp_maint         ! maintenance respiration coefficient
      real, dimension(npft) :: rresp_growth        ! growth respiration coefficient
      logical, dimension(npft) :: nfixer           ! whether plant is capable of symbiotically fixing N
      real :: reinickerp                           ! Reinicker-p for geometry
      real, dimension(npft) :: wooddens            ! wood density (gC * m-3)
      real, dimension(npft) :: allom1, allom2      ! allometry parameters
      real :: allom3                               ! allometry parameter
      real :: crownarea_max                        ! maximum crown area
      real, dimension(npft) :: latosa              ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
      real, dimension(npft) :: woodtosapw          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
      real, dimension(npft) :: lmtorm0             ! leaf to root ratio under non-water stressed conditions


      !----------------------------------------------------------------
      ! Parameter values from Makela et al., 2008, Table 1
      !----------------------------------------------------------------
#if _gpp_makela
      ! Amount of foliage capturing 50% of maximum C gain [g/m2] (given as 2500 kg-leaf-DW/ha)
      real, parameter :: K_f = 500.0

      ! non-photosynthtic or structural N [gN/gC], assuming 0.5 gC/gDW (given as 0.009 kgN/kgDW)
      real, parameter :: n_struct = 0.018

      ! reference photosynthetic N [gN/gC], assuming 0.5 gC/gDW (given as 0.002 kgN/kgDW)
      real, parameter :: n_ref = 0.004

      ! N-saturated specific rate of photosynthesis [gC/g-foliage-C/yr] (given as 8.0 kgC/kgDW/yr)
      ! xxx try: This is per year. Adjust it to get daily value: 16.0/365.0
      real, parameter :: sigma_fm0 = 0.04383562

#endif

#if _allocation_makela
       ! 'Steady-state' pipe length coefficient [m/(gC/gN)], given as 2800 m/(kg-N/kg-DW), assuming 0.5 gC/g-DW
       real, parameter :: c_H = 5600.0

       ! Sapwood weight per unit foliage and pipe length [m-1]
       real, parameter :: alpha_w = 0.8

#endif

#if _nuptake_makela
      ! Amount of roots capturing 50% of available N [g/m2] (given as 2000 kg-leaf-DW/ha)
      real, parameter :: K_r = 400.0

#endif

      contains

      subroutine getpar_modl()
        !////////////////////////////////////////////////////////////////
        !  Subroutine reads model parameters from input file
        !----------------------------------------------------------------

        implicit none

        ! local variables
        integer :: pft
        integer, parameter :: npar = 29
        real, dimension(npar) :: params_array

        integer, parameter :: npar_pft = 21
        real, dimension(npar_pft) :: pft_params_array

        real, dimension(npft) :: lu_category_prov    ! land use category associated with PFT (provisional)
        real, dimension(3) :: out_get_geometry_sapl

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

        ! canopy light extinction coefficient
        kbeer = params_array(20)

        ! additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")
        kresp = params_array(21)

        ! shape parameter of cost function of N fixation 
        ! Below parameters (MINIMUMCOSTFIX, FIXOPTIMUM, FIXWIDTH ) are based on 
        ! the assumption that the cost of symbiotic N fixation is the 
        ! inverse of nitrogenase activity. 
        ! After Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
        ! (value from Gutschik 1981)
        MINIMUMCOSTFIX = params_array(22)

        ! shape parameter of cost function of N fixation 
        FIXOPTIMUM = params_array(23)

        ! shape parameter of cost function of N fixation 
        FIXWIDTH = params_array(24)

        ! N-availability constant in cost function for active NO3 uptake
        KN_ACTIVE_NO3 = params_array(25)
        
        ! N-availability constant in cost function for active NH4 uptake
        KN_ACTIVE_NH4 = params_array(26)

        ! root exploration constant in cost function for active NO3 uptake
        KC_ACTIVE_NO3 = params_array(27)
        
        ! root exploration constant in cost function for active NH4 uptake
        KC_ACTIVE_NH4 = params_array(28)

        ! fraction of N retained at leaf abscission 
        F_NRETAIN = params_array(29)        

        ! maximum fractional plant coverage of trees (sum of all tree PFTs)
        fpc_tree_max = params_array(30)

        !----------------------------------------------------------------
        ! PFT DEPENDENT PARAMETERS
        ! read parameter input file and store values in single array
        !----------------------------------------------------------------
        open(unit=02,file='params_pft.dat',status='OLD')      
        read (02,*) pft_params_array
        close (02)

        do pft=1,npft 

          ! C:N ratio of soil organic matter [1]
          cton_soil(pft) = pft_params_array(1) 

          ! C:N ratio of new production
          cton_pro(pft) = pft_params_array(6) 

          ! land use category associated with PFT (provisional)
          lu_category_prov(pft) = pft_params_array(2)
          if (lu_category_prov(pft)==1.0) then
            lu_category(pft) = lunat
            islu(pft,lunat) = .true.
          else
            islu(pft,lunat) = .false.
          end if

          ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
          if (pft_params_array(3)==3.0) then
            grass(pft) = .true.
            tree(pft)  = .false.
          else
            tree(pft) = .true.
            grass(pft) = .false.
          end if 

          ! sapling (or grass on initialisation) LAI
          lai0(pft) = pft_params_array(4) 

          ! leaf longevity [days], read in as [years]
          long_leaf(pft) = pft_params_array(5) * ndayyear 

          ! specific leaf area [m2/gC]
          ! defined as a function of 'long_leaf'
          ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
          ! Equation based on Reich et al 1997, Fig 1f:
          ! This is a mathematically identical expression as in Sitch et al., 2003, Eq. 6
          sla(pft) = 2.0e-4 * exp( 6.15 - 0.46 * log( long_leaf(pft) * 12.0 ) ) 

          ! PFT identification code (same as in LPX)
          ! 3 = Temperate Needle Evergreen   TeNE
          ! 8 = C3 grass                     TeH
          pftcode(pft) = int(pft_params_array(7)) 

          ! sapwood turnover period [days], read in as [years]
          long_sapw(pft) = pft_params_array(8) * ndayyear 

          ! root turnover period [days], read in as [years]
          long_root(pft) = pft_params_array(9) * ndayyear 

          ! maintenance respiration coefficient
          rresp_maint(pft) = pft_params_array(10) 

          ! N-fixing plant? (0=false, 1=true)
          if (pft_params_array(11)==0.) then
            nfixer(pft) = .false.
          else
            nfixer(pft) = .true.
          endif 

          ! growth respiration coefficient
          rresp_growth(pft) = pft_params_array(12) 

          ! Reinicker-p for geometry
          reinickerp = pft_params_array(13) ! put into pft-loop to have all geometry parameters in one place 

          ! wood density (gC * m-3)
          wooddens(pft) = pft_params_array(14) 

          ! allometry parameters
          allom1(pft) = pft_params_array(15)
          allom2(pft) = pft_params_array(16)
          allom3 = pft_params_array(17) ! put into pft-loop to have all geometry parameters in one place 

          ! maximum crown area
          crownarea_max = pft_params_array(18) ! put into pft-loop to have all geometry parameters in one place 

          ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
          latosa(pft) = pft_params_array(19)

          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
          woodtosapw(pft) = pft_params_array(20)

          ! leaf to root ratio under non-water stressed conditionss
          lmtorm0(pft) = pft_params_array(21)

          ! leaf (, sapwood, and heartwood) mass of sapling
          if (grass(pft)) then
            lm_sapl(pft)%c%c12 = lai0(pft) / sla(pft)
          else
            ! call function returning (lm_sapl, sm_sapl, hm_sapl)
            out_get_geometry_sapl = get_geometry_sapl( &
                lai0(pft), sla(pft), latosa(pft), reinickerp, allom1(pft) &
              , allom2(pft), allom3, wooddens(pft), woodtosapw(pft), lmtorm0(pft) &
              )            
            lm_sapl(pft)%c%c12 = out_get_geometry_sapl(1)
            sm_sapl(pft)%c%c12 = out_get_geometry_sapl(2)
            hm_sapl(pft)%c%c12 = out_get_geometry_sapl(3)
          end if
          
          ! root mass of sapling is a function of leaf mass            
          rm_sapl(pft)%c%c12 = ( 1.0 / lmtorm0(pft) ) * lm_sapl(pft)%c%c12

          ! N content of sapling
          lm_sapl(pft)%n%n14 = lm_sapl(pft)%c%c12 / cton_pro(pft)
          rm_sapl(pft)%n%n14 = rm_sapl(pft)%c%c12 / cton_pro(pft) 
          if (tree(pft)) then
            sm_sapl(pft)%n%n14 = sm_sapl(pft)%c%c12 / cton_pro(pft)
            hm_sapl(pft)%n%n14 = hm_sapl(pft)%c%c12 / cton_pro(pft)
          endif

        end do

        end subroutine getpar_modl


        function get_geometry_sapl( lai0, sla, latosa, reinickerp, allom1, &
          allom2, allom3, wooddens, woodtosapw, lmtorm0 )
          !////////////////////////////////////////////////////////////////////
          ! Function to determine sapling tree geometry. Returns an array 
          ! containing (lm_sapl, sm_sapl, hm_sapl), given geometry parameters.
          !--------------------------------------------------------------------
          use params_core, only: pi
          implicit none

          ! arguments
          real, intent(in) :: lai0                        ! sapling (or grass on initialisation) LAI (=lai_sapl)
          real, intent(in) :: sla                         ! specific leaf area [m2/gC]
          real, intent(in) :: latosa                      ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
          real, intent(in) :: reinickerp                  ! Reinicker-p for geometry
          real, intent(in) :: allom1, allom2, allom3  ! allometry parameters
          real, intent(in) :: wooddens                    ! wood density (gC * m-3)
          real, intent(in) :: woodtosapw                  ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
          real, intent(in) :: lmtorm0                     ! leaf to root ratio under non-water stressed conditionss

          ! local variables
          real :: x
          real :: lm_sapl
          real :: sm_sapl
          real :: hm_sapl
          real :: stemdiam
          real :: height_sapl

          ! function return variable
          real, dimension(3) :: get_geometry_sapl

          ! Calculate leafmass for a sapling individual
          !  (1) lai = leafmass * sla / (crown area)
          !  (2) (leaf area) = latosa * (sapwood xs area)
          !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
          !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
          !         (Reinickes theory)
          ! From (1),
          !  (4) leafmass = lai * (crown area) / sla
          ! From (1) & (3),
          !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
          ! From (2),
          !  (6) leafmass = latosa * (sapwood xs area) / sla
          !  (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
          ! From (6) and (7),
          !  (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
          ! From (8),
          !  (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
          ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
          ! Define x,
          ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] / 
          !          (sapwood diameter)
          ! From (10) & (11),
          ! (12) (stem diameter) = x * (sapwood diameter)
          ! From (5), (9) & (12),
          ! (13) leafmass = lai * allom1 * x**reinickerp * 
          !               (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
          ! From (13),
          ! (14) leafmass = [ lai * allom1 * x**reinickerp *
          !      (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))
          x = woodtosapw
           
          lm_sapl = (lai0*allom1*x**reinickerp*  &
            (4.0*sla/pi/latosa)**(reinickerp*0.5) / &
            sla)**(2.0/(2.0-reinickerp)) !eqn 14

          ! Calculate sapling stem diameter
          ! From (9) & (12),
          ! (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5
          stemdiam = x*(4.0*lm_sapl*sla/pi/latosa)**0.5 !Eqn 15
        
          !     Calculate sapling height
          ! (16) height = allom2 * (stem diameter)**allom3 (source?)
          height_sapl = allom2*stemdiam**allom3 !Eqn 16
        
          ! Calculate sapling sapwood mass
          ! (17) (sapwood volume) = height * (sapwood xs area)
          ! (18) (sapwood xs area) = leafmass * sla / latosa
          ! From (17) & (18),
          ! (19) (sapwood volume) = height * leafmass * sla / latosa
          ! (20) (sapwood mass) = (wood density) * (sapwood volume)
          ! From (19) & (20),
          ! (21) (sapwood mass) = (wood density) * height * leafmass * sla /
          ! latosa
          sm_sapl = wooddens*height_sapl*lm_sapl*sla/latosa !Eqn 21

          ! Calculate sapling heartwood mass
          !     From (11),
          ! (22) (heartwood mass) = (x-1) * (sapwood mass)
          hm_sapl = (x-1.0)*sm_sapl !Eqn 22

          get_geometry_sapl = (/lm_sapl, sm_sapl, hm_sapl/)

        end function get_geometry_sapl


      end module params_modl
