      subroutine littersom( jpngr )
      !////////////////////////////////////////////////////////////////
      !  Litter and SOM decomposition and nitrogen mineralisation.
      !  1st order decay of litter and SOM pools, governed by temperature
      !  and soil moisture following LPJ (Sitch et al., 2003) and 
      !  Xu-Ri & Prentice (XXX).
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_core
      use classdefs
      use rates
      use dailyvars, only: drsoil, drh, dnfix_free, doc
      use params_modl
      use pools, only: pexud, plitt_af, plitt_as, plitt_bg, psoil_fs, psoil_sl, nh4

      implicit none

      ! ARGUMENTS
      integer, intent(in) :: jpngr                    ! grid cell number

      ! LOCAL VARIABLES
      integer :: lu                                   ! counter variable for landuse class
      integer :: pft                                  ! counter variable for PFT number

      ! decay constants
      real, dimension(npft) :: klitt_af               ! decay rate, above-ground fast (leaf) litter (= k_litter_leaf)
      real, dimension(npft) :: klitt_as               ! decay rate, above-ground slow (woody) litter (= k_litter_woody)
      real, dimension(npft) :: klitt_bg               ! decay rate, below-ground fast litter (= k_litter_root)
      real, dimension(npft) :: kexu                   ! decay rate, exudates (= k_exu)
      real, dimension(nlu)  :: ksoil_fs               ! decay rate, fast soil (= k_fast)
      real, dimension(nlu)  :: ksoil_sl               ! decay rate, slow soil (= k_slow)

      ! temporary pools
      type (carbon)  :: dexu                          ! exudates decomposed in time step (= exu_decom)
      type (orgpool) :: dlitt_af                      ! above-ground fast litter decomposed per time step (= litterdag_fast)
      type (orgpool) :: dlitt_as                      ! above-ground slow litter decomposed per time step (= litterdag_slow)
      type (orgpool) :: dlitt_bg                      ! below-ground slow litter decomposed per time step (= litter_decom_bg)
      type (orgpool), dimension(npft) :: dlitt        ! total litter decomposed per time step (= litter_decom_pft)

      type (orgpool) :: dsoil_fs                      ! (= cflux_fast_atmos)
      type (orgpool) :: dsoil_sl                      ! (= cflux_slow_atmos)

     ! temporary variables
      real :: eff                                     ! microbial growth efficiency 
      real :: ntoc_crit                               ! critical N:C ratio below which immobilisation occurrs  
      real :: req                                     ! N required for litter decomposition 
      real :: avl                                     ! mineral N available as NH4
      real :: grossmin                                ! gross mineralisation
      real :: netmin                                  ! net mineralisation
     
      ! xxx try
      real :: soilmoist
      real :: soiltemp


       do lu=1,nlu

         ! XXX try
         soiltemp = 10.    
         soilmoist = 0.8    

         !/////////////////////////////////////////////////////////////////////////
         ! DECAY RATES
         !-------------------------------------------------------------------------
         ! Calculate daily (monthly) decomposition rates as a function of
         ! temperature and moisture
                   
         ! k = k_10 * respir_modifier

         ! (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
         ! from (1),
         ! (2) c = c0*exp(-kt) where c0=initial pool size
         ! from (2), decomposition in any month given by
         ! (3) delta_c = c0 - c0*exp(-k)
         ! from (4)
         ! (4) delta_c = c0*(1.0-exp(-k))
         !-------------------------------------------------------------------------
         do pft=1,npft
           if (lu==lu_category(pft)) then
              
             !-------------------------------------------------------------------------
             ! LITTER TEMPERATURE AND MOISTURE MODIFIER
             ! Temperature: Lloyd & Taylor 1994, Brovkin et al., 2012
             ! Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
             ! Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012
             !-------------------------------------------------------------------------
             
              ! define decomposition rates for current soil temperature and moisture
              klitt_af(pft) = klitt_af10 * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley") ! alternative: "gerten"
              klitt_as(pft) = klitt_as10 * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley") ! alternative: "gerten"
              klitt_bg(pft) = klitt_bg10 * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley") ! alternative: "gerten"
              kexu(pft)     = kexu10     * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley") ! alternative: "gerten"

            end if
          end do

          !-------------------------------------------------------------------------
          ! SOIL TEMPERATURE AND MOISTURE MODIFIER
          ! Temperature: Lloyd & Taylor 1994
          ! Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
          !           Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012
          !-------------------------------------------------------------------------
          ksoil_fs(lu) = ksoil_fs10 * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley")     ! alternative: "gerten"
          ksoil_sl(lu) = ksoil_sl10 * ftemp(soiltemp,"lloyd_and_taylor") * fmoist(soilmoist,"foley")     ! alternative: "gerten"

          do pft=1,npft
            if (lu==lu_category(pft)) then
              
              !////////////////////////////////////////////////////////////////
              ! EXUDATES DECAY
              !----------------------------------------------------------------
              ! Calculate the exudates respiration before litter respiration.
              ! Exudates are mostly short organic compounds (poly- and mono-
              ! saccharides, amino acids, organic acids, phenolic compounds and
              ! enzymes) and are quickly respired and released as CO2.
              !----------------------------------------------------------------                
              dexu = cfrac( 1.-exp(-kexu(pft)), pexud(pft,jpngr) )
              call cmv( dexu, pexud(pft,jpngr), drsoil(lu) )
              

              !////////////////////////////////////////////////////////////////
              ! LITTER DECAY
              ! All goes to daily updated litter decomposition pool
              ! xxx what is justification for PFT-dependent decay constants? xxx
              !----------------------------------------------------------------
              dlitt_af = orgfrac( 1.-exp(-klitt_af(pft)), plitt_af(pft,jpngr) )
              dlitt_as = orgfrac( 1.-exp(-klitt_as(pft)), plitt_as(pft,jpngr) )
              dlitt_bg = orgfrac( 1.-exp(-klitt_bg(pft)), plitt_bg(pft,jpngr) )
                            

              ! Update the litter pools
              call orgmv( dlitt_af, plitt_af(pft,jpngr), dlitt(pft) )
              call orgmv( dlitt_as, plitt_as(pft,jpngr), dlitt(pft) )
              call orgmv( dlitt_bg, plitt_bg(pft,jpngr), dlitt(pft) )
              

              !////////////////////////////////////////////////////////////////
              ! ATMOSPHERIC FRACTION ~ 1 - MICROBIAL GROWTH EFFICIENCY
              ! critical C:N ratio for net mineralisation is a function of C:N
              ! ratio of decomposing litter. Eq. 9 in Xu-Ri & Prentice, 2014
              !----------------------------------------------------------------
              ntoc_crit = ntoc_crit1 * cton( dlitt(pft) ) ** ntoc_crit2
              eff = cton_microb * ntoc_crit
              

              !////////////////////////////////////////////////////////////////
              ! LITTER -> SOIL FLUX AND NET MINERALISATION/IMMOBILISATION
              ! Calculate net mineralisation/immobilisation based on Manzoni
              ! et al. (2008) and Xu-Ri & Prentice (2014).
              !----------------------------------------------------------------

              ! calculate actually required N to satisfy C:N ratio of soil
              ! = eff*delta(litter)*NCratio(soil)
              ! prescribe soil C:N ratio per PFT
              req = 1./ cton_soil(pft) * dlitt(pft)%c%c12 * eff
              
              ! gross N mineralisation 
              grossmin = (1.-eff) * dlitt(pft)%n%n14 ! defined as scalar

              ! net N mineralisation (or immobilisation if negative)
              netmin   = grossmin - req                         ! defined as scalar
              
              ! move fraction 'eff' of C and N from litter to soil
              call orgmv( orgfrac(eff*fastfrac,dlitt(pft)), dlitt(pft), psoil_fs(lu,jpngr) )
              call orgmv( orgfrac(eff*(1.-fastfrac),dlitt(pft)), dlitt(pft), psoil_sl(lu,jpngr) )

              ! move fraction '(1-eff)' of C and N to heterotrophic respiration
              call cmv( cfrac((1.-eff),dlitt(pft)%c), dlitt(pft)%c, drh(lu) )

              if (netmin>0) then                
                ! mineralised litter enters ammonium pool
                nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 + netmin 
                
              else
                ! mineral N is immobilised into soil
                ! if not enough mineral N is available imply fixation (free-living)
                req = (-1.)*netmin
                avl = nh4(lu,jpngr)%n14
                
                if (avl>=req) then
                  ! enough mineral N for immobilisation
                  psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%n%n14 + req*fastfrac
                  psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%n%n14 + req*(1.-fastfrac)             
                  nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - req
                else
                  ! not enough mineral N for immobilisation -> additional free-living fixation
                  psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%n%n14 + req*fastfrac
                  psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%n%n14 + req*(1.-fastfrac)
                  nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - avl
                  dnfix_free(lu) = dnfix_free(lu) + (req - avl)
                end if

              end if
              
            end if
          end do ! pft loop

          
          ! xxx try : un-comment these lines again

          !////////////////////////////////////////////////////////////////
          ! SOIL DECAY
          !----------------------------------------------------------------
          ! Calculate daily/monthly soil decomposition to the atmosphere
!          dsoil_fs = orgfrac( (1.-exp(-ksoil_fs(lu))), psoil_fs(lu,jpngr))
!          dsoil_sl = orgfrac( (1.-exp(-ksoil_sl(lu))), psoil_sl(lu,jpngr))

!          ! Move C from soil pool to heterotrophic respiration
!          call cmv( dsoil_fs%c, psoil_fs(lu,jpngr)%c, drh(lu) )
!          call cmv( dsoil_sl%c, psoil_sl(lu,jpngr)%c, drh(lu) )

!          ! Move N from soil pool to ammonium pool
!          call nmv( dsoil_fs%n, psoil_fs(lu,jpngr)%n, nh4(lu,jpngr) )
!          call nmv( dsoil_sl%n, psoil_sl(lu,jpngr)%n, nh4(lu,jpngr) )


          ! Record monthly (daily) soil turnover flux (labile carbon availability)
          !----------------------------------------------------------------
          ! doc(lu) = dsoil_fs%c%c12 + dsoil_sl%c%c12

        enddo                   !lu                  


       end subroutine littersom
