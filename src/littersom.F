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
      use dailyvars, only: drsoil, drhet, dnfix_free, doc
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
      type (orgpool) :: dsoil_sl                      ! (= cflux_fast_atmos)
      type (orgpool) :: dsoil_fs                      ! (= cflux_fast_atmos)

      ! temporary variables
      real :: eff                                     ! microbial growth efficiency 
      real :: ntoc_crit                               ! critical N:C ratio below which immobilisation occurrs  
      real :: Nreq_B                                  ! N required in litter decomposition to maintain microbial C:N
      real :: Nreq_S                                  ! N required in litter decomposition to maintain SOM C:N
      real :: Nfix                                    ! temporary variable, N fixation implied in litter decomposition,
      real :: rest                                    ! temporary variable
      real :: req                                     ! N required for litter decomposition 
      real :: avl                                     ! mineral N available as NH4
      real :: grossmin                                ! gross mineralisation
      real :: netmin                                  ! net mineralisation
     
      ! xxx try
      real :: soilmoist
      real :: soiltemp

      ! xxx test
      type (nitrogen):: test
      type (carbon)  :: Cbal_before
      type (carbon)  :: Cbal_after
      type (nitrogen):: Nbal_before
      type (nitrogen):: Nbal_after
      real :: Nbal2_before
      real :: Nbal2_after

      ! xxx test
      Cbal_before = cplus( plitt_af(1,1)%c, plitt_as(1,1)%c, plitt_bg(1,1)%c, psoil_fs(1,1)%c, psoil_sl(1,1)%c, drhet(1) )
      Nbal_before = nplus( plitt_af(1,1)%n, plitt_as(1,1)%n, plitt_bg(1,1)%n, psoil_fs(1,1)%n, psoil_sl(1,1)%n, nh4(1,1) )

      do lu=1,nlu

        ! XXX try
        soiltemp = 10.    
        soilmoist = 0.8   

        ! Initialise
        dlitt = orgpool( carbon(0.0), nitrogen(0.0) ) 

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
            ! Calculate the exudates respiration before litter respiration.
            ! Exudates are mostly short organic compounds (poly- and mono-
            ! saccharides, amino acids, organic acids, phenolic compounds and
            ! enzymes) and are quickly respired and released as CO2.
            ! Exudates decay goes to soil respiration 'drsoil'.
            !----------------------------------------------------------------                
            dexu = cfrac( 1.-exp(-kexu(pft)), pexud(pft,jpngr) )
            call cmv( dexu, pexud(pft,jpngr), drsoil(pft) )
            

            !////////////////////////////////////////////////////////////////
            ! LITTER DECAY
            ! All goes to daily updated litter decomposition pool
            ! xxx what is justification for PFT-dependent decay constants? xxx
            !----------------------------------------------------------------
            dlitt_af = orgfrac( 1.-exp(-klitt_af(pft)), plitt_af(pft,jpngr) )
            dlitt_as = orgfrac( 1.-exp(-klitt_as(pft)), plitt_as(pft,jpngr) )
            dlitt_bg = orgfrac( 1.-exp(-klitt_bg(pft)), plitt_bg(pft,jpngr) )

            !print*,'dlitt_af',dlitt_af
            !print*,'dlitt_as',dlitt_as
            !print*,'dlitt_bg',dlitt_bg

            ! Update the litter pools
            call orginit( dlitt(pft) )
            call orgmv( dlitt_af, plitt_af(pft,jpngr), dlitt(pft) )
            call orgmv( dlitt_as, plitt_as(pft,jpngr), dlitt(pft) )
            call orgmv( dlitt_bg, plitt_bg(pft,jpngr), dlitt(pft) )
            

            !////////////////////////////////////////////////////////////////
            ! ATMOSPHERIC FRACTION ~ 1 - MICROBIAL GROWTH EFFICIENCY
            ! critical C:N ratio for net mineralisation is a function of C:N
            ! ratio of decomposing litter. Eq. 9 in Xu-Ri & Prentice, 2014
            !----------------------------------------------------------------
            ntoc_crit = ntoc_crit1 * ntoc( dlitt(pft) ) ** ntoc_crit2  ! = rCR
            eff = ntoc_crit * cton_microb 
            !print*,'ntoc_crit1', ntoc_crit1
            !print*,'ntoc_crit2', ntoc_crit2
            !print*,'(1./cton_microb)',(1./cton_microb)
            !print*,'ntoc litt ',ntoc( dlitt(pft) )
            !print*,'ntoc_crit ',ntoc_crit
            !print*,'eff ',eff
            !stop

            !////////////////////////////////////////////////////////////////
            ! LITTER -> SOIL FLUX AND NET MINERALISATION/IMMOBILISATION
            ! Calculate net mineralisation/immobilisation based on Manzoni
            ! et al. (2008) and Xu-Ri & Prentice (2014).
            !----------------------------------------------------------------    
            ! CARBON LITTER -> SOIL TRANSFER
            !----------------------------------------------------------------    
            ! move fraction 'eff' of C from litter to soil
            call ccp( cfrac( eff*fastfrac     , dlitt(pft)%c ), psoil_fs(lu,jpngr)%c )
            call ccp( cfrac( eff*(1.-fastfrac), dlitt(pft)%c ), psoil_sl(lu,jpngr)%c )

            ! move fraction '(1-eff)' of C to heterotrophic respiration
            call ccp( cfrac( (1.-eff), dlitt(pft)%c ), drhet(lu) )

            ! NITROGEN LITTER -> SOIL TRANSFER
            !----------------------------------------------------------------    
            ! requirement to maintain rB (N:C ratio)
            Nreq_B =  dlitt(pft)%c%c12 * ntoc_crit 

            ! N requirement to maintain rS (SOM N:C ratio)
            Nreq_S = dlitt(pft)%c%c12 * eff / cton_soil(pft)  ! 1/cton_soil = rS

            ! Difference is acquired through N fixation ???
            Nfix = max( 0., Nreq_S - Nreq_B )

            ! If N supply is sufficient, mineralisation occurrs: positive (dNLit-Nreq).
            ! otherwise, immobilisation occurrs: negative (dNLit-Nreq).
            ! Thus, the balance for total organic N is:
            ! dN/dt = -(dNLit - Nreq)
            !       = -(dCLit*rL - dCLit*eff*rB)
            !       = dCLit*(rCR - rL) , rCR=eff*rB ('critical' N:C ratio)
            ! This corresponds to Eq. S3 in Manzoni et al., 2010
            !netmin = dlitt(pft)%n%n14 - Nreq_B XXX COMMENTED THIS OUT BECAUSE IT DOES NOT HOLD SOIL C:N RATIO CONSTANT!!!
            netmin = dlitt(pft)%n%n14 - Nreq_B

            if (netmin<0.0) then
              ! Immobilisation of N
              req = -1.0*netmin
              avl = nh4(lu,jpngr)%n14

              if (avl>=req) then
                ! enough mineral N for immobilisation
                nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - req
              else
                ! not enough N for immobilisation => free-living N fixation
                nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - avl
                Nfix = Nfix + (req - avl)
              endif
            endif

            ! Rest remains in the system: 
            rest = dlitt(pft)%n%n14 - netmin
            call ncp( nfrac( fastfrac     , nitrogen(rest) ), psoil_fs(lu,jpngr)%n )
            call ncp( nfrac( (1.-fastfrac), nitrogen(rest) ), psoil_sl(lu,jpngr)%n )

            ! Add Nfix to SOM ==> should maintain soil C:N ratio
            call ncp( nfrac( fastfrac     , nitrogen(Nfix) ), psoil_fs(lu,jpngr)%n )
            call ncp( nfrac( (1.-fastfrac), nitrogen(Nfix) ), psoil_sl(lu,jpngr)%n )

            ! record this for output (free-living N fixation)
            dnfix_free(lu) = dnfix_free(lu) + Nfix

            ! xxx test
            Nbal2_after = psoil_sl(lu,jpngr)%n%n14 + psoil_fs(lu,jpngr)%n%n14 + nh4(lu,jpngr)%n14 - dnfix_free(lu) 
            !print*,'ok if zero: ', Nbal2_before - Nbal2_after
            
          end if
        end do ! pft loop


        !////////////////////////////////////////////////////////////////
        ! SOIL DECAY
        !----------------------------------------------------------------

        ! xxx test
        !Nbal2_before = psoil_sl(lu,jpngr)%n%n14 + psoil_fs(lu,jpngr)%n%n14 + nh4(lu,jpngr)%n14 

        ! Calculate daily/monthly soil decomposition to the atmosphere
        dsoil_fs = orgfrac( (1.-exp(-ksoil_fs(lu))), psoil_fs(lu,jpngr))
        dsoil_sl = orgfrac( (1.-exp(-ksoil_sl(lu))), psoil_sl(lu,jpngr))

        !! Move C from soil pool to heterotrophic respiration
        call cmv( dsoil_fs%c, psoil_fs(lu,jpngr)%c, drhet(lu) )
        call cmv( dsoil_sl%c, psoil_sl(lu,jpngr)%c, drhet(lu) )

        !! Move N from soil pool to ammonium pool
        call nmv( dsoil_fs%n, psoil_fs(lu,jpngr)%n, nh4(lu,jpngr) )
        call nmv( dsoil_sl%n, psoil_sl(lu,jpngr)%n, nh4(lu,jpngr) )

        ! xxx test
        !Nbal2_after = psoil_sl(lu,jpngr)%n%n14 + psoil_fs(lu,jpngr)%n%n14 + nh4(lu,jpngr)%n14
        !print*,'ok if zero: ', Nbal2_before - Nbal2_after

        ! Record monthly (daily) soil turnover flux (labile carbon availability)
        !----------------------------------------------------------------
        ! doc(lu) = dsoil_fs%c%c12 + dsoil_sl%c%c12

      enddo                   !lu

      ! xxx test
      Cbal_after = cplus( plitt_af(1,1)%c, plitt_as(1,1)%c, plitt_bg(1,1)%c, psoil_fs(1,1)%c, psoil_sl(1,1)%c, drhet(1) )
      Nbal_after = nplus( plitt_af(1,1)%n, plitt_as(1,1)%n, plitt_bg(1,1)%n, psoil_fs(1,1)%n, psoil_sl(1,1)%n, nh4(1,1) )
      Nbal_after = nminus( Nbal_after, nitrogen(dnfix_free(lu)) )

      !print*,'ok if zero: ', cminus( Cbal_after, Cbal_before )
      !print*,'dnfix_free ', dnfix_free
      !print*,'ok if zero: ', nminus( Nbal_after, Nbal_before )


      end subroutine littersom
