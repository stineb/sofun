module littersom_lpj
  !////////////////////////////////////////////////////////////////
  ! LPJ LITTERSOM MODULE
  ! Contains the "main" subroutine 'littersom' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'littersom' must contain this list 
  ! of subroutines (names that way).
  !   - getpar_modl_littersom
  !   - initio_littersom
  !   - initoutput_littersom
  !   - getout_daily_littersom
  !   - getout_monthly_littersom
  !   - writeout_ascii_littersom
  !   - littersom
  ! Required module-independent model state variables (necessarily 
  ! updated by 'littersom') are:
  !   - xxx
  !----------------------------------------------------------------
  use params_core, only: npft, maxgrid, nmonth, nlu, recycle, ndayyear, &
    pft_start, pft_end
  use params_modl, only: islu, cton_soil
  use params_siml, only: spinup
  use classdefs
  use rates
  use statevars, only: ncpa_d
  use fluxes, only: drsoil, drhet, dcex
  use statevars, only: dtemp_soil, ddoc
  use pools, only: pexud, plitt_af, plitt_as, plitt_bg, psoil_fs, &
    psoil_sl, pninorg

  implicit none
  ! private 
  ! public getpar_modl_littersom, initio_littersom, initoutput_littersom, &
  !  getout_daily_littersom, getout_monthly_littersom, writeout_ascii_littersom, &
  !  littersom


  ! LITTERSOM PARAMETERS
  ! read from file, therefore not declared as parameter
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


  ! MODULE-SPECIFIC VARIABLES
  ! These are not required outside module 'littersom', but are used in different SRs of this module
  ! real :: xxx                    ! xxx

  ! DAILY OUTPUT VARIABLES
  ! xxx add jpngr dim
  real, dimension(nlu,ndayyear,maxgrid)  :: outdnetmin
  real, dimension(nlu,ndayyear,maxgrid)  :: outdnfixfree  ! N fixation by free-living organisms [gN/m2/d]
  ! real, dimension(nlu,ndayyear,maxgrid)  :: outdNimmo
  real, dimension(npft,ndayyear,maxgrid) :: outdcex

  real, dimension(npft,maxgrid) :: outaClitt
  real, dimension(nlu,maxgrid)  :: outaCsoil
  real, dimension(nlu,maxgrid)  :: outanreq      ! N required from litter -> soil transfer [gN/m2/d]
  real, dimension(npft,maxgrid) :: outaClit2soil
  real, dimension(npft,maxgrid) :: outaNlit2soil
  real, dimension(nlu,maxgrid)  :: outaCdsoil
  real, dimension(nlu,maxgrid)  :: outaNdsoil
  real, dimension(nlu,maxgrid)  :: outaNimmo

  ! MONTHLY OUTPUT VARIABLES
  ! xxx add jpngr dim
  real, dimension(nmonth,maxgrid) :: xxx


contains

  subroutine littersom( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Litter and SOM decomposition and nitrogen mineralisation.
    !  1st order decay of litter and SOM pools, governed by temperature
    !  and soil moisture following LPJ (Sitch et al., 2003) and 
    !  Xu-Ri & Prentice (XXX).
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------

    ! ARGUMENTS
    integer, intent(in) :: jpngr                    ! grid cell number
    integer, intent(in) :: doy                      ! day of year

    ! LOCAL VARIABLES
    integer :: lu                                   ! counter variable for landuse class
    integer :: pft                                  ! counter variable for PFT number

    ! temperature/soil moisture-modified decay constants
    real, dimension(nlu)  :: klitt_af               ! decay rate, above-ground fast (leaf) litter (= k_litter_leaf)
    real, dimension(nlu)  :: klitt_as               ! decay rate, above-ground slow (woody) litter (= k_litter_woody)
    real, dimension(nlu)  :: klitt_bg               ! decay rate, below-ground fast litter (= k_litter_root)
    real, dimension(nlu)  :: kexu                   ! decay rate, exudates (= k_exu)
    real, dimension(nlu)  :: ksoil_fs               ! decay rate, fast soil (= k_fast)
    real, dimension(nlu)  :: ksoil_sl               ! decay rate, slow soil (= k_slow)

    ! temporary pools
    type(carbon)  :: dexu                           ! exudates decomposed in time step (= exu_decom)
    type(orgpool) :: dlitt                          ! total litter decomposed per time step
    type(orgpool), dimension(npft) :: dlitt_af      ! above-ground fast litter decomposed per time step (= litterdag_fast)
    type(orgpool), dimension(npft) :: dlitt_as      ! above-ground slow litter decomposed per time step (= litterdag_slow)
    type(orgpool), dimension(npft) :: dlitt_bg      ! below-ground slow litter decomposed per time step (= litter_decom_bg)
    type(orgpool) :: dsoil_sl                       ! (= cflux_fast_atmos)
    type(orgpool) :: dsoil_fs                       ! (= cflux_fast_atmos)
   
    ! temporary variables
    real :: cton_soil_local
    real :: eff                                     ! microbial growth efficiency 
    real :: ntoc_crit                               ! critical N:C ratio below which immobilisation occurrs  
    real :: Nreq_S                                  ! N required in litter decomposition to maintain SOM C:N
    real :: Nfix                                    ! temporary variable, N fixation implied in litter decomposition,
    real :: rest                                    ! temporary variable
    real :: req                                     ! N required for litter decomposition 
    real :: avl                                     ! mineral N available as inorganic N
    real :: netmin_litt                             ! net N mineralisation from litter decomposition

    integer, save :: invocation = 0                 ! internally counted simulation year
    integer, parameter :: spinupyr_soilequil_1 = 600   ! year of analytical soil equilibration, based on mean litter -> soil input flux
    integer, parameter :: spinupyr_soilequil_2 = 1200  ! year of analytical soil equilibration, based on mean litter -> soil input flux
    integer, parameter :: spinupyr_phaseinit_2 = 900
    integer, parameter :: spinupyr_phaseinit_3 = 300

    ! real :: acc                                     ! soil equilibration acceleration factor
    ! real :: scal, hi, lo
    real, dimension(nlu,maxgrid), save :: mean_insoil_fs
    real, dimension(nlu,maxgrid), save :: mean_insoil_sl
    real, dimension(nlu,maxgrid), save :: mean_ksoil_sl
    real, dimension(nlu,maxgrid), save :: mean_ksoil_fs
    real :: ntoc_save_fs, ntoc_save_sl

#if _check_sanity
    real :: cbal_before
    real :: cbal_after
    real :: nbal_before
    real :: nbal_after

    cbal_before = plitt_af(1,1)%c%c12 + plitt_as(1,1)%c%c12 &
      + plitt_bg(1,1)%c%c12 + psoil_fs(1,1)%c%c12 &
      + psoil_sl(1,1)%c%c12 + drhet(1)%c12
    nbal_before = plitt_af(1,1)%n%n14 + plitt_as(1,1)%n%n14 &
      + plitt_bg(1,1)%n%n14 + psoil_fs(1,1)%n%n14 &
      + psoil_sl(1,1)%n%n14 + pninorg(1,1)%n14
#endif

    !-------------------------------------------------------------------------
    ! Count number of calls (one for each simulation year)
    !-------------------------------------------------------------------------
    if (doy==1) invocation = invocation + 1

    ! initialise 
    if (invocation==1 .and. doy==1) then
      mean_insoil_fs(:,:) = 0.0
      mean_insoil_sl(:,:) = 0.0
      mean_ksoil_fs(:,:)  = 0.0
      mean_ksoil_sl(:,:)  = 0.0
    end if

    !-------------------------------------------------------------------------
    ! Set soil turnover accelerator for equilibration during spinup.
    ! Evenly scale soil inputs and soil decomposition rates during
    ! first 200 years with a scalar that linearly decreases from 200
    ! to 1 over the first 200 years of the simulation.
    ! Value 200 is chosen for quick equilibration without overshooting
    ! for a temperate climate (Switzerland). May have to adjust this
    ! for improving performance with a global simulation.
    !-------------------------------------------------------------------------
    ! if (spinup) then
    !   acc = max( 1.0, 200.0 - real(invocation) ) 
    ! else
    !   acc = 1.0
    ! end if
    ! acc = 1.0

    do lu=1,nlu

      ! if ( abs( cton(psoil_fs(lu,jpngr)) - cton_soil(1) ) > 1e-5 ) stop 'A fs: C:N not ok'
      ! if ( abs( cton(psoil_sl(lu,jpngr)) - cton_soil(1) ) > 1e-5 ) stop 'A sl: C:N not ok'

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
      ! (4) delta_c = c0*(1.00-exp(-k))
      !-------------------------------------------------------------------------
      do pft=1,npft
        if (islu(pft,lu)) then
            
          !-------------------------------------------------------------------------
          ! LITTER TEMPERATURE AND MOISTURE MODIFIER
          ! Temperature: Lloyd & Taylor 1994, Brovkin et al., 2012
          ! Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
          ! Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012
          !-------------------------------------------------------------------------
           
          ! define decomposition rates for current soil temperature and moisture
          klitt_af(lu) = klitt_af10 &
            * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
            * fmoist(ncpa_d(lu),"foley") ! alternative: "gerten"
          klitt_as(lu) = klitt_as10 &
            * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
            * fmoist(ncpa_d(lu),"foley") ! alternative: "gerten"
          klitt_bg(lu) = klitt_bg10 &
            * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
            * fmoist(ncpa_d(lu),"foley") ! alternative: "gerten"
          kexu(lu)     = kexu10     &
            * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
            * fmoist(ncpa_d(lu),"foley") ! alternative: "gerten"

        end if
      end do

      !-------------------------------------------------------------------------
      ! SOIL TEMPERATURE AND MOISTURE MODIFIER
      ! Temperature: Lloyd & Taylor 1994
      ! Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
      !           Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012
      !-------------------------------------------------------------------------
      ksoil_fs(lu) = ksoil_fs10 &
        * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
        * fmoist(ncpa_d(lu),"foley")     ! alternative: "gerten"
      ksoil_sl(lu) = ksoil_sl10 &
        * ftemp(dtemp_soil(lu,jpngr),"lloyd_and_taylor") &
        * fmoist(ncpa_d(lu),"foley")     ! alternative: "gerten"

      !-------------------------------------------------------------------------
      ! Initialisation of decomposing pool 
      ! (temporary, decomposition for each LU-class).
      !-------------------------------------------------------------------------
      dlitt = orgpool( carbon(0.0), nitrogen(0.0) ) 

      !////////////////////////////////////////////////////////////////
      ! LITTER DECAY
      ! All goes to daily updated litter decomposition pool
      !----------------------------------------------------------------
      do pft=1,npft
        if (islu(pft,lu)) then
                        
          dlitt_af = orgfrac( 1.0-exp(-klitt_af(lu)), plitt_af(pft,jpngr) )
          dlitt_as = orgfrac( 1.0-exp(-klitt_as(lu)), plitt_as(pft,jpngr) )
          dlitt_bg = orgfrac( 1.0-exp(-klitt_bg(lu)), plitt_bg(pft,jpngr) )

          ! Update the litter pools
          call orgmv( dlitt_af(pft), plitt_af(pft,jpngr), dlitt )
          call orgmv( dlitt_as(pft), plitt_as(pft,jpngr), dlitt )
          call orgmv( dlitt_bg(pft), plitt_bg(pft,jpngr), dlitt )

          ! acceleration of soil equilibration 
          ! (scale soil inputs and soil decomposition rates)
          ! dlitt = orgfrac( acc, dlitt )

          ! ! tried this and it doesn't work either
          ! scal = max( 0.0, 1.0 - real(invocation)/300.0 )
          ! if (spinup.and.invocation<300) then
          !   hi = dlitt%c%c12 / cton_soil(pft)
          !   lo = dlitt%n%n14
          !   dlitt%n%n14 = scal * hi + (1.0-scal) * lo
          ! end if


          !print*,'slow litter pool, C:N', plitt_as(pft,jpngr), cton( plitt_as(pft,jpngr), default=0.0 )
          !print*,'fast litter pool, C:N', plitt_af(pft,jpngr), cton( plitt_af(pft,jpngr), default=0.0 )
          !print*,'bg   litter pool, C:N', plitt_bg(pft,jpngr), cton( plitt_bg(pft,jpngr), default=0.0 )

          !print*,'pleaf, C:N           ', pleaf(pft,jpngr), cton( pleaf(pft,jpngr), default=0.0 )
          !print*,'proot, C:N           ', proot(pft,jpngr), cton( proot(pft,jpngr), default=0.0 )
      
          end if
        end do

        !----------------------------------------------------------------
        ! Soil C:N ratio is the average PFT-specific prescribed C:N ratio
        ! weighted by the PFT-specific decomposition.
        !----------------------------------------------------------------
        cton_soil_local = sum( cton_soil(pft_start(lu):pft_end(lu)) &
          * ( dlitt_af(pft_start(lu):pft_end(lu))%c%c12  &
            + dlitt_as(pft_start(lu):pft_end(lu))%c%c12  &
            + dlitt_bg(pft_start(lu):pft_end(lu))%c%c12 ) ) / dlitt%c%c12


        !////////////////////////////////////////////////////////////////
        ! ATMOSPHERIC FRACTION ~ 1 - MICROBIAL GROWTH EFFICIENCY
        ! critical C:N ratio for net mineralisation is a function of C:N
        ! ratio of decomposing litter. Eq. 9 in Xu-Ri & Prentice, 2014
        !----------------------------------------------------------------
        ntoc_crit = ntoc_crit1 * ntoc( dlitt, default=0.0 ) ** ntoc_crit2  ! = rCR
        eff = ntoc_crit * cton_soil_local

        ! write(0,*) 'ntoc_crit1', ntoc_crit1
        ! write(0,*) 'ntoc_crit2', ntoc_crit2
        ! write(0,*) 'cton litt ', cton( dlitt )
        ! write(0,*) 'cton_crit ', 1/ntoc_crit
        ! write(0,*) 'eff ',eff
        !stop

        !////////////////////////////////////////////////////////////////
        ! LITTER -> SOIL FLUX AND NET MINERALISATION/IMMOBILISATION
        ! Calculate net mineralisation/immobilisation based on Manzoni
        ! et al. (2008) and Xu-Ri & Prentice (2014).
        !----------------------------------------------------------------    
        ! CARBON LITTER -> SOIL TRANSFER
        !----------------------------------------------------------------    

        ! record N:C ratio to override later (compensating for numerical imprecision)
        ntoc_save_fs = ntoc( psoil_fs(lu,jpngr), default=0.0 )
        ntoc_save_sl = ntoc( psoil_sl(lu,jpngr), default=0.0 )

        ! move fraction 'eff' of C from litter to soil
        call ccp( cfrac( eff*fastfrac      , dlitt%c ), psoil_fs(lu,jpngr)%c )
        call ccp( cfrac( eff*(1.0-fastfrac), dlitt%c ), psoil_sl(lu,jpngr)%c )

        ! move fraction '(1-eff)' of C to heterotrophic respiration
        call ccp( cfrac( (1.0-eff), dlitt%c ), drhet(lu) )

        ! get average litter -> soil flux for analytical soil C equilibration
        if ( spinup .and. invocation > ( spinupyr_soilequil_1 - recycle ) .and. invocation < spinupyr_soilequil_1 &
          .or. spinup .and. invocation > ( spinupyr_soilequil_2 - recycle ) .and. invocation < spinupyr_soilequil_2) then
          mean_insoil_fs(lu,jpngr) = mean_insoil_fs(lu,jpngr) + eff * fastfrac * dlitt%c%c12
          mean_insoil_sl(lu,jpngr) = mean_insoil_sl(lu,jpngr) + eff * (1.0-fastfrac) * dlitt%c%c12
        end if

        ! OUTPUT COLLECTION
        outaClit2soil(pft,jpngr) = outaClit2soil(pft,jpngr) + dlitt%c%c12 * eff
        ! outaNlit2soil(pft,jpngr) = outaNlit2soil(pft,jpngr) + dlitt%n%n14
        outaNlit2soil(lu,jpngr)  = outaNlit2soil(lu,jpngr)  + Nreq_S


        !----------------------------------------------------------------    
        ! N MINERALISATION
        !----------------------------------------------------------------    
        ! N requirement to maintain rS (SOM N:C ratio)
        Nreq_S = dlitt%c%c12 * eff / cton_soil_local  ! 1/cton_soil = rS

        ! write(0,*) 'cton_soil      ',cton_soil_local 
        ! write(0,*) 'cton_crit      ',1.0/ntoc_crit

        ! OUTPUT COLLECTION
        outanreq(lu,jpngr) = outanreq(lu,jpngr) + Nreq_S

        !----------------------------------------------------------------    
        ! If N supply is sufficient, mineralisation occurrs: positive (dNLit-Nreq).
        ! otherwise, immobilisation occurrs: negative (dNLit-Nreq).
        ! Thus, the balance for total organic N is:
        ! dN/dt = -(dNLit - Nreq)
        !       = -(dCLit*rL - dCLit*eff*rS)
        !       = dCLit*(rCR - rL) , rCR=eff*rS ('critical' N:C ratio)
        ! This corresponds to Eq. S3 in Manzoni et al., 2010, but ...
        ! rS takes the place of rB.
        !----------------------------------------------------------------    

        ! xxx try:
        ! >>>>>>>>>>>
        ! ! Trick for soil C quilibration: Add projected N mineralisation from 
        ! ! soil decomposition during spinup phase I. Like in LPX.
        ! if (spinup.and.invocation<=spinupyr_phaseinit_3 ) then
        !   netmin_litt = dlitt%c%c12 / cton_soil_local - Nreq_S
        ! else
        !   netmin_litt = dlitt%n%n14 - Nreq_S
        ! end if
        ! ===========
        netmin_litt = dlitt%n%n14 - Nreq_S
        ! <<<<<<<<<<<

        ! write(0,*) 'netmin_litt', netmin_litt
        ! write(0,*) 'immo direct', dlitt%c%c12 * ( ntoc_crit - dlitt%n%n14 / dlitt%c%c12 )

        ! OUTPUT COLLECTION
        outdnetmin(lu,doy,jpngr) = outdnetmin(lu,doy,jpngr) + netmin_litt
        ! outdNimmo(lu,doy,jpngr)  = outdNimmo(lu,doy,jpngr)  - netmin_litt  ! minus because netmin_litt < 0 for immobilisation
        outaNimmo(lu,jpngr)      = outaNimmo(lu,jpngr)      - netmin_litt  ! minus because netmin_litt < 0 for immobilisation


        if (netmin_litt>0.0) then
          !----------------------------------------------------------------    
          ! net N mineralisation
          !----------------------------------------------------------------    
          pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + netmin_litt
          ! stop 'net N mineralisation from litter decomposition'

        else

          ! xxx try:
          ! >>>>>>>>>>>
          if ( spinup .and. invocation<=spinupyr_phaseinit_2 ) then
            !----------------------------------------------------------------    
            ! N fixation by free-living bacteria in litter to prevent immo-
            ! bilisation.
            !----------------------------------------------------------------    
            req = -1.0 * netmin_litt
            Nfix = req
            req = 0.0

          !   ! ===========
          else
            ! xxx THIS LEADS TO COMPLETE DIE-OFF: equilibration without immobilisation
            !----------------------------------------------------------------    
            ! immobilisation
            !----------------------------------------------------------------    
            req = -1.0 * netmin_litt
            avl = pninorg(lu,jpngr)%n14

            ! write(0,*) 'req, avl',req,avl

            if (avl>=req) then
              ! enough mineral N for immobilisation
              pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - req
              req = 0.0
            else
              ! not enough pninorg for immobilisation
              pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - avl
              req = req - avl

              !----------------------------------------------------------------    
              ! N fixation by free-living bacteria in litter to satisfy remainder
              !----------------------------------------------------------------    
              Nfix = req
              req = 0.0

            end if
            ! <<<<<<<<<<<
          end if

        endif

        ! Nreq_S (= dlitt - netmin) remains in the system: 
        call ncp( nfrac( fastfrac      , nitrogen(Nreq_S) ), psoil_fs(lu,jpngr)%n )
        call ncp( nfrac( (1.0-fastfrac), nitrogen(Nreq_S) ), psoil_sl(lu,jpngr)%n )

        ! Prevent accumulating deviation of soil C:N ratio due to numerical imprecision.
        ! Warning: This may not strictly conserve mass!
        if (ntoc_save_fs>0.0) then
          psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%c%c12 * ntoc_save_fs
        end if
        if (ntoc_save_sl>0.0) then
          psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%c%c12 * ntoc_save_sl
        end if

        if ( abs( cton(psoil_fs(lu,jpngr)) - cton_soil(1) ) > 1e-5 ) stop 'B fs: C:N not ok'
        if ( abs( cton(psoil_sl(lu,jpngr)) - cton_soil(1) ) > 1e-5 ) stop 'B sl: C:N not ok'

        ! OUTPUT COLLECTION
        outdnfixfree(lu,doy,jpngr) = outdnfixfree(lu,doy,jpngr) + Nfix
          
        ! C:N ratio of soil influx
        if ( abs( dlitt%c%c12 * eff / Nreq_S - cton_soil(1) ) > 1e-5 ) stop 'imprecision'


        do pft=1,npft
          if (islu(pft,lu)) then

          !////////////////////////////////////////////////////////////////
          ! EXUDATES DECAY
          ! Calculate the exudates respiration before litter respiration.
          ! Exudates are mostly short organic compounds (poly- and mono-
          ! saccharides, amino acids, organic acids, phenolic compounds and
          ! enzymes) and are quickly respired and released as CO2.
          ! Exudates decay goes to soil respiration 'drsoil'.
          ! This is executed after litter mineralisation as N fixation by 
          ! free-living bacteria is driven by exudates availability.
          !----------------------------------------------------------------                
          dexu = cfrac( 1.0-exp(-kexu(lu)), pexud(pft,jpngr) )
          call cmv( dexu, pexud(pft,jpngr), drsoil(pft) )

        end if
      end do ! pft loop


      !////////////////////////////////////////////////////////////////
      ! SOIL DECAY
      !----------------------------------------------------------------
      ! Calculate daily/monthly soil decomposition to the atmosphere

      ! record N:C ratio to override later (compensating for numerical imprecision)
      ntoc_save_fs = ntoc( psoil_fs(lu,jpngr), default=0.0 )
      ntoc_save_sl = ntoc( psoil_sl(lu,jpngr), default=0.0 )

      ! xxx debug
      dsoil_fs%c%c12 = (1.0-exp(-ksoil_fs(lu))) * psoil_fs(lu,jpngr)%c%c12
      dsoil_sl%c%c12 = (1.0-exp(-ksoil_sl(lu))) * psoil_sl(lu,jpngr)%c%c12

      dsoil_fs%n%n14 = dsoil_fs%c%c12 / cton_soil_local
      dsoil_sl%n%n14 = dsoil_sl%c%c12 / cton_soil_local

      if ( abs( cton(dsoil_fs) - cton_soil_local ) > 1e-5 ) stop 'dsoil imprecision fs'
      if ( abs( cton(dsoil_sl) - cton_soil_local ) > 1e-5 ) stop 'dsoil imprecision sl'

      ! carbon
      psoil_fs(lu,jpngr)%c%c12 = psoil_fs(lu,jpngr)%c%c12 - dsoil_fs%c%c12
      psoil_sl(lu,jpngr)%c%c12 = psoil_sl(lu,jpngr)%c%c12 - dsoil_sl%c%c12
      
      ! to heterotrophic respiration
      drhet(lu)%c12 = drhet(lu)%c12 + dsoil_fs%c%c12 + dsoil_sl%c%c12

      ! nitrogen
      psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%n%n14 - dsoil_fs%n%n14
      psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%n%n14 - dsoil_sl%n%n14

      if ( abs( cton( psoil_fs(lu,jpngr)) - cton_soil_local ) > 1e-5 ) stop 'C fs: C:N not ok'
      if ( abs( cton( psoil_sl(lu,jpngr)) - cton_soil_local ) > 1e-5 ) stop 'C sl: C:N not ok'

      ! Prevent accumulating deviation of soil C:N ratio due to numerical imprecision.
      ! Warning: this does not strictly conserve mass!
      psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%c%c12 / cton_soil_local
      psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%c%c12 / cton_soil_local


      ! xxx try:
      ! >>>>>>>>>>>
      ! ! to inorganic N pool
      ! pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + dsoil_fs%n%n14 + dsoil_sl%n%n14
      ! ===========
      if ( spinup .and. invocation <= spinupyr_phaseinit_3 ) then
        pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + eff * dlitt%c%c12 / cton_soil_local
        ! write(0,*) 'fraction immobilised', (-1.0)*netmin_litt / ( eff * dlitt%c%c12 / cton_soil_local )
        ! ! xxx try:
        ! pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - 0.02 * (eff * dlitt%c%c12 / cton_soil_local)
      else  
        pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + dsoil_fs%n%n14 + dsoil_sl%n%n14
        ! write(0,*) 'fraction immobilised', (-1.0)*netmin_litt / ( dsoil_fs%n%n14 + dsoil_sl%n%n14 )
        ! ! xxx try:
        ! pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - 0.02 * ( dsoil_fs%n%n14 + dsoil_sl%n%n14 )
      end if
      ! <<<<<<<<<<<      



      ! print*,'cton(dsoil_sl) ',cton(dsoil_sl)
      ! print*,'cton(dsoil_fs) ',cton(dsoil_fs)
      ! print*,'cton_soil_local',cton_soil_local

      ! get average litter -> soil flux for analytical soil C equilibration
      if ( spinup .and. invocation > ( spinupyr_soilequil_1 - recycle ) .and. invocation<spinupyr_soilequil_1 &
        .or. spinup .and. invocation > ( spinupyr_soilequil_2 - recycle ) .and. invocation<spinupyr_soilequil_2 ) then
        mean_ksoil_fs(lu,jpngr) = mean_ksoil_fs(lu,jpngr) + ksoil_fs(lu)
        mean_ksoil_sl(lu,jpngr) = mean_ksoil_sl(lu,jpngr) + ksoil_sl(lu)
      end if

      ! analytical soil C equilibration
      if ( spinup .and. invocation==spinupyr_soilequil_1 .and. doy==ndayyear &
        .or. spinup .and. invocation==spinupyr_soilequil_2 .and. doy==ndayyear ) then
        psoil_fs(lu,jpngr)%c%c12 = mean_insoil_fs(lu,jpngr) / mean_ksoil_fs(lu,jpngr)
        psoil_sl(lu,jpngr)%c%c12 = mean_insoil_sl(lu,jpngr) / mean_ksoil_sl(lu,jpngr)
        psoil_fs(lu,jpngr)%n%n14 = psoil_fs(lu,jpngr)%c%c12 * ntoc_save_fs
        psoil_sl(lu,jpngr)%n%n14 = psoil_sl(lu,jpngr)%c%c12 * ntoc_save_sl
        mean_insoil_fs(lu,jpngr) = 0.0
        mean_insoil_sl(lu,jpngr) = 0.0
        mean_ksoil_fs(lu,jpngr)  = 0.0
        mean_ksoil_sl(lu,jpngr)  = 0.0
      end if

      
      ! OUTPUT COLLECTION
      outdnetmin(lu,doy,jpngr) = outdnetmin(lu,doy,jpngr) + dsoil_fs%n%n14 + dsoil_sl%n%n14
      outaCdsoil(lu,jpngr)     = outaCdsoil(lu,jpngr)     + dsoil_fs%c%c12 + dsoil_sl%c%c12
      outaNdsoil(lu,jpngr)     = outaNdsoil(lu,jpngr)     + dsoil_fs%n%n14 + dsoil_sl%n%n14

      ! if (doy==ndayyear) then
      !   print*,'outaNdsoil(lu,jpngr)                   ',outaNdsoil(lu,jpngr)
      !   print*,'outaClit2soil(lu,jpngr)/cton_soil_local',outaClit2soil(lu,jpngr)/cton_soil_local
      ! end if     

      ! print*,'eff * dlitt%c', eff * dlitt%c%c12
      ! print*,'dsoil%c      ', dsoil_fs%c%c12 + dsoil_sl%c%c12
      ! print*,'cton_soil_local',cton_soil_local
      ! print*,'doy',doy
      ! print*,'cton(soil)     ',cton( orgplus(psoil_fs(lu,jpngr),psoil_sl(lu,jpngr)))
      ! print*,'cton(soil fast)',cton( psoil_fs(lu,jpngr) )
      ! print*,'cton(soil slow)',cton( psoil_sl(lu,jpngr) )

      ! Record monthly (daily) soil turnover flux (labile carbon availability)
      !----------------------------------------------------------------
      ddoc(lu) = dsoil_fs%c%c12 + dsoil_sl%c%c12

    enddo                   !lu

#if _check_sanity
    cbal_after = plitt_af(1,1)%c%c12 + plitt_as(1,1)%c%c12 &
      + plitt_bg(1,1)%c%c12 + psoil_fs(1,1)%c%c12 &
      + psoil_sl(1,1)%c%c12 + drhet(1)%c12
    nbal_after = plitt_af(1,1)%n%n14 + plitt_as(1,1)%n%n14 &
      + plitt_bg(1,1)%n%n14 + psoil_fs(1,1)%n%n14 &
      + psoil_sl(1,1)%n%n14 + pninorg(1,1)%n14
    nbal_after = Nbal_after - dnfix_free(1)
    if (abs(cbal_after-cbal_before)>1.0d9) then
      print*,'cbal_before, cbal_after ', cbal_before, cbal_after
      stop 'C balance violated in LITTERSOM'
    endif
    if (abs(nbal_after-nbal_before)>1.0d9) then
      print*,'nbal_before, nbal_after ', nbal_before, nbal_after
      stop 'N balance violated in LITTERSOM'
    endif
#endif

  end subroutine littersom


  subroutine getpar_modl_littersom( paramfilnam )
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads littersom module-specific parameters 
    ! from input file
    !----------------------------------------------------------------

    ! ARGUMENTS
    character(len=*) paramfilnam

    ! LOCAL VARIABLES
    integer, parameter :: npar = 10
    real, dimension(npar) :: params_array

    open(unit=06,file=trim(paramfilnam),status='OLD')      
    read (06,*) params_array
    close (06)

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

  end subroutine getpar_modl_littersom


  subroutine initio_littersom( prefix )
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------

    ! ARGUMENTS
    character(len=*) :: prefix

    ! LOCAL VARIABLES
    character(len=256) :: filnam


    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------

    ! C EXUDATION
    filnam=trim(prefix)//'.d.cex.out'
    open(105,file=filnam,err=999,status='unknown')

    ! NET N MINERALISATION
    filnam=trim(prefix)//'.d.netmin.out'
    open(106,file=filnam,err=999,status='unknown')

    ! BIOLOGICAL NITROGEN FIXATION OF FREE-LIVING  ORGANISMS
    filnam=trim(prefix)//'.d.nfixfree.out'
    open(108,file=filnam,err=999,status='unknown')



    !----------------------------------------------------------------
    ! MONTHLY OUTPUT
    !----------------------------------------------------------------


    !----------------------------------------------------------------
    ! ANNUAL OUTPUT
    !----------------------------------------------------------------

    ! LITTER C
    filnam=trim(prefix)//'.a.clitt.out'
    open(301,file=filnam,err=999,status='unknown')

    ! SOIL C
    filnam=trim(prefix)//'.a.csoil.out'
    open(302,file=filnam,err=999,status='unknown')

    ! N REQUIRED FOR LITTER -> SOIL TRANSFER
    filnam=trim(prefix)//'.a.nreq.out'
    open(304,file=filnam,err=999,status='unknown')

    ! C LITTER -> SOIL TRANSFER
    filnam=trim(prefix)//'.a.clit2soil.out'
    open(305,file=filnam,err=999,status='unknown')

    ! N LITTER -> SOIL TRANSFER
    filnam=trim(prefix)//'.a.nlit2soil.out'
    open(306,file=filnam,err=999,status='unknown')

    ! C MINERALISATION FROM SOIL DECOMPOSITION
    filnam=trim(prefix)//'.a.cdsoil.out'
    open(313,file=filnam,err=999,status='unknown')

    ! N MINERALISATION FROM SOIL DECOMPOSITION
    filnam=trim(prefix)//'.a.ndsoil.out'
    open(314,file=filnam,err=999,status='unknown')

    ! N IMMOBILISATION FROM LITTER DECOMPOSITION
    filnam=trim(prefix)//'.a.nimmo.out'
    open(315,file=filnam,err=999,status='unknown')

    return

  888  stop 'INITIO_littersom: error opening output files'
  999  stop 'INITIO: error opening output files'

  end subroutine initio_littersom


  subroutine initoutput_littersom
    !////////////////////////////////////////////////////////////////
    !  Initialises littersomance-specific output variables
    !----------------------------------------------------------------

    ! xxx remove their day-dimension
    outaClitt(:,:)      = 0.0
    outaCsoil(:,:)      = 0.0
    outanreq(:,:)       = 0.0
    outaClit2soil(:,:)  = 0.0
    outaNlit2soil(:,:)  = 0.0
    outaNdsoil(:,:)     = 0.0
    outaCdsoil(:,:)     = 0.0
    outaNimmo(:,:)      = 0.0

    outdcex(:,:,:)      = 0.0
    outdnetmin(:,:,:)   = 0.0
    ! outdNimmo(:,:,:)    = 0.0
    outdnfixfree(:,:,:) = 0.0


  end subroutine initoutput_littersom


  subroutine getout_daily_littersom( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use pools, only: dwtot

    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    ! Save the daily totals:
    ! xxx add lu-dimension and jpngr-dimension
    outdcex(:,doy,jpngr)   = dcex(:)

  end subroutine getout_daily_littersom


  ! subroutine getout_monthly_littersom( jpngr, moy )
  !   !////////////////////////////////////////////////////////////////
  !   !  SR called daily to sum up output variables.
  !   !----------------------------------------------------------------
  !   implicit none

  !   ! ARGUMENTS
  !   integer, intent(in) :: jpngr
  !   integer, intent(in) :: moy

  !   ! Calculate other monthly totals:
  !   outmxxx(:,moy,jpngr) = outxxx(:,moy,jpngr)

  !   return  

  ! end subroutine getout_monthly_littersom


  subroutine getout_annual_littersom( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use params_modl, only: npft, lu_category
    use pools

    ! ARGUMENTS
    integer, intent(in) :: jpngr

    outaClitt(:,jpngr) = plitt_af(:,jpngr)%c%c12 + plitt_as(:,jpngr)%c%c12 + plitt_bg(:,jpngr)%c%c12
    outaCsoil(:,jpngr) = psoil_sl(:,jpngr)%c%c12 + psoil_fs(:,jpngr)%c%c12

  end subroutine getout_annual_littersom


  subroutine writeout_ascii_littersom( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE littersom-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use params_core, only: ndayyear, nmonth
    use params_siml, only: firstyeartrend, spinupyears

    ! Arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! Local variables
    real :: itime
    integer :: day, moy, jpngr
    
    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1


    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if (.not.spinup) then
      ! Write daily output only during transient simulation
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

        if (nlu>1) stop 'writeout_ascii_littersom: write out lu-area weighted sum'

        ! xxx lu-area weighted sum if nlu>0
        write(105,999) itime, sum(outdcex(:,day,jpngr))
        write(106,999) itime, sum(outdnetmin(:,day,jpngr))
        ! write(109,999) itime, sum(outdNimmo(:,day,jpngr))
        write(108,999) itime, sum(outdnfixfree(:,day,jpngr))

      end do
    end if


    !-------------------------------------------------------------------------
    ! MONTHLY OUTPUT
    ! Write monthly value
    !-------------------------------------------------------------------------
    if (.not.spinup) then
      do moy=1,nmonth

        itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(moy-1)/real(nmonth)

      end do
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(year) + real(firstyeartrend) - real(spinupyears)

    write(301,999) itime, sum(outaClitt(:,jpngr))
    write(302,999) itime, sum(outaCsoil(:,jpngr))
    write(304,999) itime, sum(outanreq(:,jpngr))
    write(305,999) itime, sum(outaClit2soil(:,jpngr))
    write(306,999) itime, sum(outaNlit2soil(:,jpngr))
    write(313,999) itime, sum(outaCdsoil(:,jpngr))
    write(314,999) itime, sum(outaNdsoil(:,jpngr))
    write(315,999) itime, sum(outaNimmo(:,jpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_littersom

end module littersom_lpj