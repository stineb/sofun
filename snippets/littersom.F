c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE LITTERSOM
c     Litter and soil decomposition
c     Incorporates analytical solution for soil pool sizes once litter inputs
c     are (assumed to be) at equilibrium, reducing spin-up time for carbon
c     fluxes due to soil respiration.
c     ------------------------------------------------------------------

      subroutine littersom(year,jpngr,m,dm,initcell)
      
      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "soil.inc"
#include "flux.inc"
#include "pools.inc"
#include "pftpar.inc"
#include "climate.inc"
#include "output.inc"
#include "litterdecomp.inc"
#include "soildecomp.inc"
      include 'steer.inc'
      include 'landuse.inc'
      include 'coordinates.inc'
      include 'permafrost.inc'

c     ARGUMENTS
      INTEGER year,jpngr,m,dm
      LOGICAL initcell

c     LOCAL PARAMETERS
      INTEGER npools
c     1=litter_leaf, 2=litter_woody, 3=exu&litter_root, 4=soil_fast, 5=soil_slow
      PARAMETER(npools=5)
      INTEGER MIDX              !mixed layer index
      PARAMETER(MIDX=airl+snowl+mixedl)
#if _sublayer_hydhd
      INTEGER CIDX              !first catotelm layer index
      PARAMETER(CIDX=airl+snowl+mixedl+nacrotelm+1)
#endif
c     LOCAL VARIABLES 
      INTEGER SLIDX         ! loops
      INTEGER pft,ppft,lu,n,i
      INTEGER kindex,lindex     !time step counter (either months or days)
      REAL*8 respir_modifier(npools)      
      REAL*8 k_litter_leaf(npft) !monthly leaf litter decomposition rate (/month)
      REAL*8 k_litter_woody(npft) !monthly woody litter decomposition rate (/month)
      REAL*8 k_litter_root(npft) !monthly root litter decomposition rate (/month)
      REAL*8 k_exu(npft)        !monthly exudates decomposition rate (/month)

      REAL*8 soilwater(nlu,nsl+mixedl) !auxiliary variable
      REAL*8 soiltemp(nlu,nsl+mixedl)  !auxiliary variable
      REAL*8 thawdepth(nlu)     !auxiliary variable
      REAL*8 watertpos          !auxiliary variable
      REAL*8 whc
      REAL*8 totsoilC           !sum of slow and fast C pools

      REAL*8 cflux_litter_soil  !litter decomposition flux to soil
      REAL*8 cflux_litter_atmos !litter decomposition flux to atmosphere
      REAL*8 cflux_litter_soil_pft !auxiliary variable in PFT loop
      REAL*8 cflux_litter_atmos_pft !auxiliary variable in PFT loop
      REAL*8 cflux_fast_atmos(nlu)   !soil fast pool decomposition flux to atmosphere
      REAL*8 cflux_slow_atmos   !soil slow pool decomposition flux to atmosphere

      REAL*8 litterdag_fast     !above-ground component of fast litter decomposition
      REAL*8 litterdag_slow     !above-ground component of slow litter decomposition
      REAL*8 litter_decom_bg    !below-ground component of litter decomposition
      REAL*8 exu_decom(ncvar)   !exudates component of litter decomposition
      REAL*8 tmp_ntoc_fast      ! local C/N ratio
      REAL*8 tmp_ntoc_slow      ! local C/N ratio


c     NEW LOCAL VARIABLES/PARAMETERS FOR PEATLANDS
c     AND MONTHLY LITTER UPDATE
      
#if _peatlands_on
c     Peatland Parameters
      REAL*8 acro_por           ! acrotelm porosity minus gas fraction
      PARAMETER (acro_por=acrotelm_por-Fgas)
c      REAL*8 cpool_acro         ! carbon pool of norm acrotelm
c      PARAMETER (cpool_acro=7500.0d0) ! 300mm, dens=25kg C/m3=25g C/m2 /mm
c     New density: Belyea and Climo, 2001
c      PARAMETER (cpool_acro=5616.0d0) ! 300mm, dens=18.72kg C/m3=18.72g C/m2 /mm 
#endif

c     estimate equilibration time for SPINUP with DyN before 
c     analytical solution of carbon pools is applied
c     testing shows that a typical time scale is around 50 years
c     must be at least 1 recycle period for averaging
      INTEGER N_equil_time      
      PARAMETER(N_equil_time = max(recycle+1,48)) ! years
#if _DyN
      REAL*8 nh4_inc
      REAL*8 ntoc_litter_pft
      REAL*8 nsource
#if _peatlands_on
      REAL*8 litter_N_diag      
      REAL*8 cton_soil_peat      
#endif
#endif

      SAVE whc,watertpos,cflux_fast_atmos

#if _litter_monthly
c     Add a monthly (daily) fraction of litter that was produced in last year's
c     annually executed SRs (kill, ... , light). Note that in GPP, litter depth
c     is calculated and used by heat diffusion. Calling update_litter SR here
c     causes the first month to lack this litter addition, however has the
c     advantage of being able to use when _sublayer_hydhd==1.
c     -------------------------------------------------------------------------
#if _fixedlandmask==0
      if (.not.newmill) call update_litter(jpngr,m,dm)
#else
      call update_litter(jpngr,m,dm)
#endif
#endif

c     -------------------------------------------------------------------------
c     LU LOOP
c     -------------------------------------------------------------------------


      if (m.eq.1.and.dm.eq.1) then
c     /////////////////////////////////////////////////////////////////////////
c     INITIALISATION IN DAY OF THE YEAR
c     -------------------------------------------------------------------------
c     These variables are saved or common, thus available at next call to SR.
c     -------------------------------------------------------------------------
        cflux_fast_atmos(:) = 0.0d0
        litter_decom(:,:,:) = 0.0d0
        litter_decom_pft(:,:,:) = 0.0d0
        
        arh(:,:) = 0.0d0
        mrh(:,:,:) = 0.0d0
#if _output_monthly
        outrh(jpngr,:,:,:) = 0.0d0
#if _DyN
        mlitterc(jpngr,:,:) = 0.0d0
        mlittern(jpngr,:,:) = 0.0d0
        nmob(jpngr,:,:)     = 0.0d0
#endif          
#endif
        k_slow_out(jpngr,:) = 0.0d0
        input_slow_out(jpngr,:) = 0.0d0
#if _peatlands_on
        acro_to_cato_ann(jpngr) = 0.0d0
#endif
c     Water holding capacity initialisation
        whc=soilpar(3,jpngr)/d1
c     Water table position initialisation
        watertpos = 0.0d0
#if _permafrostc_on
c     ----------------------------------------------------------------
c     PERMAFROST C TRANSFER
c     ----------------------------------------------------------------
c     In LPX permafrost occures in a grid cell when a soil layer (1-8)
c     is completely frozen for more than 2 years. The depth to this
c     layer is called thawdepth, which is calculated daily.
c     If the soil freezes from top down, but not completely down to 2m,
c     the depth is called freezing depth.
c     Thawdepth then decides wheter carbon from the fast and the slow
c     soil pools is moved to the permafrost carbon pools or vice versa.
c
c     Assumptions based on soil carbon profiles of ???:
c     1) Surface density in grid cell with permafrost is 60 kgC/m³
c     2) Permafrost density in grid cell with permafrost is 21 kgC/m³
c     3) Permafrost density is constant with depth
c     4) Carbon density in thaw layer (0cm to thawdepth) is decreasing
c        exponentially from 60 to 21 kgC/m³
c     5) Carbon exceeding the total possible amount in thaw layer is
c        moved to permafrost carbon pool
c     6) Carbon from the both soil carbon pools are moved in
c        proportional amounts to permafrost carbon pools
c     7) If thawdepth increases, carbon from the permafrost carbon
c        pools are moved back to the soil carbon pools
c     8) If thawdepth increases to 2m depth, permafrost carbon is
c        transfered quickly back to soil as a constant fraction (10%/yr)

        soil_to_permaf_fast(:,:) = 0.0d0
        soil_to_permaf_slow(:,:) = 0.0d0

        do lu=1,nlu
          if (permafrostyr(jpngr,lu).ge.2.and.
     $         athawdepth_old(jpngr,lu,1).lt.d1+d2) then

            if (lu.ne.lupeat.and.lu.ne.lupeatold) then

              if (athawdepth_old(jpngr,lu,1).le.athawdepth_old(jpngr,lu,2)) then
                if (cpool_fast(jpngr,lu,1).gt.0.0d0.or.
     $               cpool_slow(jpngr,lu,1).gt.0.0d0) then
                  call build_permafrost_carbon(jpngr,lu,athawdepth_old(jpngr,lu,1),
     $                 soil_to_permaf_fast(lu,:),soil_to_permaf_slow(lu,:))
                endif
              else
                if (cpool_permaf_fast(jpngr,lu,1).gt.0.0d0.or.
     $               cpool_permaf_slow(jpngr,lu,1).gt.0.0d0) then
                  call thaw_permafrost_carbon(jpngr,lu,athawdepth_old(jpngr,lu,1),
     $                 athawdepth_old(jpngr,lu,1),soil_to_permaf_fast(lu,:),
     $                 soil_to_permaf_slow(lu,:))
                endif
              endif

            endif ! not lupeat or lupeatold
            
          else
c     Case where permafrost is gone (thawdepth greater equal 2m)
c     -------------------------------------------------------------------------

c     Check if permafrost C is left to be transferred to fast soil pool
            if (cpool_permaf_fast(jpngr,lu,1).gt.0.0d0) then
              soil_to_permaf_fast(lu,1)=-k_thaw*cpool_permaf_fast(jpngr,lu,1)
     $             /dble(ntimesteps)
              soil_to_permaf_fast(lu,2:ncvar)=cpool_permaf_fast(jpngr,lu,2:ncvar)
#if _DyN
              soil_to_permaf_fast(lu,n_id)=-k_thaw*cpool_permaf_fast(jpngr,lu,n_id)
     $             /dble(ntimesteps)
#endif
            endif
            
c     Check if permafrost C is left to be transferred to slow soil pool
            if (cpool_permaf_slow(jpngr,lu,1).gt.0.0d0) then
              soil_to_permaf_slow(lu,1)=-k_thaw*cpool_permaf_slow(jpngr,lu,1)
     $             /dble(ntimesteps)
              soil_to_permaf_slow(lu,2:ncvar)=cpool_permaf_slow(jpngr,lu,2:ncvar)
#if _DyN
              soil_to_permaf_slow(lu,n_id)=-k_thaw*cpool_permaf_slow(jpngr,lu,n_id)
     $             /dble(ntimesteps)
#endif
            endif
            
          endif
          
        enddo
#endif
      endif                     ! if (m.eq.1.and.dm.eq.1)


c     /////////////////////////////////////////////////////////////////////////
c     DO EVERY DAY (MONTH)
c     -------------------------------------------------------------------------

#if (_output_monthly==1) && (_DyN==1)
c     Record monthly litter C and N (sum of fast, slow and below-ground) at the
c     beginnig of the month. 
c     -------------------------------------------------------------------------
      if (dm.eq.1) then
        do lu=1,nlu      
          do pft=1,npft
            if (lu_category(pft).eq.lu) then
              mlitterc(jpngr,m,lu) = mlitterc(jpngr,m,lu)
     $             +litter_ag_fast(pft,jpngr,1)
     $             +litter_ag_slow(pft,jpngr,1)
     $             +litter_bg(pft,jpngr,1)
              mlittern(jpngr,m,lu) = mlittern(jpngr,m,lu)
     $             +litter_ag_fast(pft,jpngr,n_id)
     $             +litter_ag_slow(pft,jpngr,n_id)
     $             +litter_bg(pft,jpngr,n_id)
            endif
          enddo                 !pft
        enddo                   !lu
      endif                     !dm=1
#endif
      
      
#if _dailymode
c     Use daily updated values for soil temperature and soil moisture.
c     -------------------------------------------------------------------------
      kindex=dm
      if (m.gt.1) then
        do i=1,m-1
          kindex=kindex+ndaymonth(i)
        enddo
      endif

      do i=mixedl+1,nlayers-airl-snowl
#if _sublayer_hydhd
#if _soil25cm_resp
c     Set always layer at 25cm depth for respiration modifiers
        SLIDX=airl+snowl+mixedl+3
#else
c     Set individual layers for respiration modifiers
        SLIDX=i+airl+snowl
#endif
        do lu=1,nlu
          soiltemp(lu,i) = Tsoil(kindex,SLIDX,lu)
          if (lu.ne.lupeat) then
            soilwater(lu,i) = max((Fwater(kindex,SLIDX,lu)
     $           -Fpwp(SLIDX,jpngr))/whc,0.0d0)
#if _peatlands_on
          else
            if (SLIDX.lt.CIDX) then
c     select 3 acrotelm layers in peatlands, use different porosity
              soilwater(lupeat,i)=max(Fwater(kindex,SLIDX,lupeat)
     $             -Fpwp(SLIDX,jpngr),0.0d0)
     $             /(acro_por-Fpwp(SLIDX,jpngr))
            else
              soilwater(lupeat,i)=1.0d0 ! catotelm is always saturated
            endif
#endif
          endif
        enddo
#else
        soiltemp(:,i) = dtemp_soil(kindex,:)
        soilwater(:,i) = dw1(kindex,:)
#endif
      enddo
#if _sublayer_hydhd
#if _soil25cm_resp
      soiltemp(:,1) = Tsoil(kindex,SLIDX,:)
#else
      soiltemp(:,1) = Tsoil(kindex,MIDX,:)
#endif
      thawdepth(:) = dthawdepth(kindex,:)
#else
      soiltemp(:,1) = dtemp_soil(kindex,:)
      thawdepth(:) = d1+d2
#endif
      soilwater(:,1) = soilwater(:,2)

#if _peatlands_on
      if (lupeat.gt.0.and.ptgridcell(jpngr)) then
        watertpos = wtp(kindex)
        soilwater(lupeat,1) = 1.0d0
      endif
#endif        
#else
c     Use monthly means for soil temperature and soil moisture.
c     -------------------------------------------------------------------------
      if (dm.eq.ndaymonth(m)) then
        kindex = m

        do i=mixedl+1,nlayers-airl-snowl
#if _sublayer_hydhd
#if _soil25cm_resp
c     Set always layer at 25cm depth for respiration modifiers
          SLIDX=airl+snowl+mixedl+3
#else
c     Set individual layers for respiration modifiers
          SLIDX=i+airl+snowl
#endif
          do lu=1,nlu
            soiltemp(lu,i) = mTsoil(kindex,SLIDX,lu)
            if (lu.ne.lupeat) then
              soilwater(lu,i) = max((mFwater1(kindex,lu,jpngr)
     $             -Fpwp(SLIDX,jpngr))/whc,0.0d0)
#if _peatlands_on
            else
              if (SLIDX.lt.CIDX) then
c     select 3 acrotelm layers in peatlands, use different porosity
                soilwater(lupeat,i)=max(mFwater_i(kindex,SLIDX,jpngr)
     $               -Fpwp(SLIDX,jpngr),0.0d0)
     $               /(acro_por-Fpwp(SLIDX,jpngr))
              else
                soilwater(lupeat,i)=1.0d0 ! catotelm is always saturated
              endif
#endif
            endif
          enddo
#else
          soiltemp(:,i) = mtemp_soil(jpngr,kindex,:)
          soilwater(:,i) = mw1(kindex,:,jpngr)
#endif
        enddo
#if _sublayer_hydhd
#if _soil25cm_resp
        soiltemp(:,1) = mTsoil(kindex,SLIDX,:)
#else
        soiltemp(:,1) = mTsoil(kindex,MIDX,:)
#endif
        thawdepth(:) = mthawdepth(kindex,:)
#else
        soiltemp(:,1) = mtemp_soil(jpngr,kindex,:)
        thawdepth(:) = d1+d2
#endif
        soilwater(:,1) = soilwater(:,2)

#if _peatlands_on
        if (lupeat.gt.0.and.ptgridcell(jpngr)) then
          watertpos=mwtp(kindex,jpngr)
          soilwater(lupeat,1) = 1.0d0
        endif
#endif        
#endif

c     Initialise modifier
        respir_modifier(:)=0.0d0

        do lu=1,nlu

c     /////////////////////////////////////////////////////////////////////////
c     DECAY RATES
c     -------------------------------------------------------------------------
c     Calculate daily (monthly) decomposition rates as a function of
c     temperature and moisture
          
c     k = k_10 * respir_modifier

c     (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
c     from (1),
c     (2) c = c0*exp(-kt) where c0=initial pool size
c     from (2), decomposition in any month given by
c     (3) delta_c = c0 - c0*exp(-k)
c     from (4)
c     (4) delta_c = c0*(1.0-exp(-k))
c     -------------------------------------------------------------------------

!!!      call set_C_balance(jpngr,testveg,testlitter,testsoil,testflux)          

          do pft=1,npft
            if (lu.eq.lu_category(pft)) then
              
c     LITTER TEMPERATURE AND MOISTURE MODIFIER
c     -------------------------------------------------------------------------
c     Temperature: Lloyd & Taylor 1994, Brovkin et al., 2012
c     Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
c               Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012

c     leaf litter
              call litter_modifier(.false.,lu,pft,soiltemp(lu,1),
     $             soilwater(lu,1),respir_modifier(1))
c     woody litter
              call litter_modifier(.true.,lu,pft,soiltemp(lu,1),
     $             soilwater(lu,1),respir_modifier(2))
c     root litter & exudates
              call litter_modifier(.false.,lu,pft,soiltemp(lu,2),
     $             soilwater(lu,2),respir_modifier(3))

c     define pft dependent decomposition rates
#if _sitch_hydrology || _soil25cm_resp
              k_litter_leaf(pft)  = k_litter_fast10*respir_modifier(1)
     $             /dble(ntimesteps)
              k_litter_woody(pft) = k_litter_slow10*respir_modifier(2)
     $             /dble(ntimesteps)
              k_litter_root(pft)  = k_litter_fast10*respir_modifier(3)
     $             /dble(ntimesteps)
#else
              k_litter_leaf(pft)  = k_litter_leaf10(pft)*respir_modifier(1)
     $             /dble(ntimesteps)
              k_litter_woody(pft) = k_litter_woody10(pft)*respir_modifier(2)
     $             /dble(ntimesteps)
              k_litter_root(pft)  = k_litter_root10(pft)*respir_modifier(3)
     $             /dble(ntimesteps)
#endif              
              k_exu(pft)          = k_exu10*respir_modifier(3)
     $             /dble(ntimesteps)

            endif
          enddo

c     //////////////////////////////////////////////////////////////////////////
c     SOIL TEMPERATURE AND MOISTURE MODIFIER
c
c     Temperature: Lloyd & Taylor 1994
c     Moisture: Foley, 1995; Fang and Moncrieff, 1999; Gerten et al., 2004;
c               Wania et al., 2009; Frolking et al., 2010; Spahni et al., 2012
#if _sublayer_hydhd
          totsoilC = cpool_fast(jpngr,lu,1)+cpool_slow(jpngr,lu,1)

c     fast soil pool and acrotelm
          call soil_modifier_layers(jpngr,lu,1,soiltemp(lu,2:nsl+mixedl),
     $         soilwater(lu,2:nsl+mixedl),watertpos,totsoilC,
     $         thawdepth(lu),respir_modifier(4))

c     slow soil pool and catotelm
          call soil_modifier_layers(jpngr,lu,2,soiltemp(lu,2:nsl+mixedl),
     $         soilwater(lu,2:nsl+mixedl),watertpos,totsoilC,
     $         thawdepth(lu),respir_modifier(5))
#else
c     fast soil pool
          call soil_modifier(soiltemp(lu,1),soilwater(lu,1),respir_modifier(4))

c     slow soil pool
          call soil_modifier(soiltemp(lu,1),soilwater(lu,1),respir_modifier(5))

#endif
c     copy decomposition rate at timestep kindex
          k_fast(kindex,lu) = k_soil_fast10(lu)*respir_modifier(4)
     $         /dble(ntimesteps)
          k_slow(kindex,lu) = k_soil_slow10(lu)*respir_modifier(5)
     $         /dble(ntimesteps)

c     Record mean decomposition rate of slow soil pool for output
c     ----------------------------------------------------------------
          k_slow_out(jpngr,lu) = k_slow_out(jpngr,lu) + k_slow(kindex,lu)

          do pft=1,npft
            if (lu.eq.lu_category(pft)) then

c     ----------------------------------------------------------------
c     Avoid numerical instability implying exploding CN ratios.
c     This is not a very nice solution. Don't understand why
c     litter-C is not reduced any further once it's smaller than 1.d-16,
c     while litter-N still is still reduced to numerical limit (~1.d-305).
c     ----------------------------------------------------------------
              if (present(pft,jpngr)) everpresent(pft,jpngr) = .true.

              if (everpresent(pft,jpngr).and.(.not.present(pft,jpngr))
     $             .and.(litter_ag_fast(pft,jpngr,1).lt.1.0d-16)) then
                litter_ag_fast(pft,jpngr,:) = 0.0d0
              endif
              if (everpresent(pft,jpngr).and.(.not.present(pft,jpngr))
     $             .and.(litter_ag_slow(pft,jpngr,1).lt.1.0d-16)) then
                litter_ag_slow(pft,jpngr,:) = 0.0d0
              endif
              if (everpresent(pft,jpngr).and.(.not.present(pft,jpngr))
     $             .and.(litter_bg(pft,jpngr,1).lt.1.0d-16)) then
                litter_bg(pft,jpngr,:) = 0.0d0
              endif

              
c     ////////////////////////////////////////////////////////////////
c     EXUDATES DECAY
c     ----------------------------------------------------------------
c     Calculate the exudates respiration before litter respiration.
c     Exudates are mostly short organic compounds (poly- and mono-
c     saccharides, amino acids, organic acids, phenolic compounds and
c     enzymes) and are quickly respired and released as CO2.
c     ----------------------------------------------------------------                
              exu_decom(1) = exu_pool(pft,jpngr,1)*
     $             (1.0d0-exp(-k_exu(pft)))
              exu_decom(2:ncvar) = exu_pool(pft,jpngr,2:ncvar)

              call addc(jpngr,
     $             mrh(lu_category(pft),m,:),
     $             exu_decom(1),
     $             exu_decom(2:ncvar)
     $             )
              exu_pool(pft,jpngr,1) = exu_pool(pft,jpngr,1)-exu_decom(1)

c     ////////////////////////////////////////////////////////////////
c     LITTER DECAY
c     ----------------------------------------------------------------
              litterdag_fast=litter_ag_fast(pft,jpngr,1)
     $             * (1.0d0-exp(-k_litter_leaf(pft))) 
              litterdag_slow=litter_ag_slow(pft,jpngr,1)
     $             * (1.0d0-exp(-k_litter_woody(pft))) 
              litter_decom_bg=litter_bg(pft,jpngr,1)
     $             * (1.0d0-exp(-k_litter_root(pft)))
              
              call addc(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litterdag_fast,
     $             litter_ag_fast(pft,jpngr,2:ncvar)
     $             )
              call addc(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litterdag_slow,
     $             litter_ag_slow(pft,jpngr,2:ncvar)
     $             )
              call addc(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litter_decom_bg,
     $             litter_bg(pft,jpngr,2:ncvar)
     $             )
              
c     Update the litter pools
              litter_ag_fast(pft,jpngr,1)=litter_ag_fast(pft,jpngr,1)
     $             -litterdag_fast
              litter_ag_slow(pft,jpngr,1)=litter_ag_slow(pft,jpngr,1)
     $             -litterdag_slow
              litter_bg(pft,jpngr,1)=litter_bg(pft,jpngr,1)
     $             -litter_decom_bg
              
              call addc(jpngr,
     $             litter_decom(kindex,lu_category(pft),:),
     $             litter_decom_pft(kindex,pft,1),
     $             litter_decom_pft(kindex,pft,2:ncvar)
     $             )

#if _DyN
c     ----------------------------------------------------------------
c     Mineralization of litter-N goes with the same rate as C (somdynam.cpp:326)
c     Use the same variables (litterdag_*) also for N (unused further below)
c     ----------------------------------------------------------------
              litterdag_fast = litter_ag_fast(pft,jpngr,n_id)
     $             * (1.0d0-exp(-k_litter_leaf(pft))) !eqn 4
              litterdag_slow = litter_ag_slow(pft,jpngr,n_id)
     $             * (1.0d0-exp(-k_litter_woody(pft))) !eqn 4
              litter_decom_bg = litter_bg(pft,jpngr,n_id)
     $             * (1.0d0-exp(-k_litter_root(pft)))
              
              call addn(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litterdag_fast
     $             )
              call addn(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litterdag_slow
     $             )
              call addn(jpngr,
     $             litter_decom_pft(kindex,pft,:),
     $             litter_decom_bg
     $             )
              
              litter_ag_fast(pft,jpngr,n_id) =
     $             litter_ag_fast(pft,jpngr,n_id)-litterdag_fast
              litter_ag_slow(pft,jpngr,n_id) =
     $             litter_ag_slow(pft,jpngr,n_id)-litterdag_slow
              litter_bg(pft,jpngr,n_id) =
     $             litter_bg(pft,jpngr,n_id)-litter_decom_bg
              
              call addn(jpngr,
     $             litter_decom(kindex,lu_category(pft),:),
     $             litter_decom_pft(kindex,pft,n_id)
     $             )

c     ----------------------------------------------------------------
c     Mineralized N of litter is converted to NH4+
c     ----------------------------------------------------------------
c     SPINUP: speed up N-cycle by prescribing C related N input to NH4+

              if (spinup.and.year.le.soil_equil_year-N_equil_time) then

                nh4_inc = atmosfrac(lu)*litter_decom_pft(kindex,pft,1)
     $               /cton_soil(pft)
                call addn1(jpngr,
     $               nh4(jpngr,lu,:),
     $               nh4_inc
     $               )
                call addn1(jpngr,
     $               nmob(jpngr,lu,:),
     $               nh4_inc
     $               )
                
              else

c     NON-SPINUP: mineralize actual N litter pool to increase NH4+
                nh4_inc = atmosfrac(lu)*litter_decom_pft(kindex,pft,n_id)
                call addn1(jpngr,
     $               nh4(jpngr,lu,:),
     $               nh4_inc
     $               )
                call addn1(jpngr,
     $               nmob(jpngr,lu,:),
     $               nh4_inc
     $               )
              endif
c     ----------------------------------------------------------------
c     N FLUX LITTER->SOIL
c     ----------------------------------------------------------------
c     Use fixed soil C/N ratio (Tab.4, XP08) for flux per PFT. This
c     should yield the same results as using pft-specific soil pools
c     and calculating the soil nh4 production from SOM mineralization
c     divided by the pft-specific C/N ratio (see somdynam.cpp:384)
c     The present calculation implies that soil-N pools have to be
c     simulated explicitly.
c     ----------------------------------------------------------------
c     Forced soil N-input = soilfrac*litter_decom_pft(kindex,pft,1)/cton_soil(pft)
c     Explicit soil N-input = soilfrac*litter_decom_pft(kindex,pft,n_id)
c     => Use greater source of each
c     ----------------------------------------------------------------
#if _fixnsrc
              if (spinup.and.litter_decom_pft(kindex,pft,1)/cton_soil(pft).gt.
     $             litter_decom_pft(kindex,pft,n_id)) then
#else
              if (litter_decom_pft(kindex,pft,1)/cton_soil(pft).gt.
     $             litter_decom_pft(kindex,pft,n_id)) then
#endif
                cflux_litter_soil_pft =
     $               soilfrac(lu)*litter_decom_pft(kindex,pft,1)              

                call addn(jpngr,
     $               cpool_fast(jpngr,lu,:),
     $               cflux_litter_soil_pft/cton_soil(pft)
     $               *fastfrac
     $               )
                call addn(jpngr,
     $               cpool_slow(jpngr,lu,:),
     $               cflux_litter_soil_pft/cton_soil(pft)
     $               *(1.d0-fastfrac)
     $               )

c     ----------------------------------------------------------------
c     EVALUATE IMPLICIT N SOURCE
c     ----------------------------------------------------------------
c     Implicit N source is a result of holding the soil CN-ratio
c     constant (and calculating soil-N = soil-C/CN-ratio).
c     Only this amount of N "should" enter soil from litter decomposition:
c     litter_decom_pft(kindex,pft,n_id)*soilfrac
c     But, in fact, this amount "implicitly" enters soil:
c     cflux_litter_soil_pft/cton_soil(pft)
c     The difference is the implicit N source:
c     ----------------------------------------------------------------              
                if (lu.eq.lubuilt.or.lu.eq.luburied) then
                  nsource = 0.0d0
                else
                  nsource = cflux_litter_soil_pft/cton_soil(pft)
     $                 - litter_decom_pft(kindex,pft,n_id)*soilfrac(lu)
                endif

                anflux_nsource_tot(jpngr,lu,1) =
     $               anflux_nsource_tot(jpngr,lu,1)
     $               + nsource

c     ----------------------------------------------------------------

              else
c     ----------------------------------------------------------------              
c     In this case (explicit soil N-input is large enough),
c     the implicit N source is zero!
c     ----------------------------------------------------------------              
                
                call addn(jpngr,
     $               cpool_fast(jpngr,lu,:),
     $               litter_decom_pft(kindex,pft,n_id)*soilfrac(lu)
     $               *fastfrac
     $               )
                call addn(jpngr,
     $               cpool_slow(jpngr,lu,:),
     $               litter_decom_pft(kindex,pft,n_id)*soilfrac(lu)
     $               *(1.d0-fastfrac)
     $               )

              endif
#endif
              
            endif               !(lu_category(pft).eq.lu)
          enddo                 !pft

#if _DyN && _fixnsrc          
          if (.not.spinup.and.lu.eq.lunat) then

            nsource = nsrc_buf(mod(year,climbuf)+1,jpngr)/dble(ntimesteps)

            anflux_nsource_tot(jpngr,lu,1) = anflux_nsource_tot(jpngr,lu,1)
     $           + nsource
            
            call addn(jpngr,
     $           cpool_fast(jpngr,lu,:),
     $           nsource*fastfrac
     $           )
            call addn(jpngr,
     $           cpool_slow(jpngr,lu,:),
     $           nsource*(1.d0-fastfrac)
     $         )

          endif
#endif

c     ----------------------------------------------------------------
c     C FLUX LITTER -> SOIL/ATMOSPHERE
c     ----------------------------------------------------------------
c     Calculate carbon flux to atmosphere and soil
          cflux_litter_atmos=atmosfrac(lu)*litter_decom(kindex,lu,1) 
          cflux_litter_soil=soilfrac(lu)*litter_decom(kindex,lu,1)

c     Further subdivide soil fraction between fast and slow soil pools
          call addc(jpngr,
     $         cpool_fast(jpngr,lu,:),
     $         cflux_litter_soil*fastfrac,
     $         litter_decom(kindex,lu,2:ncvar)
     $         )
          call addc(jpngr,
     $         cpool_slow(jpngr,lu,:),
     $         cflux_litter_soil*(1.d0-fastfrac),
     $         litter_decom(kindex,lu,2:ncvar)
     $         )

c     Record input to slow soil pool (catotelm) for output
          input_slow_out(jpngr,lu) = input_slow_out(jpngr,lu)
     $         + cflux_litter_soil*(1.d0-fastfrac)

#if _peatlands_on
c     ----------------------------------------------------------------
c     PEATLAND C TRANSFER
c     ----------------------------------------------------------------
c     In peatlands, C input happens at the top of
c     the peat column. The organic material gets compressed over time
c     as more material accumulates on the top. At some point, the
c     material reaches the boundary between the partially oxic,
c     partially anoxic layer (acrotelm) and the continuously anoxic
c     layer (catotelm).  Decomposition rates are much slower in the
c     catotelm then in the acrotelm because most of the quickly
c     decomposable material has been decomposed already and because of
c     the lack of oxygen.
c     Here, I assume that cpool_fast represents the acrotelm and
c     cpool_slow the catotelm.  To represent the process of C
c     transfer, a constant fraction of the cpool_fast is transferred
c     to the cpool_slow, once enough organic material has accumulated.
c     The value 7500 comes from assuming a bulk density of 0.05 g/cm3,
c     which is 0.025g C/cm3, which is 25kg C/m3 soil. The acrotelm is
c     0.3m thick, hence 25*0.3 is 7.5kg C is 7500g C/m2 in the
c     acrotelm.  Assume that a fraction of the acrotelm is transfered
c     to the catotelm through the growth of the acrotelm. See Clymo
c     1984. 
c     Following Belyea and Clymo, 2001, bulk density at acro-cato 
c     transition ranges from 0.011664 to 0.077841 with an average of
c     0.03744 gm/cm3 thus 0.01872 gC/cm3 = 5616 gC/m2
c     ----------------------------------------------------------------
          if (lu.eq.lupeat.and.ptgridcell(jpngr)) then

c     Constant transfer rate
c     ----------------------------------------------------------------
c     12 g C/m2/yr Clymo 1984, p.617 (tunable), test: 24 g C/m2/yr            
c          acro_to_cato(1) = acrotocato / dble(ntimesteps)

c     Variable transfer rate
c     ----------------------------------------------------------------
c            print*,'cpool_acro_old(jpngr)',cpool_acro_old(jpngr)
              acro_to_cato(1) = (cpool_acro_old(jpngr)-cpool_acro_norm(jpngr))
     $             /dble(ntimesteps)

#if _DyN
c     Calculate average C/N for litter flux into peatland soils based on fixed C/N
            litter_N_diag = 0.0d0
            do pft=1,npft
              if (lu_category(pft).eq.lupeat) then
                litter_N_diag = litter_N_diag
     $               +litter_decom_pft(kindex,pft,1)/cton_soil(pft)
              endif
            enddo
            if (litter_N_diag.gt.0.0d0) then
              cton_soil_peat = litter_decom(kindex,lupeat,1)/litter_N_diag
            else
              cton_soil_peat=0.0d0
              do pft = 1,npft
                if (lu_category(pft).eq.lupeat) then
                  cton_soil_peat = cton_soil_peat+cton_soil(pft)/dble(npftpeat)
                end if
              end do
            endif
#endif

            if (acro_to_cato(1).ge.0.0d0) then
c     Prevent negative fast carbon pool
              acro_to_cato(1)=min(acro_to_cato(1),cpool_fast(jpngr,lu,1))

#if _DyN
              if (cpool_fast(jpngr,lu,1).gt.0.0d0) then

                  acro_to_cato(n_id)=min(
     $                 acro_to_cato(1)
     $                 *cpool_fast(jpngr,lu,n_id)/cpool_fast(jpngr,lu,1),
     $                 cpool_fast(jpngr,lu,n_id))
c                endif

              else
                acro_to_cato(1) = 0.0d0
                acro_to_cato(n_id) = 0.0d0
              endif
#endif
#if compc13
            acro_to_cato(2:ncvar)=cpool_fast(jpngr,lu,2:ncvar)
#endif
            else
c     Prevent negative slow carbon pool
              acro_to_cato(1)=max(acro_to_cato(1),-cpool_slow(jpngr,lu,1))
#if _DyN
              if (cpool_slow(jpngr,lu,1).gt.0.0d0) then
                acro_to_cato(n_id)=max(
     $               acro_to_cato(1)
     $               *cpool_slow(jpngr,lu,n_id)/cpool_slow(jpngr,lu,1),
     $               -cpool_slow(jpngr,lu,n_id))
              else
                acro_to_cato(1) = 0.0d0
                acro_to_cato(n_id) = 0.0d0
             endif
#endif
#if compc13
              acro_to_cato(2:ncvar)=cpool_slow(jpngr,lu,2:ncvar)
#endif
            endif


            call addc(jpngr,
     $           cpool_fast(jpngr,lu,:),
     $           -acro_to_cato(1),
     $           acro_to_cato(2:ncvar)                 !fractionation?
     $           )

            if (spinup.and.year.le.soil_equil_year-N_equil_time) then

              call addc(jpngr,
     $             cpool_slow(jpngr,lu,:),
     $             (exp(1.5d0*(1.0d0-year/(soil_equil_year-N_equil_time)))
     $             -1.0d0+fastfrac)*cflux_litter_soil,
     $             acro_to_cato(2:ncvar)
     $             )
              
            else

              call addc(jpngr,
     $             cpool_slow(jpngr,lu,:),
     $             acro_to_cato(1),
     $             acro_to_cato(2:ncvar)
     $             )
              
            endif

c     Record input to slow soil pool (catothelm) for output
            input_slow_out(jpngr,lu) = input_slow_out(jpngr,lu)
     $           + acro_to_cato(1)

c     Record annual input to slow soil pool (catotelm) for equilibrium
            acro_to_cato_ann(jpngr) = acro_to_cato_ann(jpngr)+acro_to_cato(1)
            
#if _DyN            
            call addn(jpngr,
     $           cpool_fast(jpngr,lu,:),
     $           -acro_to_cato(n_id)
     $           )
            
            if (spinup.and.year.le.soil_equil_year-N_equil_time) then

              call addn(jpngr,
     $             cpool_slow(jpngr,lu,:),
     $             (exp(1.5d0*(1.0d0-year/(soil_equil_year-N_equil_time)))
     $             -1.0d0+fastfrac)*cflux_litter_soil/cton_soil_peat
     $             )
              
            else
              
              call addn(jpngr,
     $             cpool_slow(jpngr,lu,:),
     $             acro_to_cato(n_id)
     $             )

            endif
#endif            
          endif
#endif
#if _permafrostc_on
c     ----------------------------------------------------------------
c     PERMAFROST C TRANSFER (only for non-peatlands)
c     ----------------------------------------------------------------
c     In LPX permafrost occures in a grid cell when a soil layer (1-8)
c     is completely frozen for more than 2 years. The depth to this
c     layer is called thawdepth, which is calculated daily.
c     If the soil freezes from top down, but not completely down to 2m,
c     the depth is called freezing depth.
c     Thawdepth then decides wheter carbon from the fast and the slow
c     soil pools is moved to the permafrost carbon pools or vice versa.

c     Carbon is transferred from soil to permafrost or vice versa

          if (soil_to_permaf_fast(lu,1).ge.0.0d0) then
c     Prevent negative fast soil carbon pool
            soil_to_permaf_fast(lu,1)=min(soil_to_permaf_fast(lu,1),
     $           cpool_fast(jpngr,lu,1))
#if _DyN
            soil_to_permaf_fast(lu,n_id)=min(soil_to_permaf_fast(lu,n_id),
     $           cpool_fast(jpngr,lu,n_id))
#endif
          else
c     Prevent negative fast permafrost carbon pool
            soil_to_permaf_fast(lu,1)=max(soil_to_permaf_fast(lu,1),
     $           -cpool_permaf_fast(jpngr,lu,1))
#if _DyN
            soil_to_permaf_fast(lu,n_id)=max(soil_to_permaf_fast(lu,n_id),
     $           -cpool_permaf_fast(jpngr,lu,n_id))
#endif        
          endif
          
          if (soil_to_permaf_slow(lu,1).ge.0.0d0) then
c     Prevent negative slow soil carbon pool
            soil_to_permaf_slow(lu,1)=min(soil_to_permaf_slow(lu,1),
     $           cpool_slow(jpngr,lu,1))
#if _DyN
            soil_to_permaf_slow(lu,n_id)=min(soil_to_permaf_slow(lu,n_id),
     $           cpool_slow(jpngr,lu,n_id))
#endif
          else
c     Prevent negative slow permafrost carbon pool
            soil_to_permaf_slow(lu,1)=max(soil_to_permaf_slow(lu,1),
     $           -cpool_permaf_slow(jpngr,lu,1))
#if _DyN
            soil_to_permaf_slow(lu,n_id)=max(soil_to_permaf_slow(lu,n_id),
     $           -cpool_permaf_slow(jpngr,lu,n_id))
#endif        
          endif
          
          
c     Transfer carbon

          call addc(jpngr,
     $         cpool_permaf_fast(jpngr,lu,:),
     $         soil_to_permaf_fast(lu,1),
     $         soil_to_permaf_fast(lu,2:ncvar)
     $         )
          
          call addc(jpngr,
     $         cpool_permaf_slow(jpngr,lu,:),
     $         soil_to_permaf_slow(lu,1),
     $         soil_to_permaf_slow(lu,2:ncvar)
     $         )
          
          call addc(jpngr,
     $         cpool_fast(jpngr,lu,:),
     $         -soil_to_permaf_fast(lu,1),
     $         soil_to_permaf_fast(lu,2:ncvar)
     $         )
          
          call addc(jpngr,
     $         cpool_slow(jpngr,lu,:),
     $         -soil_to_permaf_slow(lu,1),
     $         soil_to_permaf_slow(lu,2:ncvar)
     $         )
          
#if _DyN
          call addn(jpngr,
     $         cpool_permaf_fast(jpngr,lu,:),
     $         soil_to_permaf_fast(lu,n_id)
     $         )
          
          call addn(jpngr,
     $         cpool_permaf_slow(jpngr,lu,:),
     $         soil_to_permaf_slow(lu,n_id)
     $         )
          
          call addn(jpngr,
     $         cpool_fast(jpngr,lu,:),
     $         -soil_to_permaf_fast(lu,n_id)
     $         )
          
          call addn(jpngr,
     $         cpool_slow(jpngr,lu,:),
     $         -soil_to_permaf_slow(lu,n_id)
     $         )
#endif
#endif

          
c     ////////////////////////////////////////////////////////////////
c     SOIL C (AND N) DECAY
c     ----------------------------------------------------------------
c     Calculate daily/monthly soil decomposition to the atmosphere
c     1st: fluxes for non-lu
          cflux_fast_atmos(lu) =
     $         cpool_fast(jpngr,lu,1)*(1.0d0-exp(-k_fast(kindex,lu))) !eqn 4
          cflux_slow_atmos =
     $         cpool_slow(jpngr,lu,1)*(1.0d0-exp(-k_slow(kindex,lu))) !eqn 4
          
c     Update the soil pools
c     1st: updates for non-lu
c     ----------------------------------------------------------------
          cpool_fast(jpngr,lu,1)=cpool_fast(jpngr,lu,1)-cflux_fast_atmos(lu)
          cpool_slow(jpngr,lu,1)=cpool_slow(jpngr,lu,1)-cflux_slow_atmos

#if _DyN
c     Record monthly (daily) soil turnover flux (labile carbon availability)
c     ----------------------------------------------------------------
          doc(lu) = cflux_fast_atmos(lu) + cflux_slow_atmos

c     ----------------------------------------------------------------
c     Turnover rate of soil N pool is equal to tor of soil C pool.
c     The following statement is analogous to somdynam.cpp:384. See
c     comment above.
c     ----------------------------------------------------------------
          nh4_inc = cpool_fast(jpngr,lu,n_id)*(1.0d0-exp(-k_fast(kindex,lu)))
          cpool_fast(jpngr,lu,n_id) = cpool_fast(jpngr,lu,n_id)-nh4_inc
          call addn1(jpngr,
     $         nh4(jpngr,lu,:),
     $         nh4_inc
     $         )

          nh4_inc = cpool_slow(jpngr,lu,n_id)*(1.0d0-exp(-k_slow(kindex,lu)))
          cpool_slow(jpngr,lu,n_id) = cpool_slow(jpngr,lu,n_id)-nh4_inc
          call addn1(jpngr,
     $         nh4(jpngr,lu,:),
     $         nh4_inc
     $         )

#endif

c     ////////////////////////////////////////////////////////////////
c     MONTHLY AND ANNUAL HETEROTROPHIC RESPIRATION
c     ----------------------------------------------------------------
          call addc(jpngr,
     $         mrh(lu,m,:),
     $         cflux_litter_atmos,
     $         litter_decom(kindex,lu,2:ncvar)
     $         )
          call addc(jpngr,
     $         mrh(lu,m,:),
     $         cflux_fast_atmos(lu),
     $         cpool_fast(jpngr,lu,2:ncvar)
     $         )
          call addc(jpngr,
     $         mrh(lu,m,:),
     $         cflux_slow_atmos,
     $         cpool_slow(jpngr,lu,2:ncvar)
     $         )

          if (dm.eq.ndaymonth(m)) then
            call addc(jpngr,
     $           arh(lu,:),
     $           mrh(lu,m,1),
     $           mrh(lu,m,2:ncvar)
     $           )
#if _output_monthly
            outrh(jpngr,m,lu,:)=mrh(lu,m,:)
#endif
          endif
          
c     do not empty soil pools below a minimum threshold (in contrast to sitch version!)
        enddo                   !lu                  

#if (_dailymode==0)
      endif                     ! (dm.eq.ndaymonth(m))
#endif


c     /////////////////////////////////////////////////////////////////////////
c     DO AT THE END OF THE YEAR
c     --------------------------------------------------------------------------
      if (m.eq.nmonth.and.dm.eq.ndaymonth(m)) then

c     Save litterpools to know how much was added to litter in SRs that are exe-
c     cuted annually (calculated in update_litter SR).
c     --------------------------------------------------------------------------
        litter_ag_fast_save(:,jpngr,1) = litter_ag_fast(:,jpngr,1)
        litter_ag_slow_save(:,jpngr,1) = litter_ag_slow(:,jpngr,1)
        litter_bg_save(:,jpngr,1) = litter_bg(:,jpngr,1)
#if _DyN
        litter_ag_fast_save(:,jpngr,2) = litter_ag_fast(:,jpngr,n_id)
        litter_ag_slow_save(:,jpngr,2) = litter_ag_slow(:,jpngr,n_id)
        litter_bg_save(:,jpngr,2) = litter_bg(:,jpngr,n_id)
#endif

c     /////////////////////////////////////////////////////////////////////////
c     ANALYTICAL CALCULATION OF SOIL C and N EQUILIBRIUM
c     --------------------------------------------------------------------------
c     Analytical solution of differential flux equations for fast and slow
c     soil carbon pools.  Implemented after (soil_equil_year) simulation
c     years, when annual litter inputs should be close to equilibrium.  Assumes
c     average climate (temperature and soil moisture) from all years up to
c     soil_equil_year.

c     Fill up permafrost soils with 25% of maximum carbon content in thaw layer.
c     This will help spinup to speed up considerably for permafrost soils.
c     Also fill up acrotelm soil carbon pool.
        if (spinup.and.year.eq.soil_initial_year) then
          do lu=1,nlu

            if (lu.eq.lupeat) then
#if _peatlands_on
              cpool_fast(jpngr,lu,1) = cpool_acro_norm(jpngr)
#if _DyN
              cpool_fast(jpngr,lu,1) = cpool_fast(jpngr,lu,1)/cton_soil_peat
#endif
#endif
            else
#if _permafrostc_on
              if (permafrostyr(jpngr,lu).ge.2.and.
     $             athawdepth_old(jpngr,lu,1).lt.d1+d2) then
#if _DyN
c     Store N/C ratio before soil equilibration to impose same ratio after equilibration.
                if (cpool_fast(jpngr,lu,1).gt.0.d0) then
                  tmp_ntoc_fast = cpool_fast(jpngr,lu,n_id)
     $                 /cpool_fast(jpngr,lu,1)
                else
                  tmp_ntoc_fast = 0.0d0
                endif
                if (cpool_slow(jpngr,lu,1).gt.0.d0) then
                  tmp_ntoc_slow = cpool_slow(jpngr,lu,n_id)
     $                 /cpool_slow(jpngr,lu,1)
                else
                  tmp_ntoc_slow = 0.0d0
                endif
#endif
                cpool_fast(jpngr,lu,1) = thawCDens*athawdepth_old(jpngr,lu,1)*0.25d0
                cpool_slow(jpngr,lu,1) = thawCDens*athawdepth_old(jpngr,lu,1)*0.25d0
#if _DyN
                cpool_fast(jpngr,lu,n_id) = cpool_fast(jpngr,lu,1)*tmp_ntoc_fast
                cpool_slow(jpngr,lu,n_id) = cpool_slow(jpngr,lu,1)*tmp_ntoc_slow
#endif
              endif
#endif
            endif
          enddo
        endif

        if (spinup.and.year.eq.soil_equil_year+1) then

          call calc_eq_pools(jpngr,.true.)
          
c     re-initialise averages
          k_fast_ave(jpngr,:)=0.d0
          k_slow_ave(jpngr,:)=0.d0
          litter_decom_ave(jpngr,:,:)=0.d0
#if _DyN
          litter_decom_ave_pft(jpngr,:)=0.d0
#endif
#if _peatlands_on
          acro_to_cato_ave(jpngr,:,:)=0.d0
#endif
#if _permafrostc_on
          soil_to_permaf_fast_ave(jpngr,: ,:)=0.0d0
          soil_to_permaf_slow_ave(jpngr,: ,:)=0.0d0
#endif          
        endif

        if (spinup.and.year.eq.soil_equil_year_two+1) then
          
          call calc_eq_pools(jpngr,.false.)
          
        endif

#if _soil_decomp_recycle
c     --------------------------------------------------------------------------
c     AVERAGE OVER RECYCLE PERIOD
c     --------------------------------------------------------------------------
c     use only the last recycle period for averaging of k

        if (spinup.and.(year.le.soil_equil_year).and.
     $       (year.gt.soil_equil_year-recycle)) then
          
          call calc_flux_ave(jpngr,recycle,lu_category)

        endif                                          

        if (spinup.and.(year.le.soil_equil_year_two).and.
     $       (year.gt.soil_equil_year_two-recycle)) then

          call calc_flux_ave(jpngr,recycle,lu_category)

        endif                                          

#else
c     --------------------------------------------------------------------------
c     AVERAGE OVER ENTIRE SPINUP
c     --------------------------------------------------------------------------
c     averaging k over whole spinup period
c     this is generally unnecessary but practically equivalent to 'recycle'

        if (spinup.and.year.le.soil_equil_year) then

          call calc_flux_ave(jpngr,soil_equil_year,lu_category)

        endif
#endif
          
#if _soil_diseq
c     --------------------------------------------------------------------------
c     TEST SOIL DISEQUILIBRIUM
c     --------------------------------------------------------------------------
c     use this code to check equilibrium C content / soil disequilibrium
c     based on last climbuf years
c     (without 13C/14C)
c     -> must be after spinup, of course
        if (year.le.soil_diseq_year) then
            
          if(year.gt.soil_diseq_year-climbuf) then

            call calc_flux_ave(jpngr,climbuf,lu_category)

          elseif(year.eq.soil_diseq_year-climbuf) then
c     re-initialise averages and pseudo-pools
            k_fast_ave(jpngr,:)=0.d0
            k_slow_ave(jpngr,:)=0.d0
            litter_decom_ave(jpngr,:,:)=0.d0
#if _peatlands_on
            acro_to_cato_ave(jpngr,:,:)=0.d0
#endif
#if _permafrostc_on
            soil_to_permaf_fast_ave(jpngr,: ,:)=0.0d0
            soil_to_permaf_slow_ave(jpngr,: ,:)=0.0d0
#endif
            cpool_fast_eq(jpngr,:)=0.d0
            cpool_slow_eq(jpngr,:)=0.d0
          endif
        elseif (year.eq.soil_diseq_year+1) then
          do lu=1,nlu        
c     Analytically calculate pool sizes this year only
            
c     Rate of change of soil pool size = litter input - decomposition
c     (5) dc/dt = litter_decom - kc
c     At equilibrium,
c     (6) dc/dt = 0
c     From (5) & (6),
c     (7) c = litter_decom / k

#if _permafrostc_on
c     Prevent negative soil pools due to numerical and averaging issues
            soil_to_permaf_fast_ave(jpngr,lu,1)=min(soil_to_permaf_fast_ave(jpngr,lu,1),
     $           soilfrac(lu)*fastfrac*litter_decom_ave(jpngr,lu,1)*0.99d0)
            soil_to_permaf_slow_ave(jpngr,lu,1)=min(soil_to_permaf_slow_ave(jpngr,lu,1),
     $           soilfrac(lu)*(1.0d0-fastfrac)*litter_decom_ave(jpngr,lu,1)*0.99d0)
#endif
            
c     d(12Cs)/dt=lit_dec-atoc-12Cs*k => 12Cseq=(lit_dec-atoc)/k
            cpool_fast_eq(jpngr,lu)=(soilfrac(lu)*fastfrac
     $           * litter_decom_ave(jpngr,lu,1)
#if _peatlands_on
     $           - acro_to_cato_ave(jpngr,lu,1)
#endif
#if _permafrostc_on
     $           - soil_to_permaf_fast_ave(jpngr,lu,1)
#endif
     $           ) / k_fast_ave(jpngr,lu)                     !eqn 7

            cpool_slow_eq(jpngr,lu)=(soilfrac(lu)*(1.0d0-fastfrac)
     $           * litter_decom_ave(jpngr,lu,1)
#if _peatlands_on
     $           + acro_to_cato_ave(jpngr,lu,1)
#endif
#if _permafrostc_on
     $           - soil_to_permaf_slow_ave(jpngr,lu,1)
#endif
     $           ) / k_slow_ave(jpngr,lu)                     !eqn 7

          enddo                 !lu
        endif
#if _DyN
        stop '# soil_diseq not yet implemented for _DyN
     $(see littersom SR)'
#endif
          
#endif
#if _peatlands_on
!     update C content of acrotelm and save cpool_acro_old
        call update_acrotelm(jpngr,fpc_grid(mosspft,jpngr),
     $       cpool_fast(jpngr,lupeat,1))  
#endif 
#if _mineralsoils_on
c     Save thawdepth for next year
        athawdepth_old(jpngr,:,2) = athawdepth_old(jpngr,:,1)
        athawdepth_old(jpngr,:,1) = athawdepth(jpngr,:)
#endif        
      endif                     !(m.eq.12)         
     
      return
      end
      

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE LLOYD_AND_TAYLOR
c     Calculates decomposition rate modifier for a given temperature
c     according to Lloyd & Taylor 1994
c     ------------------------------------------------------------------
      subroutine lloyd_and_taylor(pool_temp,ref_temp,modifier)
      
      implicit none

c     ARGUMENTS
      REAL*8 pool_temp          ! actual temperature in °C
      REAL*8 ref_temp           ! reference temperature in °C, normal: 10°C
      REAL*8 modifier           

c     LOCAL VARIABLES
      REAL*8 E0
      REAL*8 T0
      REAL*8 Tzero

      PARAMETER (E0=308.56d0)   ! Activation Energy
      PARAMETER (T0=227.13d0)   ! calibration temperature [K]
      PARAMETER (Tzero=273.15d0) ! 0°C = 273.15 K 

c     CHECK ARGUMENT
c     --------------------------------------------------------------------------
c      if (ref_temp.lt.-40.0d0) then 
c        print*,'Error: Reference temperature in Lloyd and Taylor'
c        print*,'       parametrisation is set below -40°C'
c        stop
c      endif


c     TEMPERATURE MODIFIER
c     --------------------------------------------------------------------------
c     Temperature response function is a modified Q10 relationship
c     (Lloyd & Taylor 1994)
      if (pool_temp.ge.-40.0d0) then ! avoid numerical errors

        modifier=exp(E0*((1.0d0/(ref_temp+Tzero-T0))
     $       -(1.0d0/(pool_temp+Tzero-T0))))

      else
c     set temperature response to a constant at value of -40°C

        modifier=exp(E0*((1.0d0/(ref_temp+Tzero-T0))
     $       -(1.0d0/(-40.0d0+Tzero-T0))))

      endif

      return
      end

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE BROVKIN_Q10
c     Calculates decomposition rate modifier for a given temperature
c     using Q10 according to Brovkin et al., 2012 for woody litter
c     ------------------------------------------------------------------
      subroutine brovkin_Q10(pool_temp,Q10,modifier)
      
      implicit none

c     ARGUMENTS
      REAL*8 pool_temp, Q10, modifier !temperature in °C

c     TEMPERATURE MODIFIER
c     --------------------------------------------------------------------------
c     Temperature response function is a Q10 relationship
c     (Brovkin et al., BG, 2012)
        modifier=Q10**((pool_temp-10.0d0)/10.0d0) 

      return
      end


c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE FOLEY HYDROLOGY
c     Calculates decomposition rate modifier for a given water fraction
c     according to Foley 1995
c     ------------------------------------------------------------------
      subroutine foley(pool_water,modifier)
      
      implicit none

c     ARGUMENTS
      REAL*8 pool_water, modifier !water as fraction of whc

c     Foley Hydrology
c     --------------------------------------------------------------------------
c     Moisture response based on soil layer 1 moisture content (Foley 1995)
      modifier=(0.25d0+(0.75d0*pool_water))

      return
      end


c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE Gerten HYDROLOGY
c     Calculates decomposition rate modifier for a given water fraction
c     according to Gerten et al., 2004
c     ------------------------------------------------------------------
      subroutine gerten(pool_water,modifier)
      
      implicit none

c     ARGUMENTS
      REAL*8 pool_water, modifier !water as fraction of whc


c     Gerten Hydrology
c     --------------------------------------------------------------------------
c     new soil moisture response function: Wania et al., Integrating peatlands 
c     and permafrost into a dynamic global vegetation model: 2. Evaluation and 
c     sensitivity of vegetation and carbon cycle processes, GBC, vol. 23, 2009
c     --------------------------------------------------------------------------
c     Rita: I don't know where the new Eqn below comes from. Compared
c     to Foley's eqn, it slows down the moisture response at low water
c     contents, which seems reasonable.
c     --------------------------------------------------------------------------
c     spahni 2013-01-21: looks like Fang and Moncrieff 1999 approach to me:
c
c     Fang, C. and Moncrieff, J.: A model for soil CO2 production and
c     transport 1: model development, Agr. Forest Meterol., 95, 225–
c     236, 1999.
c     --------------------------------------------------------------------------

      if (pool_water.gt.1.0d-5) then
        modifier=(1.0d0-exp(-pool_water))
     $       / (1.0d0-exp(-1.0d0))
      else
        
c     spahni 2009-09-23: set moisture response to a constant (=1.0d-5) for 
c     very low monthly whc fraction 
        modifier=1.0d-5
        
      endif

      return
      end
      

#if _peatlands_on
c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE FROLKING HYDROLOGY
c     Calculates decomposition rate modifier for a given water fraction
c     according to Frolking et al., 2004 fro peatlands
c     ------------------------------------------------------------------
      subroutine frolking(pool_water,watertpos,modifier)
      
      implicit none
#include "para.inc"
#include "soil.inc"
#include "soildecomp.inc"

c     ARGUMENTS
      REAL*8 pool_water(nacrotelm), watertpos !water as fraction of whc
      REAL*8 modifier    

c     PARAMETERS
      REAL*8 wopt               ! optimal degree of saturation for decomposition
      PARAMETER (wopt=0.45d0)   ! Frolking et al., 2010: 0.45
      REAL*8 cmoist             ! Frolking et al., 2010: 2.31
      REAL*8 csat               ! Frolking et al., 2010: 300mm
      PARAMETER (csat=z1+z2/2.0d0)

c     VARIABLES
      INTEGER acrol               ! loops
      REAL*8 fmin               ! Frolking et al., 2010: 0.3
      REAL*8 wtd                ! monthly water table depth
      REAL*8 fdepth(nacrotelm)  ! fraction of layer to wtd
      REAL*8 fmoist             ! acrotelm moisture response function
      REAL*8 zmid               ! mid depth between wtd and acrotelm bottom
      REAL*8 fsat               ! acrotelm saturated response function

c     Frolking hydrology
c     --------------------------------------------------------------------------
c     Frolking, S., Roulet, N. T., Tuittila, E., Bubier, J. L., Quillet, A., Talbot,
c     J., and Richard, P. J. H.: A new model of Holocene peatland net primary
c     production, decomposition, water balance, and peat accumulation, Earth Syst.
c     Dynam., 1, 1–21, doi:10.5194/esd-1-1-2010, 2010.
c     --------------------------------------------------------------------------

c     cmoist parameter
      cmoist=(1.0d0-Rmoist_soil)/(1.0d0-wopt)**2.0d0

c     peatland moisture response
      modifier = Rmoist_soil
      fmin=k_soil_cato10/k_soil_acro10

c     calculate fraction of saturation for 3 acrotelm layers
      if (watertpos.lt.0.0d0) then
        wtd=-1.0d0*watertpos
        if (wtd.ge.z1) then
          fdepth(:)=Dz_acro/z1
        else if (wtd.ge.2.0d0*Dz_acro) then
          fdepth(1)=Dz_acro/wtd
          fdepth(2)=Dz_acro/wtd
          fdepth(3)=(wtd-2.0d0*Dz_acro)/wtd
        else if (wtd.ge.Dz_acro) then
          fdepth(1)=Dz_acro/wtd
          fdepth(2)=(wtd-Dz_acro)/wtd
          fdepth(3)=0.0d0
        else
          fdepth(1)=1.0d0
          fdepth(2)=0.0d0
          fdepth(3)=0.0d0                
        endif

        fmoist=0.0d0
        do acrol=1,nacrotelm
c     Moisture response in undersaturated acrotelm according to Frolking, 2010
c     --------------------------------------------------------------------------
c     fmoist=1.0d0-cmoist*(pool_water(acrol)-wopt)**2
          fmoist=fmoist+(1.0d0-cmoist*(pool_water(acrol)-wopt)**2.0d0)
     $         *fdepth(acrol)                
        enddo
c     Variable scaling for saturated acrotelm
        zmid=(z1+wtd)/2.0d0
        fsat=fmin+(Rmoist_soil-fmin)*exp(-1.0d0*(zmid-wtd)/csat)
c     Weighted mean
        modifier=wtd/z1*fmoist+(z1-wtd)/z1*fsat
      endif

      return
      end
#endif

#if _sublayer_hydhd
c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE SOIL_CARBON_PROFILE
c     Calculates the soil carbon distribution in the top 3m of soil for
c     mineral soils in the presence or absence of permafrost
c     -------------------------------------------------------------------------
      subroutine soil_carbon_profile(jpngr,lu,totsoilC,thawdepth,Cdistr)
      
      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "soil.inc"
      include 'permafrost.inc'
      
c     ARGUMENTS
      INTEGER jpngr, lu
      REAL*8 totsoilC           !total soil carbon [gC/m2]
      REAL*8 thawdepth          !daily or monthly thawdepth [mm]
      REAL*8 Cdistr(nsl)        !carbon distribution

c     VARIABLES
      INTEGER ns
c     INTEGER pflay
      REAL*8 Cdistrm(nsl)       !carbon distribution in mineral unfrozen soils
      REAL*8 alscale            !scaling length for active layer
      REAL*8 pfdepth            !permafrost depth
      REAL*8 pfdeep             !deepest permafrost depth

c     ARRAYS
      DATA (Cdistrm(ns), ns=1,nsl) /0.108677d0, 0.082055d0, 0.069677d0,
     $     0.116576d0, 0.096925d0, 0.120067d0, 0.153384d0, 0.252638d0/
c     from:
c     
c     Yugang Wang, Yan Li, Xuehua Ye, Yu Chu, Xinping Wang, Profile storage of
c     organic/inorganic carbon in soil: From forest to desert, Science of The
c     Total Environment, Volume 408, Issue 8, 15 March 2010, Pages 1925-1931,
c     ISSN 0048-9697, 10.1016/j.scitotenv.2010.01.015.
c     
c     and LPX layer depth intervals: Dz

c     Calibrate soil C distribution
c     -------------------------------------------------------------------------
c     Decision tree:
c     
c     1) Is permafrost present and thus thawdepth < 2m?
c     2) Is permafrost all through to the surface?
c     3) Is permafrost depth greater than 2m?
c     4) Is C allocatable in the active layer?
c     5) Is there C left for allocation in permafrost?

c            if (totsoilC.gt.20000.0d0.and.permafrostyr(jpngr,lu).ge.2) then
c              print*,'its here'
c            endif

#if _permafrostc_on      
c     1) yes
      if (permafrostyr(jpngr,lu).ge.2.and.thawdepth.lt.d1+d2) then
c     2) yes
        if (thawdepth.lt.soil_middepth(1)) then
          Cdistr(:)=0.0d0
          Cdistr(1)=1.0d0
c          pfdepth=totsoilC/PermDens ![gC/m2/(kgC/m3)=mm]
cc     3) yes
c          if (pfdepth.gt.soil_middepth(nsl)) then
c            Cdistr(:)=Dz_soil(:)/pfdepth
c            Cdistr(nsl)=1.0d0-sum(Cdistr(1:nsl-1))
cc     3) no
c          else
c            Cdistr(:)=0.0d0
c            Cdistr(1)=Dz_soil(1)
c            do ns=2,nsl
c              if (pfdepth.gt.soil_middepth(ns)) then
c                Cdistr(ns)=Dz_soil(ns)
c              else
c                exit
c              endif
c            enddo
c            Cdistr(:)=Cdistr(:)/sum(Cdistr(:))
c          endif
cc     2) no
        else
          alscale=ethaw/thawdepth
          Cdistr(:)=0.0d0
          Cdistr(1)=SurfDens*exp(alscale*soil_middepth(1))
     $         *Dz_soil(1)
c          pflay=nsl
          do ns=2,nsl
c     4) yes
            if (thawdepth.gt.soil_middepth(ns)) then
              Cdistr(ns)=SurfDens*exp(alscale*soil_middepth(ns))
     $             *Dz_soil(ns)
c     4) no
            else
c              pflay=ns
              exit
            endif
          enddo
c          pfdepth=(totsoilC-sum(Cdistr(1:pflay-1)))/PermDens
cc     5) no
c          if (pfdepth.le.Dz_soil(pflay)/2.0d0) then
c            Cdistr(:)=Cdistr(:)/sum(Cdistr(:))
cc     5) yes
c          else
c            Cdistr(1:pflay-1)=Cdistr(1:pflay-1)
c     $           /sum(Dz_soil(1:pflay-1))
c            pfdeep=soil_middepth(pflay)+pfdepth-Dz_soil(pflay)/2.0d0
cc     3) yes
c            if (pfdeep.gt.soil_middepth(nsl)) then
c              Cdistr(pflay:nsl)=Dz_soil(pflay:nsl)/pfdepth
c              if (pflay.eq.nsl) then
c                Cdistr(nsl)=1.0d0
c              else
c                Cdistr(nsl)=1.0d0-sum(Cdistr(pflay:nsl-1))
c              endif
c              Cdistr(pflay:nsl)=Cdistr(pflay:nsl)*PermDens
cc     3) no
c            else
c              Cdistr(pflay:nsl)=0.0d0
c              do ns=pflay,nsl
c                if (pfdeep.gt.soil_middepth(ns)) then
c                  Cdistr(ns)=Dz_soil(ns)*PermDens
c                else
c                  Cdistr(pflay:ns-1)=Cdistr(pflay:ns-1)
c     $                 /sum(Dz_soil(pflay:ns-1))
c                  exit
c                endif
c              enddo
c            endif
            Cdistr(:)=Cdistr(:)/sum(Cdistr(:))
c          endif  
        endif
c     1) no
      else
        Cdistr(:)=Cdistrm(:)
      endif
#else
      Cdistr(:)=Cdistrm(:)
#endif      
      return
      end
#endif

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE LITTER MODIFIER
c     Calculates decomposition rate modifier for temperature and moisture
c     -------------------------------------------------------------------------
      subroutine litter_modifier(woody,lu,pft,soilt,soilw,respir_modifier)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "litterdecomp.inc"

c     ARGUMENTS
      LOGICAL woody
      INTEGER lu,pft
      REAL*8 soilt
      REAL*8 soilw
      REAL*8 respir_modifier
      
c     VARIABLES
      REAL*8 temp_respir        ! monthly temperature response of decomposition
      REAL*8 moist_respir       ! monthly moisture response of decomposition

      temp_respir=0.0d0
      moist_respir=0.0d0

c     TEMPERATURE MODIFIER
c     --------------------------------------------------------------------------
#if _soil25cm_resp
c     For 25cm temperature always use Lloyd & Taylor:
c     Temperature response function is a modified Q10 relationship
c     (Lloyd & Taylor 1994)
c     reference temperature: 10°C
      call lloyd_and_taylor(soilt,10.0d0,temp_respir)
#else
      if (woody) then
c     Replace Lloyd and Taylor with Q10 parametrisation for woody PFT litter
c     (Brovkin et al., 2012)
        call brovkin_Q10(soilt,Q10_woody(pft),temp_respir)        
      else
c     Temperature response function is a modified Q10 relationship
c     (Lloyd & Taylor 1994)
c     reference temperature: 10°C
        call lloyd_and_taylor(soilt,10.0d0,temp_respir)
      endif
#endif

c     MOISTURE MODIFIER
c     --------------------------------------------------------------------------
#if _sublayer_hydhd
c     Gerten Hydrology
c     --------------------------------------------------------------------------
c     new soil moisture response function: Wania et al., Integrating peatlands 
c     and permafrost into a dynamic global vegetation model: 2. Evaluation and 
c     sensitivity of vegetation and carbon cycle processes, GBC, vol. 23, 2009
      if (lu.eq.lupeat) then
        moist_respir=Rmoist_litter
      else
        call gerten(soilw,moist_respir)
      endif
#else
#if _sitch_hydrology
c     Foley Hydrology
c     --------------------------------------------------------------------------
c     Moisture response based on soil layer 1 moisture content (Foley 1995)
      call foley(soilw,moist_respir)
#else
      call gerten(soilw,moist_respir)
#endif
#endif

      respir_modifier=temp_respir*moist_respir
        
      return
      end


#if _sublayer_hydhd
c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE Soil Temperature and Moisture MODIFIER
c     Calculates decomposition rate modifier
c     ------------------------------------------------------------------
      subroutine soil_modifier_layers(jpngr,lu,np,soilt,soilw,watertpos,
     $     totsoilC,thawdepth,respir_modifier)
      
      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "soil.inc"
#include "soildecomp.inc"
      
c     ARGUMENTS
      INTEGER jpngr,lu, np
      REAL*8 soilt(nsl)
      REAL*8 soilw(nsl)
      REAL*8 watertpos
      REAL*8 totsoilC
      REAL*8 thawdepth
      REAL*8 respir_modifier

c     VARIABLES
      INTEGER l
      REAL*8 temp_respir        !monthly temperature response of decomposition
      REAL*8 moist_respir       !monthly moisture response of decomposition
      REAL*8 temp_respir_ave    !monthly temperature response of decomposition
      REAL*8 moist_respir_ave   !monthly moisture response of decomposition
      REAL*8 Cdistr(nsl)        !carbon distribution

c     loop over pools: fast/slow soil pool or acro-/catotelm pool
      
c     apply weighted mean temperature and moisture
c     --------------------------------------------------------------------------         
      temp_respir=0.0d0
      temp_respir_ave=0.0d0
      moist_respir=0.0d0
      moist_respir_ave=0.0d0

#if _peatlands_on
      if (lu.eq.lupeat) then
c     acrotelm: top 3 soil layers
        if (np.eq.1) then
c     temperature modifier
          do l=1,nacrotelm
c     reference temperature: 10°C            
            call lloyd_and_taylor(soilt(l),10.0d0,temp_respir)
            temp_respir_ave=temp_respir_ave+temp_respir
     $           *Dz_soil(l)/z1
          enddo

c     MOISTURE MODIFIER
c     --------------------------------------------------------------------------
c     soil moisture response function for peatland acrotelm described in
c
c     Spahni et al., 2012, Transient simulations of the carbon and nitrogen
c     dynamics in northern peatlands: from the Last Glacial Maximum to the 21st
c     century, Clim. Past Discuss., 8, 5633-5685, doi:10.5194/cpd-8-5633-2012,
c     2012.
c     --------------------------------------------------------------------------
c     following an approach by Frolking et al., 2010
          call frolking(soilw(1:nacrotelm),watertpos,moist_respir_ave)
          
c     catotelm: lower 5 soil layers
        else                    !np=2
c     temperature modifier
          do l=nacrotelm+1,nacrotelm+ncatotelm
c     reference temperature: 10°C
            call lloyd_and_taylor(soilt(l),10.0d0,temp_respir)
            temp_respir_ave=temp_respir_ave+temp_respir
     $           *Dz_soil(l)/z2
          enddo
          
c     moisture modifier
          moist_respir_ave=Rmoist_soil
        endif
        
      else                      ! not lupeat
#endif

c     Calculate weight function from soil carbon profile, sum(Cdistr)=1.0d0
        Cdistr(:)=0.0d0
        call soil_carbon_profile(jpngr,lu,totsoilC,thawdepth,Cdistr)
        
c        if (abs(sum(Cdistr(:))-1.0d0).gt.10.0d-5) then
c          print*,'bad distribution'
c        endif

        if (np.eq.1) then
          do l=1,nsl
c     reference temperature: 10°C
            call lloyd_and_taylor(soilt(l),10.0d0,temp_respir)
            temp_respir_ave=temp_respir_ave+temp_respir
     $           *Cdistr(l)
            call gerten(soilw(l),moist_respir)
            moist_respir_ave=moist_respir_ave+moist_respir
     $           *Cdistr(l)
          enddo
        else
c     do excactly the same for slow C pool in lu.ne.lupeat
          do l=1,nsl
c     reference temperature: 10°C
            call lloyd_and_taylor(soilt(l),10.0d0,temp_respir)
            temp_respir_ave=temp_respir_ave+temp_respir
     $           *Cdistr(l)
            call gerten(soilw(l),moist_respir)
            moist_respir_ave=moist_respir_ave+moist_respir
     $           *Cdistr(l)
          enddo
        endif
#if _peatlands_on
      endif
#endif
      respir_modifier=temp_respir_ave*moist_respir_ave
      
      return
      end

#else

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE Soil Temperature and Moisture MODIFIER
c     Calculates decomposition rate modifier
c     ------------------------------------------------------------------
      subroutine soil_modifier(soilt,soilw,respir_modifier)
      
      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "soil.inc"
#include "soildecomp.inc"
      
c     ARGUMENTS
      REAL*8 soilt
      REAL*8 soilw
      REAL*8 respir_modifier

c     VARIABLES
      REAL*8 temp_respir        !monthly temperature response of decomposition
      REAL*8 moist_respir       !monthly moisture response of decomposition
      
c     temperature and moisture
c     --------------------------------------------------------------------------         
      temp_respir=0.0d0
      moist_respir=0.0d0

c     reference temperature: 10°C
      call lloyd_and_taylor(soilt,10.0d0,temp_respir)
#if _sitch_hydrology
      call foley(soilw,moist_respir)
#else
      call gerten(soilw,moist_respir)
#endif
      respir_modifier=temp_respir*moist_respir
      
      return
      end
#endif

#if _permafrostc_on
c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE BUILD PERMAFROST CARBON
c     Calculates if permfrost C is sequestered
c     -------------------------------------------------------------------------
      subroutine build_permafrost_carbon(jpngr,lu,thawdepth,STPF,STPS)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "pools.inc"
      include 'permafrost.inc'

c     ARGUMENTS
      INTEGER jpngr             ! Grid cell
      INTEGER lu                ! Land unit
      REAL*8 thawdepth          ! Thawdepth in permafrost area [mm]
      REAL*8 STPF(allvar)       ! soil to permafrost transfer rate for fast pool
      REAL*8 STPS(allvar)       ! soil to permafrost transfer rate for slow pool

c     VARIABLES
      REAL*8 TSC                ! Total Soil Carbon [kgC/m2]
      REAL*8 THLC               ! Total Carbon in thaw layer [kgC/m2]
      REAL*8 bf                 ! buildfraction of C
         
c     Calculate max carbon content in thaw layer from annual thaw depth
      THLC = thawCDens*thawdepth/1000.0d0
c     Calculate total actual soil carbon in kgC/m2
      TSC = (cpool_fast(jpngr,lu,1) + cpool_slow(jpngr,lu,1))/1000.0d0
c     Calculate min permafrost formation from loess deposition
      if (TSC.gt.0.0d0) then
        bf = ldep * PermDens / TSC
      else
        bf = 0.0d0
      endif

      if (TSC.gt.THLC) then
c     Bury additionally 35% of excess carbon in permafrost
        bf = bf + (TSC - THLC) / TSC * 0.35d0 
      endif
      
      STPF(1)=bf*cpool_fast(jpngr,lu,1)/dble(ntimesteps)
      STPF(2:ncvar)=cpool_fast(jpngr,lu,2:ncvar)
#if _DyN
      STPF(n_id)=bf*cpool_fast(jpngr,lu,n_id)/dble(ntimesteps)
#endif
      STPS(1)=bf*cpool_slow(jpngr,lu,1)/dble(ntimesteps)
      STPS(2:ncvar)=cpool_slow(jpngr,lu,2:ncvar)
#if _DyN
      STPS(n_id)=bf*cpool_slow(jpngr,lu,n_id)/dble(ntimesteps)
#endif

      return
      end
#endif
#if _permafrostc_on
c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE THAW PERMAFROST CARBON
c     Calculates if permafrost C is thawed and moved back to soil pool
c     -------------------------------------------------------------------------
      subroutine thaw_permafrost_carbon(jpngr,lu,thawdepth,thawdepth_old,STPF,STPS)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "pools.inc"
      include 'permafrost.inc'

c     ARGUMENTS
      INTEGER jpngr             ! Grid cell
      INTEGER lu                ! Land unit
      REAL*8 thawdepth          ! Thawdepth in permafrost area [mm]
      REAL*8 thawdepth_old      ! Thawdepth in permafrost area in previous time step [mm]
      REAL*8 STPF(allvar)       ! soil to permafrost transfer rate for fast pool
      REAL*8 STPS(allvar)       ! soil to permafrost transfer rate for slow pool

c     VARIABLES
      REAL*8 TPC                ! Total Permafrost Carbon [kgC/m2]
      REAL*8 dthd               ! Difference in thawdepth > 0 [m]
      REAL*8 pfd                ! Permafrost depth > 0 [m]
      REAL*8 tf                 ! thawfraction of C
     
      TPC = (cpool_permaf_fast(jpngr,lu,1)+cpool_permaf_slow(jpngr,lu,1))/1000.0d0
      dthd = (thawdepth - thawdepth_old)/1000.0d0
      pfd = max(TPC/PermDens,0.001d0)

      if (dthd.lt.pfd) then
        tf=dthd/pfd
      else
        tf=1.0d0
      endif
      
      STPF(1)=tf*cpool_permaf_fast(jpngr,lu,1)/dble(ntimesteps)
      STPF(2:ncvar)=cpool_permaf_fast(jpngr,lu,2:ncvar)
#if _DyN
      STPF(n_id)=tf*cpool_permaf_fast(jpngr,lu,n_id)/dble(ntimesteps)
#endif
      STPS(1)=tf*cpool_permaf_slow(jpngr,lu,1)/dble(ntimesteps)
      STPS(2:ncvar)=cpool_permaf_slow(jpngr,lu,2:ncvar)
#if _DyN
      STPS(n_id)=tf*cpool_permaf_slow(jpngr,lu,n_id)/dble(ntimesteps)
#endif

      return
      end
#endif

#if _peatlands_on
!//////////////////////////////////////////////////////////////////////////////
!     SUBROUTINE CALCULATE ACROTELM DENSITY 
!     Calculates the average acrotelm density
!     -------------------------------------------------------------------------
      subroutine update_acrotelm(jpngr,fpc_moss,cpool_fast)

      implicit none
#include "para.inc"
#include "soil.inc"

!     ARGUMENTS
      INTEGER jpngr
      REAL*8 fpc_moss
      REAL*8 cpool_fast

!     Bulk peat density (OM density) is calculated from Belyea and Clymo, 2001
!     based on density (y) curve for acrotelm depth (x):
!
!     y = (0.108+0.0057*x)^2
!
!     This is averaged from 0 to 30cm depth:
!     
!     y_mid = 1/30cm * INTEGRAL(y dx)[30][0] 
!           = 1/30 * 1/3 * 1/0.0057 * ((0.108+0.0057*30)^3 - (0.108+0.0057*0)^3)
!           = 0.039879d0
!     

!     LOCAL PARAMETERS
!     acrotelm ash-free bulk density (g OM/cm³)
      REAL*8 bulk_density            
      PARAMETER(bulk_density=0.039879d0)
!     C content in OM for sphagnum (C%/100%) from Loisel et al., 2014
      REAL*8 c_cont_sphagnum
      PARAMETER(c_cont_sphagnum=0.420d0) 
!     C content in OM for other plants (C%/100%)
      REAL*8 c_cont_other
      PARAMETER(c_cont_other=0.502d0) 
    
!     LOCAL VARIABLES
      REAL*8 added_C_cont       !memory variable: C content of added OM (C%/100%)
      REAL*8 cpool_acro_diff    !memory variable: added mass of acrotelm C (gC/m²)
      REAL*8 acro_C_content_old !memory variable: old C content of OM (C%/100%)

!     Calculat the C content of added OM depending on moss fpc
      added_C_cont = c_cont_sphagnum*fpc_moss+c_cont_other*(1.0d0-fpc_moss)

!     Calculate added C in acrotelm this year:
!     This is an approximation for strongly varying acrotelm size, but exact
!     in steady state. 
      cpool_acro_diff = (cpool_acro_old(jpngr)-cpool_acro_norm(jpngr))

!     Calculate weighted average C content in acrotelm for positiv mass balance
      if (cpool_acro_diff.gt.0.0d0) then
        acro_C_content_old = acro_C_content(jpngr)
        acro_C_content(jpngr) = (cpool_acro_norm(jpngr)*acro_C_content_old
     $       +cpool_acro_diff*added_C_cont)/cpool_acro_old(jpngr)
      endif

!     Calculate acrotelm C reference pool (gC/m²) with following multiplications
!     bulk density: g OM/cm³ -> g OM/m²/cm -> *10000
!     acrotelm depth (z1);  300mm -> 30cm -> /10
!     C content (-): -> *1
!     total g OM/cm³ -> gC/m² -> *1000
      cpool_acro_norm(jpngr) = bulk_density*z1*acro_C_content(jpngr)*1000.0d0
      
!     Save fast pool in peatlands (=acrotelm) for next year
      cpool_acro_old(jpngr)=cpool_fast
      
      end
#endif

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE CALCULATE EQUILIBRIUM CARBON POOLS
c     Calculates the equilibrium of C pools
c     -------------------------------------------------------------------------
      subroutine calc_eq_pools(jpngr,eq_peatlands)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "pools.inc"
#include "litterdecomp.inc"
#include "soildecomp.inc"
      include 'permafrost.inc'

c     ARGUMENTS
      INTEGER jpngr             ! Grid cell
      LOGICAL eq_peatlands

c     VARIABLES
      INTEGER lu                ! Land unit
      REAL*8 tmp_ntoc_fast      ! C/N ratio in fast soil pool
      REAL*8 tmp_ntoc_slow ! C/N ratio in slow soil pool

      REAL*8 f_in_slow          ! average input flux into fast soil pool (acrotelm), auxiliary for readability
      REAL*8 f_in_fast          ! average input flux into slow soil pool (acrotelm), auxiliary for readability

c     Analytically calculate pool sizes this year only
            
c     Rate of change of soil pool size = litter input - decomposition
c     (5) dc/dt = litter_decom - kc
c     At equilibrium,
c     (6) dc/dt = 0
c     From (5) & (6),
c     (7) c = litter_decom / k
          
c     spahni 2010-09-24: additionally for peatlands the slow cpool recieves r=12 gC/m2/yr
c     -> (5) dc/dt = litter_decom - kc + r
c     -> (7) c = (litter_decom+r) / k
c     resp. for slow
c     -> (7) c = (litter_decom-r) / k
      do lu=1,nlu
        if (eq_peatlands.or.(lu.ne.lupeat.and.lu.ne.lupeatold)) then
#if _DyN
c     Store N/C ratio before soil equilibration to impose same ratio after equilibration.
          if (cpool_fast(jpngr,lu,1).gt.0.d0) then
            tmp_ntoc_fast = cpool_fast(jpngr,lu,n_id)
     $           /cpool_fast(jpngr,lu,1)
          else
            tmp_ntoc_fast = 0.0d0
          endif
          if (cpool_slow(jpngr,lu,1).gt.0.d0) then
            tmp_ntoc_slow = cpool_slow(jpngr,lu,n_id)
     $           /cpool_slow(jpngr,lu,1)
          else
            tmp_ntoc_slow = 0.0d0
          endif
#endif
          
#if _permafrostc_on
c     Prevent negative soil pools due to numerical and averaging issues
          soil_to_permaf_fast_ave(jpngr,lu,1)=min(soil_to_permaf_fast_ave(jpngr,lu,1),
     $         soilfrac(lu)*fastfrac*litter_decom_ave(jpngr,lu,1)*0.99d0)
          soil_to_permaf_slow_ave(jpngr,lu,1)=min(soil_to_permaf_slow_ave(jpngr,lu,1),
     $         soilfrac(lu)*(1.0d0-fastfrac)*litter_decom_ave(jpngr,lu,1)*0.99d0)
#endif
        
c     --------------------------------------------------------------------------
c     SOIL CARBON - C12
c     --------------------------------------------------------------------------
c     d(12Cs)/dt=lit_dec+atoc-12Cs*k => 12Cseq=(lit_dec+atoc)/k

!     Calculate mean input flux into slow and fast soil pool
          f_in_fast = soilfrac(lu)*fastfrac
     $         * litter_decom_ave(jpngr,lu,1)
#if _peatlands_on
     $         - acro_to_cato_ave(jpngr,lu,1)
#endif
#if _permafrostc_on
     $         - soil_to_permaf_fast_ave(jpngr,lu,1)
#endif

          f_in_slow = soilfrac(lu)*(1.0d0-fastfrac)
     $         * litter_decom_ave(jpngr,lu,1)
#if _peatlands_on
     $         + acro_to_cato_ave(jpngr,lu,1)
#endif
#if _permafrostc_on
     $         - soil_to_permaf_slow_ave(jpngr,lu,1)
#endif

#if _peatlands_on

!     Prevent negative inflows leading to negative eq. pools
          f_in_fast=max(0.0d0,f_in_fast)
          f_in_slow=max(0.0d0,f_in_slow)

          if (lu.eq.lupeat.and.pt_equil.ge.0) then
!     Calculate analytically the pool size after N yr of constant input flux
            cpool_fast(jpngr,lu,1) =
     $           (1.d0-exp(-k_fast_ave(jpngr,lu)*dble(pt_equil)))
     $           * f_in_fast / k_fast_ave(jpngr,lu)  
            cpool_slow(jpngr,lu,1) =
     $           (1.d0-exp(-k_slow_ave(jpngr,lu)*dble(pt_equil)))
     $           * f_in_slow / k_slow_ave(jpngr,lu)
          else
#endif
!     Calculate analytical soil C equilibrium (inf. number of constant input flux)
            cpool_fast(jpngr,lu,1) = f_in_fast / k_fast_ave(jpngr,lu)
            cpool_slow(jpngr,lu,1) = f_in_slow / k_slow_ave(jpngr,lu)
#if _peatlands_on
          endif
#endif
          
#if compc13
c     --------------------------------------------------------------------------
c     SOIL CARBON - C13
c     --------------------------------------------------------------------------
#if _peatlands_on
c     d(13C)/dt = R13lit*lit_dec - R13soil*12Csoil*k
c     => R13soil_eq = R13lit*lit_dec / (12Csoil*k) = R13lit
c     => d13soil_eq = d13lit
c     for fast pool the loss through atoc has the identical signature as litter input
          cpool_fast(jpngr,lu,2)=litter_decom_ave(jpngr,lu,2)
c     for the slow pool the carbon input has the weighted signal of litter and atoc input
          if (litter_decom_ave(jpngr,lu,1)+acro_to_cato_ave(jpngr,lu,1)
     $         .gt.0.0d0) then
            cpool_slow(jpngr,lu,2)=(litter_decom_ave(jpngr,lu,2)
     $           * litter_decom_ave(jpngr,lu,1)
     $           + acro_to_cato_ave(jpngr,lu,2)
     $           * acro_to_cato_ave(jpngr,lu,1))
     $           / (litter_decom_ave(jpngr,lu,1)
     $           + acro_to_cato_ave(jpngr,lu,1))
          else
            cpool_slow(jpngr,lu,2)=litter_decom_ave(jpngr,lu,2)
          endif
#else          
c     d(13C)/dt = R13lit*lit_dec - R13soil*12Csoil*k
c     => R13soil_eq = R13lit*lit_dec / (12Csoil*k) = R13lit
c     => d13soil_eq = d13lit
c     Assumption: Permafrost C has no impact on 13C isotopes during spinup
          cpool_fast(jpngr,lu,2)=litter_decom_ave(jpngr,lu,2)
          cpool_slow(jpngr,lu,2)=litter_decom_ave(jpngr,lu,2)
#endif
#endif
        
#if compc14
c     --------------------------------------------------------------------------
c     SOIL CARBON - C14
c     --------------------------------------------------------------------------
#if _peatlands_on
c     d(14C)/dt = R14lit*lit_dec - R14soil*12Csoil*k - R14soil*12Csoil*1/tau
c     => R14soil_eq = R14lit*lit_dec/(12Csoil*(k+1/tau)) = R14lit/(1+1/(k*tau))
c     => d14soil_eq = (d14lit + 1000) / (1+1/(k*tau)) - 1000
c     for fast pool the loss through atoc has the identical signature as litter input
          cpool_fast(jpngr,lu,3)=(litter_decom_ave(jpngr,lu,3)+1000.d0)/
     $         (1.d0+1.d0 / (k_fast_ave(jpngr,lu)*tau14)) - 1000.d0
c     for the slow pool the carbon input has the weighted signal of litter and atoc input
          if (litter_decom_ave(jpngr,lu,1)+acro_to_cato_ave(jpngr,lu,1)
     $         .gt.0.0d0) then
            cpool_slow(jpngr,lu,3)=((litter_decom_ave(jpngr,lu,3)
     $           *litter_decom_ave(jpngr,lu,1)
     $           +acro_to_cato_ave(jpngr,lu,3)
     $           *acro_to_cato_ave(jpngr,lu,1))
     $           /(litter_decom_ave(jpngr,lu,1)
     $           +acro_to_cato_ave(jpngr,lu,1)) + 1000.d0)/
     $           (1.d0+1.d0 / (k_slow_ave(jpngr,lu)*tau14)) - 1000.d0
          else
            cpool_slow(jpngr,lu,3)=(litter_decom_ave(jpngr,lu,3)
     $           +1000.d0)/
     $           (1.d0+1.d0 / (k_slow_ave(jpngr,lu)*tau14)) - 1000.d0
          endif
#else          
c     d(14C)/dt = R14lit*lit_dec - R14soil*12Csoil*k - R14soil*12Csoil*1/tau
c     => R14soil_eq = R14lit*lit_dec/(12Csoil*(k+1/tau)) = R14lit/(1+1/(k*tau))
c     => d14soil_eq = (d14lit + 1000) / (1+1/(k*tau)) - 1000
c     Assumption: Permafrost C has no impact on 14C isotopes during spinup
          cpool_fast(jpngr,lu,3)=(litter_decom_ave(jpngr,lu,3)+1000.d0)/
     $         (1.d0+1.d0 / (k_fast_ave(jpngr,lu)*tau14)) - 1000.d0
          cpool_slow(jpngr,lu,3)=(litter_decom_ave(jpngr,lu,3)+1000.d0)/
     $         (1.d0+1.d0 / (k_slow_ave(jpngr,lu)*tau14)) - 1000.d0
#endif
#endif
          
#if _DyN
c     --------------------------------------------------------------------------
c     SOIL ORGANIC NITROGEN
c     --------------------------------------------------------------------------
c     Use same ratio as at current
          cpool_fast(jpngr,lu,n_id) = cpool_fast(jpngr,lu,1)
     $         * tmp_ntoc_fast
          cpool_slow(jpngr,lu,n_id) = cpool_slow(jpngr,lu,1)
     $         * tmp_ntoc_slow
#endif 

#if _permafrostc_on
c     --------------------------------------------------------------------------
c     PERMAFROST CARBON - C12 (not in equilibrium)
c     --------------------------------------------------------------------------
c     12Cseq=s_to_p*PermBuildTime
          if (lu.ne.lupeat.and.lu.ne.lupeatold) then
            if (soil_to_permaf_fast_ave(jpngr,lu,1).gt.0.0d0) then
              cpool_permaf_fast(jpngr,lu,1)=soil_to_permaf_fast_ave(jpngr,lu,1)
     $             *PermBuildTime
            else
              cpool_permaf_fast(jpngr,lu,1)=0.0d0
            endif
            
            if (soil_to_permaf_slow_ave(jpngr,lu,1).gt.0.0d0) then
              cpool_permaf_slow(jpngr,lu,1)=soil_to_permaf_slow_ave(jpngr,lu,1)
     $             *PermBuildTime
            else
              cpool_permaf_fast(jpngr,lu,1)=0.0d0
            endif
#if compc13
c     --------------------------------------------------------------------------
c     PERMAFROST CARBON - C13
c     --------------------------------------------------------------------------
            cpool_permaf_fast(jpngr,lu,2)=cpool_fast(jpngr,lu,2)
            cpool_permaf_slow(jpngr,lu,2)=cpool_slow(jpngr,lu,2)
#endif
#if compc14
c     --------------------------------------------------------------------------
c     PERMAFROST CARBON - C14
c     --------------------------------------------------------------------------
c     d(14C)/dt = R14soil*s_to_p - R14permaf*12Cpermaf*1/tau
c     => R14permaf_eq = d(14C)/dt*PBT/12Cpermaf 
c                     = R14soil*s_to_p/12Cpermaf*PBT - R14permaf_eq*PBT/tau
c                     = R14soil*s_to_p/12Cpermaf*PBT/(1+PBT/tau)
c     => d14permaf_eq = (d14soil + 1000)*s_to_p/12Cpermaf*PBT/(1+PBT/tau) - 1000
c                     = (d14soil + 1000)/(1+PBT/tau) - 1000
            cpool_permaf_fast(jpngr,lu,3)=(cpool_fast(jpngr,lu,3)+1000.0d0)
     $           /(1.0d0+PermBuildTime/tau14) - 1000.0d0
            cpool_permaf_slow(jpngr,lu,3)=(cpool_slow(jpngr,lu,3)+1000.0d0)
     $           /(1.0d0+PermBuildTime/tau14) - 1000.0d0
#endif
#if _DyN
c     --------------------------------------------------------------------------
c     PERMAFROST ORGANIC NITROGEN
c     --------------------------------------------------------------------------
c     Use same ratio as at current in soil pools
            cpool_permaf_fast(jpngr,lu,n_id) = cpool_permaf_fast(jpngr,lu,1)
     $           * tmp_ntoc_fast
            cpool_permaf_slow(jpngr,lu,n_id) = cpool_permaf_slow(jpngr,lu,1)
     $           * tmp_ntoc_slow
#endif
          
          endif            
#endif
        endif
      enddo
      
      return
      end

c//////////////////////////////////////////////////////////////////////////////
c     SUBROUTINE AVERAGE C FLUXES OVER A GIVEN PERIODO
c     Calculates flux averages
c     -------------------------------------------------------------------------
      subroutine calc_flux_ave(jpngr,period,lu_category)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "pools.inc"
#include "litterdecomp.inc"
#include "soildecomp.inc"
      include 'coordinates.inc'
      include 'permafrost.inc'

c     ARGUMENTS
      INTEGER jpngr             ! Grid cell
      INTEGER period            ! Averaging period
      INTEGER lu_category(npft) ! Look up mask for lu class

c     VARIABLES
      INTEGER lindex,ppft       ! local indexes
      INTEGER lu                ! Land unit

c     Update running average respiration rates and litter input

      do lu=1,nlu
        do lindex=1,ntimesteps
c     --------------------------------------------------------------------------
          k_fast_ave(jpngr,lu)=k_fast_ave(jpngr,lu)+k_fast(lindex,lu)/
     $         dble(period)
          k_slow_ave(jpngr,lu)=k_slow_ave(jpngr,lu)+k_slow(lindex,lu)/
     $         dble(period)
          call addc(jpngr,
     $         litter_decom_ave(jpngr,lu,:),
     $         litter_decom(lindex,lu,1)/dble(period),
     $         litter_decom(lindex,lu,2:ncvar)
     $         )
c     --------------------------------------------------------------------------
#if _peatlands_on
          if (lu.eq.lupeat.and.ptgridcell(jpngr)) then

            acro_to_cato_ave(jpngr,lu,1) = acro_to_cato_ave(jpngr,lu,1)
     $           + acro_to_cato_ann(jpngr)/dble(period)/dble(ntimesteps)
            acro_to_cato_ave(jpngr,lu,2:ncvar) = acro_to_cato(2:ncvar)
            
          endif
#endif
#if _permafrostc_on
          if (lu.ne.lupeat.and.lu.ne.lupeatold) then
            if (soil_to_permaf_fast_ave(jpngr,lu,1).gt.
     $           -soil_to_permaf_fast(lu,1)/dble(period)) then
              call addc(jpngr,
     $             soil_to_permaf_fast_ave(jpngr,lu,:),
     $             soil_to_permaf_fast(lu,1)/dble(period),
     $             soil_to_permaf_fast(lu,2:ncvar)
     $             )
            else
              soil_to_permaf_fast_ave(jpngr,lu,1)=0.0d0
            endif
            
            if (soil_to_permaf_slow_ave(jpngr,lu,1).gt.
     $           -soil_to_permaf_slow(lu,1)/dble(period)) then
              call addc(jpngr,
     $             soil_to_permaf_slow_ave(jpngr,lu,:),
     $             soil_to_permaf_slow(lu,1)/dble(period),
     $             soil_to_permaf_slow(lu,2:ncvar)
     $             )
            else
              soil_to_permaf_slow_ave(jpngr,lu,1)=0.0d0
            endif

          endif
#endif
#if _DyN
          do ppft=1,npft
            if (lu_category(ppft).eq.lu) then
              litter_decom_ave_pft(jpngr,ppft) =
     $             litter_decom_ave_pft(jpngr,ppft)+
     $             litter_decom_pft(lindex,ppft,1)/dble(period)
            endif
          enddo
#endif
        enddo
      enddo                     !lu
      
      return
      end
