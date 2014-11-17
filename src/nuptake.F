      subroutine nuptake( jpngr, pft )
      !******************************************************************************
      ! SUBROUTINE N_UPTAKE FOR FUN APPROACH
      !-------------------------------------------------------------------------
      ! This model calculates first the passive uptake of N via the transpiration stream
      ! then, if that uptake is insufficient to satisfy the demand from NPP, the remaining
      ! NPP is used to pay for N uptake by one or more of fixation and active uptake. The N
      ! extracted by this method must be the same as the N used
      ! for growth (n_grow=n_uptake) hence we calculate the optimum amount of carbon
      ! supplied to the roots/nodules to satisfy this constraint.
      ! Adopted from Fischer et al., 2010 by Beni Stocker, July 2012
      !-------------------------------------------------------------------------
      use classdefs
      use params_core
      use params_modl, only: lu_category, cton_pro, nfixer
      use statevars, only: dnpp, dtransp, dnup, dnup_pas, dnup_act, dnup_fix, dnup_ret, nind, ispresent, dcex
      use pools, only: nh4, no3, pexud, plabl, proot, wtot ! xxx try: wtot included here

      implicit none
      
      ! ARGUMENTS
      integer, intent(in) :: jpngr
      real :: soilmoist          ! is either dw1 or mw1, depending on temporal mode
      real :: soiltemp           ! is either dtemp_soil, or mtemp_soil, depending on temporal mode
      
      ! LOCAL VARIABLES
      integer, parameter :: nsl_eff = 1           ! effective number of soil layers for vertical distribution of N inorg
      integer, parameter :: icostActiveNH4 = 1    ! process number for active uptake of NH4
      integer, parameter :: icostActiveNO3 = 2    ! process number for active uptake of NO3
      integer, parameter :: icostFix = 3          ! process number for symbiotic N-fixation
      integer, parameter :: MAXSTEP = 100         ! maximum number of sub-time steps in N uptake optimisation
      integer, parameter :: PARTS = 100           ! number of parts of split up of inorganic N for uptake

      real, parameter :: SMALLVALUE = 1.e-9       ! to prevent numerical instability
      real, parameter :: BIG_COST = 999999.       ! dummy value

      logical unmetDemand
      logical ranked(nsl_eff)
      logical empty
      logical option(3)                           ! whether N-fixation is an option

      integer :: pft, iz, istep, lu   ! counter variables
      integer :: iminCost(nsl_eff)                ! info which source (-ID) is cheapest in each layer
      integer :: mloc(1)                          ! used to hold result of minloc
      integer :: cheapRank(nsl_eff)               ! ranking of layers: e.g. cheapRank(1) contains layer number with cheapest N
      integer :: icost
      integer :: nstep

      real :: nh4_nsl(nlu,nsl_eff)                ! vertically resolved NH4 pool (auxiliary variable!)
      real :: no3_nsl(nlu,nsl_eff)                ! vertically resolved NO3 pool (auxiliary variable!)
      real :: avail_nh4(nsl_eff)                  ! available NH4 in soil layer (gN/m2)
      real :: avail_no3(nsl_eff)                  ! available NO3 in soil layer (gN/m2)
      real :: nh4_conc(nsl_eff)                   ! NH4 concentration (gN/gH2O)
      real :: no3_conc(nsl_eff)
      real :: nh4_uptake(nsl_eff)                 ! (gN) 
      real :: no3_uptake(nsl_eff)                 ! (gN)
      real :: n_uptake_pass                       ! (gN)
      real :: npp_remaining
      real :: npp_remaining_step
      real :: Nacq_active
      real :: Nacq_fix
      real :: Cacq_active
      real :: Cacq_fix
      real :: n_demand
      real :: n_demand_remaining
      real :: n_uptake_retrans
      real :: epsilon_wtot       ! minimum water content [gH2O] for active uptake
      real :: cost(nsl_eff,3)    ! cost of N for each layer and source (-ID) [gC/gN]
      real :: cost_of_n(nsl_eff) ! cheapest cost of N for each layer [gC/gN]
      real :: Cacq                                           ! C spent for N acquisition
      real :: Nacq               ! N acquired by spending Cacq
      real :: excess
            
      ! xxx try
      soilmoist = 0.4
      soiltemp = 18.0
      wtot(:,:) = 12.
      

      ! LU category to which this PFT is associated
      lu = lu_category(pft)
      
      !-------------------------------------------------------------------------
      ! FUN treats different soil layers explicitly. In LPX N-pools represent a
      ! single soil N pool (100 cm). The code below is adopted for a
      ! formulation with different soil layers and soil N pools. To be consistent
      ! with (current) LPX, nsl_eff (effective number of soil layers) is
      ! set to 1. Thus, the model does not distinguish between N acquisition costs
      ! in different layers. An implementation of layer-specific N
      ! uptake, N-inorg pools, root turnover, N mineralization would be nice to have.
      !-------------------------------------------------------------------------


      !//////////////////////////////////////////////////////////////////////////
      ! INITIALIZATION
      !-------------------------------------------------------------------------
      if ( nsl_eff /= 1 ) then
        write(0,*) 'Think of something more clever wrt soil layering'
        stop
      else
        nh4_nsl(lu,:) = nh4(lu,jpngr)%n14
        no3_nsl(lu,:) = no3(lu,jpngr)%n14
      endif
      
      ! Absolute water content corresponding to permanent wilting point [gH2O].
      !epsilon_wtot = Fpwp( nlayers-nsl+1, jpngr ) * sum( Dz_soil(:) ) * 1000.
      !-------------------------------------------------------------------------

      ! LU-PFT association
      lu = lu_category(pft)

      ! xxx think about this: order in which PFTs get access to Ninorg matters!
      
      ! xxx try:
      epsilon_wtot = 0.05
      dtransp(pft) = wtot(1,1)*0.02

      ! ! Change the order in which pfts get access to N stores daily.
      ! do ppft=1,npft
      !   pft = modulo(ppft+day-2,npft)+1

      if (dnpp(pft)%c12 < 0.) then
        !-------------------------------------------------------------------------
        ! In case of negative daily NPP, no N uptake is necessary, no C is exuded
        ! npp has been calculated assuming no N limitation. On exit, npp has been reduced to account
        ! for expenditure on N uptake.
        !-------------------------------------------------------------------------
        dcex(pft) = 0.
       
      else

        !//////////////////////////////////////////////////////////////////////////
        ! INITIALIZATION
        !-------------------------------------------------------------------------
        unmetDemand    = .true.
        n_uptake_pass  = 0.
        nh4_uptake(:)  = 0.                           !sum of passive and active uptake. =0 if demand covered by N-stores
        no3_uptake(:)  = 0.                           !sum of passive and active uptake. =0 if demand co vered by N-stores

        Nacq_active    = 0.
        Nacq_fix       = 0.
        Cacq_active    = 0.
        Cacq_fix       = 0.
         
        no3_conc(:)    = 9999.e27                     ! assign dummy value before variable is passed to SR
        nh4_conc(:)    = 9999.e27                     ! assign dummy value before variable is passed to SR
        avail_no3(:)   = 9999.e27                     ! assign dummy value before variable is passed to SR
        avail_nh4(:)   = 9999.e27                     ! assign dummy value before variable is passed to SR
         
        option(:)      = .true. 

        !//////////////////////////////////////////////////////////////////////////
        ! N DEMAND
        !--------------------------------------------------------------------------
        ! As opposed to FUN, where N-demand is assumed to maintain the current
        ! CN-ratio of the whole plant, given NPP, here, N-demand is assumed is
        ! driven by the fixed CN-ratio of PFT-specific new production, given NPP.
        ! This is the same approach as in _DyN-LPJ (Xu-Ri & Prentice, 2008).
        ! Note, that N-fixers have higher CN-ratios of new production, implying a
        ! competitive disadvantage against non-fixers.
        !--------------------------------------------------------------------------
        n_demand = max( dnpp(pft)%c12/cton_pro(pft), 0. )  !in gN/m2/(timestep), timestetp=day or month 

        !tmp = dnpp(pft)

        !-------------------------------------------------------------------------
        ! Keep track of how much of the NPP-driven N-demand is not yet covered.
        !--------------------------------------------------------------------------
        n_demand_remaining = n_demand
         
        !//////////////////////////////////////////////////////////////////////////
        ! USE STORED N (RETRANSLOCATION)
        !--------------------------------------------------------------------------
        ! As opposed to original FUN model, in which N is retranslocated at a
        ! variable cost during leaf fall (turnover), a fraction of N is retained here
        ! from turnover. It is stored at the end of the last year and available to
        ! cover N demand during next year.
        !--------------------------------------------------------------------------
        n_uptake_retrans = min( n_demand, plabl(pft,jpngr)%n%n14 )
        
        ! just reduce the demand by amount retranslocated, not the labile N pool itself
        !plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - n_uptake_retrans

        !  if ( n_demand_remaining .gt. 0. .and. pft .eq. 3 ) then
        !  print*,'day, pft, frac. cov. by n_uptake_retrans', day, pft, n_uptake_retrans/n_demand_remaining
        !  endif

        n_demand_remaining = n_demand_remaining - n_uptake_retrans


        if ( n_demand_remaining > 0. ) then
          !//////////////////////////////////////////////////////////////////////////
          ! NO3 AND NH4 CONCENTRATION IN SOIL WATER
          !-------------------------------------------------------------------------
          ! N is assumed immobile and unavailable in dry soils. Using the total soil
          ! water and soil N to calculate concentration. Should be equivalent to
          ! assuming concentration of N as being equal in liquid and ice and assuming
          ! only N in liquid is available. nh4_conc and no3_conc are in gN/gH2O.
          !-------------------------------------------------------------------------
          do iz=1,nsl_eff

            call update_Ninorg( lu, nh4_nsl(lu,iz), epsilon_wtot, avail_nh4(iz), nh4_conc(iz) )
            call update_Ninorg( lu, no3_nsl(lu,iz), epsilon_wtot, avail_no3(iz), no3_conc(iz) )

          enddo
           
          !--------------------------------------------------------------------------
          ! Comment from FUN:
          ! Soil N is assumed to be uniform across the gridbox (we have a gridbox mean
          ! value), so all PFTs have access to the same resource. In theory this means
          ! that any PFT can access and extract all the N, even although the fractional
          ! cover of the PFT might be relatively small.
          ! Done this way, the earlier PFTs get preferential access to N (which, once
          ! extracted, will not be available to later PFTs). So as to avoid favouring
          ! these PFTs, the order in which the PFTs are considered changes between
          ! calls, using a counter of number of calls. (An alternative would be to use
          ! random number generator.)
          ! Note that if the sub-step loop was on the outside (and number of substeps>1),
          ! the order of the PFTs would be less important...but then all costs need to
          ! be in land_pts arrays, and generally takes more work to keep track of status
          ! within step loop.
          ! Calculations are carried out using fluxes for the fractional area in
          ! question, then converted to effect on gridbox mean by multiplying by frac.
          ! 
          ! I think, as long the time steps are small and N-demand(PFT)<<N-availability
          ! in any given time step, this should be not a problem. Giving a different
          ! PFT "first crack" at Ninorg in every time step would have to be taken
          ! care of in calling SR.
          !--------------------------------------------------------------------------
                 
          ! *** Why is passive uptake inside sub-loop in original FUN? ***

                 
          !//////////////////////////////////////////////////////////////////////////
          ! PASSIVE UPTAKE
          !--------------------------------------------------------------------------
          do iz=1,nsl_eff
             
            nh4_uptake(iz) = nh4_conc(iz) * dtransp(pft)
            no3_uptake(iz) = no3_conc(iz) * dtransp(pft)

            ! *** Associate cost of passive NO3 uptake due to pH equilibration by plant (Gutschik, 1981)? ***
             
            nh4_uptake(iz) = min( nh4_uptake(iz), avail_nh4(iz) )
            no3_uptake(iz) = min( no3_uptake(iz), avail_no3(iz) )
             
            n_uptake_pass = n_uptake_pass + nh4_uptake(iz) + no3_uptake(iz)
             
          enddo

          !--------------------------------------------------------------------------
          ! Passive uptake exceeds N-demand. Store excess in temporary storage.
          !--------------------------------------------------------------------------
          if ( n_demand_remaining <= n_uptake_pass ) then

            unmetDemand = .false.
            n_demand_remaining = 0.
            excess = n_uptake_pass - n_demand_remaining

            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 + excess
            excess = 0.

          endif

           
          do iz=1,nsl_eff
            !--------------------------------------------------------------------------
            ! Update N taken up and available N in soil.
            !--------------------------------------------------------------------------
            nh4_nsl(lu,iz) = nh4_nsl(lu,iz) - nh4_uptake(iz)
            no3_nsl(lu,iz) = no3_nsl(lu,iz) - no3_uptake(iz)
             
            call update_Ninorg( lu, nh4_nsl(lu,iz), epsilon_wtot, avail_nh4(iz), nh4_conc(iz) )
            call update_Ninorg( lu, no3_nsl(lu,iz), epsilon_wtot, avail_no3(iz), no3_conc(iz) )

            ! XXX May pools be emptied here?

          enddo

          if (unmetDemand) then
            !--------------------------------------------------------------------------
            ! If passive uptake is insufficient, consider fixation or active uptake
            !--------------------------------------------------------------------------

            !//////////////////////////////////////////////////////////////////////////
            ! COST ON SYMBIOTIC N-FIXATION
            !--------------------------------------------------------------------------
            ! This cost is independent of N availability and will thus not change in the
            ! course of the process of N uptake. Therefore, this cost is calculated out-
            ! side the sub-timestep loop, but updated every day in order to account for
            ! varying soil temperatures.
            !--------------------------------------------------------------------------
            if ( nfixer(pft) ) then
              do iz=1,nsl_eff
                cost(iz,icostFix) = fun_cost_fix( soiltemp )
              enddo
            else
              cost(:,icostFix) = BIG_COST
              option(icostFix) = .false.
            endif

            !--------------------------------------------------------------------------
            ! Comment from FUN:
            ! Start the sub-timestep integration loop. This loop was introduced
            ! as a conceptually simple way in which to account for the fact that
            ! the cost of N from a pool increases continuously as N is extracted
            ! from the pool, which can result in another pool becoming cheaper
            ! part way through the timestep. Each step in the loop attempts to
            ! satisfy part of the demand by extracting from the cheapest store.
            ! At the end of the step, the stores and costs are updated, for use
            ! in the following step.
            ! The number of sub-steps could be made a function of available store
            ! since a small extraction from a large store occurs at relatively
            ! constant cost, and it is probably unnecessary to use the sub-steps
            ! to update the cost. e.g. if demand<10% of store, use a single step.
            ! Also...if passive uptake is sufficient to meet demand, we could again
            ! use a single step...but probably rarely occurs.
            ! FOR ME NOT LOGICAL WHY PASSIVE UPTAKE IS INSIDE SUB-LOOP IN ORIGINAL
            ! CODE.
            !--------------------------------------------------------------------------

            !--------------------------------------------------------------------------
            ! Calculate the number of sub-steps so that the N demand to be covered per step
            ! (n_demand_remaining/nstep) is equal to a fraction of 1/PARTS of the inorganic
            ! N pool (taking the larger of NO3 or NH4):
            !          Ninorg/PARTS = n_demand_remaining/nstep
            ! => solve for nstep...
            !--------------------------------------------------------------------------
            !print*,'MAX(no3(lu,jpngr)%n14,nh4(lu,jpngr)%n14))',MAX(no3(lu,jpngr)%n14,nh4(lu,jpngr)%n14)
            !nstep = MIN(INT(n_demand_remaining*real(PARTS)/MAX(no3(lu,jpngr)%n14,nh4(lu,jpngr)%n14)),MAXSTEP)
            nstep = 100
            
            !--------------------------------------------------------------------------
            ! Calculate NPP still to pay for after passive uptake and retranslocation
            !--------------------------------------------------------------------------
            npp_remaining = n_demand_remaining*cton_pro(pft)
             
            !--------------------------------------------------------------------------
            ! Split NPP to pay for into steps
            !--------------------------------------------------------------------------
            npp_remaining_step = npp_remaining/real(nstep)


            !//////////////////////////////////////////////////////////////////////////
            ! SUB-TIMESTEP LOOP STARTS HERE
            !--------------------------------------------------------------------------
            do istep=1,nstep

              !--------------------------------------------------------------------------
              ! Initializations
              !--------------------------------------------------------------------------
              Nacq = 0.
              Cacq = 0.
              empty = .false.

              !//////////////////////////////////////////////////////////////////////////
              ! COST OPTIMIZATION
              !--------------------------------------------------------------------------
              ! Calculate the costs of the various possible sources and identify the
              ! cheapest. The cost of nitrogen fixation does not vary within a timestep,
              ! and has been calculated outside this loop.
              !--------------------------------------------------------------------------
               do iz=1,nsl_eff
                 
                !--------------------------------------------------------------------------
                ! Cost of active uptake of NH4, NO3 and organic N
                !--------------------------------------------------------------------------
                cost(iz,icostActiveNH4) = fun_cost_active_nh4( avail_nh4(iz), proot(pft,jpngr)%c%c12*nind(pft,jpngr) )
                cost(iz,icostActiveNO3) = fun_cost_active_no3( avail_no3(iz), proot(pft,jpngr)%c%c12*nind(pft,jpngr) )
                ! cost(iz,icostActiveOrg) = fun_cost_active_org(litterN)
                 
                !--------------------------------------------------------------------------
                ! Find the lowest cost in each soil layer and each sub-timestep
                !--------------------------------------------------------------------------
                mloc(1:1) =  MINLOC( cost(iz,:), mask=option(:) )
                iminCost(iz) = mloc(1)
                cost_of_n(iz)= cost(iz,iminCost(iz))
                 
              enddo
               
              !--------------------------------------------------------------------------
              ! Rank layers in terms of cost of N uptake, cheapest first
              !--------------------------------------------------------------------------
              ranked(:) = .false.
              do iz=1,nsl_eff
                mloc(1:1) = MINLOC(cost_of_n(:),.not.ranked(:))
                cheapRank(iz) = mloc(1)
                ranked(mloc(1)) = .true.
              enddo
               
              !--------------------------------------------------------------------------
              ! Spend C on extracting N
              !--------------------------------------------------------------------------

              ! Use the layer with the cheapest N
              iz = cheapRank(1)

              if ( cost_of_n(iz) - BIG_COST > SMALLVALUE ) then
                ! If the cheapest N is very expensive, stores are empty                 
                empty = .true.
                 
              else
                ! Identify the process involved with the cheapest N
                icost = iminCost(iz)
                 
                !--------------------------------------------------------------------------
                ! Calculate NPP decrement and N increment (for derivation, see bottom of this
                ! file). Cacq: C spent for N acquisition. Nacq: N acquired
                !--------------------------------------------------------------------------
                Cacq = npp_remaining_step/(cton_pro(pft)/cost(iz,icost)+1.)
                Nacq = Cacq/cost(iz,icost)
                 
                ! print*,'pft, cost(iz,icost)', pft, cost(iz,icost)

                !--------------------------------------------------------------------------
                ! Load values into process-specific variables and limit
                ! increments if store exhausted. Note that the fluxes and
                ! increments calculated so far are for this fraction, and need
                ! to be scaled to get effect on gridbox scale.
                !--------------------------------------------------------------------------

                if (icost==icostActiveNH4) then
                  !//////////////////////////////////////////////////////////////////////////
                  ! ACTIVE NH4 UPTAKE
                  !--------------------------------------------------------------------------
                  ! Decrement available store and update concentration.
                  ! If store insufficient, reduce uptake and NPP increment.
                  !--------------------------------------------------------------------------
                  if ( Nacq > avail_nh4(iz) ) then
                    Nacq = avail_nh4(iz)
                    Cacq = Nacq*cost(iz,icost)
                  endif

                  nh4_nsl(lu,iz) = nh4_nsl(lu,iz) - Nacq

                  call update_Ninorg( lu, nh4_nsl(lu,iz), epsilon_wtot, avail_nh4(iz), nh4_conc(iz) )

                  Nacq_active = Nacq_active + Nacq
                  Cacq_active = Cacq_active + Cacq
                  nh4_uptake(iz) = nh4_uptake(iz) + Nacq
                   
                elseif (icost==icostActiveNO3) then
                  !//////////////////////////////////////////////////////////////////////////
                  ! ACTIVE NO3 UPTAKE
                  !--------------------------------------------------------------------------
                  ! Decrement available store and update concentration.
                  ! If store insufficient, reduce uptake and NPP increment.
                  !--------------------------------------------------------------------------
                  if ( Nacq > avail_no3(iz) ) then
                    Nacq = avail_no3(iz)
                    Cacq = Nacq*cost(iz,icost)
                  endif

                  no3_nsl(lu,iz) = no3_nsl(lu,iz) - Nacq

                  call update_Ninorg( lu, no3_nsl(lu,iz), epsilon_wtot, avail_no3(iz), no3_conc(iz) )

                  Nacq_active = Nacq_active + Nacq
                  Cacq_active = Cacq_active + Cacq
                  no3_uptake(iz) = no3_uptake(iz) + Nacq
                   
                  ! elseif ( icost == icostActiveOrg ) then
                  !//////////////////////////////////////////////////////////////////////////
                  ! ACTIVE UPTAKE OF ORGANIC N
                  !--------------------------------------------------------------------------
                  ! Does it deplete any explecit store?
                  ! XXX To be implemented XXX

                else

                  !//////////////////////////////////////////////////////////////////////////
                  ! SYMBIOTIC N FIXATION
                  !--------------------------------------------------------------------------
                  if ( .not. nfixer(pft) ) then
                    write(0,*) 'cost:', cost(:,:)
                    write(0,*) 'Non-N-fixer tries to fix N! '//'BIG_COST too small?'
                    stop
                  endif
                   
                  Nacq_fix = Nacq_fix + Nacq
                  Cacq_fix = Cacq_fix + Cacq
                   
                  ! XXX Note that fixing will remain the cheapest store for all
                  ! subsequent steps (effectively it is accessing an infinite
                  ! store of gaseous N2), so explicit integration over each of
                  ! those individual steps could be avoided and NPP summed over
                  ! all remaining steps. But this has not been coded.

                endif
                 
                !--------------------------------------------------------------------------
                ! Decrement unmet NPP over sub-timestep. This is necessary to know
                ! how much is exuded (the remaing npp) when pools get empty (during
                ! loop.)
                !--------------------------------------------------------------------------
                npp_remaining = npp_remaining - Cacq
                 
                !--------------------------------------------------------------------------
                ! Update costs in order to update ranking of soil layers w.r.t. cost 
                !--------------------------------------------------------------------------
                if (istep < nstep) then
                   
                  if (icost==icostActiveNH4) then

                    cost(iz,icostActiveNH4) = fun_cost_active_nh4( avail_nh4(iz), proot(pft,jpngr)%c%c12*nind(pft,jpngr) )

                  elseif (icost==icostActiveNO3) then

                    cost(iz,icostActiveNO3) = fun_cost_active_no3( avail_no3(iz), proot(pft,jpngr)%c%c12*nind(pft,jpngr) )
                     
                    ! If N-org uptake depletes stores of N-org, update costs here
                    ! elseif (icost==icostActiveOrg) then

                    ! cost(iz,icostActiveOrg) = fun_cost_active_org(litterN)

                  endif

                  !--------------------------------------------------------------------------
                  ! Decide which source of N is now cheapest in this layer. If cheapest layer
                  ! is now not the same as the cheapest layer before, re-rank layers.
                  !--------------------------------------------------------------------------
                  mloc(1:1) = MINLOC( cost(iz,:) )
                  iminCost(iz) = mloc(1)
                  cost_of_n(iz) = cost(iz,iminCost(iz))
                   
                  ! Only used when nsl_eff > 1. 
                  !                if ( nsl_eff .gt. 1 ) then
                  !                  if ( cost_of_n(iz) .gt. cost_of_n(cheapRank(2)) ) then
                  ! ! We only have to insert a single value...but sort whole list!
                  !                    ranked(:) = .false.
                  !                    do i=1,nsl_eff
                  !                      mloc(1:1) = MINLOC( cost_of_n(:), .not. ranked(:) )
                  !                      cheapRank(i) = mloc(1)
                  !                      ranked(mloc(1)) = .true.
                  !                    enddo
                  !                  endif
                  !                endif
                   
                endif                                    !  istep<nstep
                 
              endif                                      !  cost>big_cost (empty)
               
              !//////////////////////////////////////////////////////////////////////////
              ! UPDATE FLUXES
              !--------------------------------------------------------------------------
              ! If pools got empty during integration loop, that is, before n_demand has
              ! fully been covered, 'empty' is true. In this case, all remaining npp is
              ! exuded. And the integration loop is exited.
              !--------------------------------------------------------------------------
              if ( empty ) then
                 
                ! Do not reduce "real" npp here, but substract exu in npp-SR
                ! npp(pft,1) = npp(pft,1) - npp_remaining
                dcex(pft) = dcex(pft) + npp_remaining

                ! print*,'everything exuded'

                ! Leave integration loop (steps)
                exit
                 
              else                                       !not empty

                !--------------------------------------------------------------------------
                ! Decrement unmatched NPP over sub-timestep. This is necessary to know
                ! how much is exuded (the remaing npp) when pools get empty (during
                ! loop.) npp_remaining is the NPP still to pay for. Cacq is the cost of
                ! N-acquisition in this integration step. 
                !--------------------------------------------------------------------------
                npp_remaining = npp_remaining - Cacq
                 
                ! Do not reduce "real" npp here, but substract exu in npp-SR
                ! npp(pft,1) = npp(pft,1)-Cacq

                dcex(pft) = dcex(pft) + Cacq

              endif                                      !empty
               
            enddo                                        !steps
             
            !  print*,'fraction of npp paid for N uptake',
            !  $             (Cacq_fix+Cacq_active)/tmp

            !--------------------------------------------------------------------------
            ! Update annual N-uptake of this PFT. N-retranslocation is not considered
            ! N-uptake.
            !--------------------------------------------------------------------------
            dnup(pft)%n14 = dnup(pft)%n14 + n_uptake_pass  + Nacq_active + Nacq_fix  ! n_uptake_retrans is not considered uptake
            dnup_pas(pft) = dnup_pas(pft) + n_uptake_pass
            dnup_act(pft) = dnup_act(pft) + Nacq_active                   
            dnup_fix(pft) = dnup_fix(pft) + Nacq_fix  
            dnup_ret(pft) = dnup_ret(pft) + n_uptake_retrans                 


            ! ! XXX debug
            !          if (pft.eq.1) then
            ! !                 print*,'N uptake',(n_uptake_pass+Nacq_active+Nacq_fix+n_uptake_retrans)
            ! !                 print*,'Cgrowth',(npp(pft)-exu(pft))
            ! !                 print*,'N uptake by BNF',Nacq_fix
            ! !                 print*,'N uptake by act',Nacq_active
            !            print*,'actual CN ratio',
            !$                (npp(pft)-exu(pft))
            !$                /(n_uptake_pass+Nacq_active+Nacq_fix+n_uptake_retrans)
            ! !                 if (Nacq_fix.gt.0.) print*,'daily cost of BNF', Cacq_fix/Nacq_fix
            ! !                 if (Nacq_active.gt.0.) print*,'daily cost of act', Cacq_active/Nacq_active
            !          endif

               
          endif                                          !unmetDemand
           
           
        else                                             !(n_demand_remaining.lt.0.)

        !  XXX no n_uptake at cost > 0, n-store can cover demand  
        !  any executions?
           
        endif                                            !(n_demand_remaining.lt.0.)

      endif                                              ! (npp(pft).lt.0.)
          
      ! Update common-variable-Ninorg-pools again
      nh4(lu,jpngr)%n14 = sum(nh4_nsl(lu,:))
      no3(lu,jpngr)%n14 = sum(no3_nsl(lu,:))


      contains
      !//////////////////////////////////////////////////////////////////////////////

      subroutine update_Ninorg( lu, Ninorg, epsilon_wtot, avail_Ninorg, Ninorg_conc )
      !******************************************************************************
      ! Update N taken up and available N in soil.
      !--------------------------------------------------------------------------
      use pools, only: wtot
      implicit none

      integer, intent(in) :: lu
      real, intent(in)  :: Ninorg
      real, intent(in)  :: epsilon_wtot
      real, intent(out) :: avail_Ninorg
      real, intent(out) :: Ninorg_conc
                    
      Ninorg_conc = Ninorg/wtot(lu,jpngr)
      
      if ( wtot(lu,jpngr) > epsilon_wtot ) then 
        avail_Ninorg = Ninorg - epsilon_wtot*Ninorg_conc
      else
        avail_Ninorg = 0.
        Ninorg_conc  = 0.
      endif
    
      end subroutine update_Ninorg


      function fun_cost_active_nh4( avail_nh4, croot )
      !******************************************************************************
      ! Cost of active NH4 uptake by exudation of labile C. From Fisher er al., 2010.
      !--------------------------------------------------------------------------      
      use params_modl
      implicit none

      real :: fun_cost_active_nh4         ! gC/gN
      real, intent(in) :: avail_nh4
      real, intent(in) :: croot

      if ( croot > SMALLVALUE ) then
        if ( avail_nh4 > SMALLVALUE ) then
          fun_cost_active_nh4 = KN_ACTIVE_NH4/avail_nh4 * KC_ACTIVE_NH4/croot
        else
          fun_cost_active_nh4 = BIG_COST
        endif

      else
        fun_cost_active_nh4 = BIG_COST
      endif

      end function fun_cost_active_nh4


      function fun_cost_active_no3( avail_no3, croot )
      !******************************************************************************
      ! Cost of active NO3 uptake by exudation of labile C. From Fisher er al., 2010.
      !--------------------------------------------------------------------------      
      use params_modl
      implicit none

      real :: fun_cost_active_no3         ! gC/gN
      real, intent(in) :: avail_no3
      real, intent(in) :: croot

      if (croot > SMALLVALUE ) then
        if ( avail_no3 > SMALLVALUE ) then
          fun_cost_active_no3 = KN_ACTIVE_NO3/avail_no3 * KC_ACTIVE_NO3/croot ! XXX +x [gC/gN for pH equilibration]
        else
          fun_cost_active_no3 = BIG_COST
        endif

      else
        fun_cost_active_no3 = BIG_COST
      endif
      
      end function fun_cost_active_no3


      !******************************************************************************

      !     function fun_cost_active_org(litterN)
      !     implicit none
      !     include XXX
      !     real*8 litterN


      !     return
      !     end

      function fun_cost_fix( soiltemp )
      !******************************************************************************
      ! Cost of symbiotic N fixation is the inverse of nitrogenase activity
      ! after Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
      ! (value from Gutschik 1981)
      !--------------------------------------------------------------------------      
      use params_modl
      implicit none

      real, intent(in) :: soiltemp
      real :: fun_cost_fix                 ! function return variable


      fun_cost_fix = MINIMUMCOSTFIX + exp((soiltemp-FIXOPTIMUM)**2/(2*FIXWIDTH**2))    ! inverse gauss function  (take WARMEST layer)

      end function fun_cost_fix
      

      !******************************************************************************
      ! Derivation of Cacq (C spent to cover cost of N-uptake) after
      ! Fisher et al., 2010 (Equation numbers from paper)
      ! 
      !    C_growth = C_npp - C_acq                (eq.6b)
      !    N_acq    = C_acq / Cost_acq             (eq.6c)
      !    r_cton   = C_growth / (N_passive+N_acq) (eq.6d)  [equation presented in paper is incorrect!]

      ! Using 6b and 6c, eq.6d becomes
      !    r_cton   = (C_npp - C_acq) / (N_passive + C_acq/Cost_acq)

      ! Solving for C_acq yields
      !    C_acq    = (C_npp - r_cton * N_pass)/(r_cton/Cost_acq + 1)

      ! Identify terms with variables in code:
      ! (C_npp - r_cton * N_pass) <=> npp_remaining_step
      ! C_acq <=> Cacq
      ! N_acq <=> Nacq   [rest is obvious]
      ! 
      !******************************************************************************

      ! REFERENCES
      ! Fisher 
      ! Gutschik
      ! Houlton

      end subroutine nuptake

      
c$$$      subroutine n_fixation_cryptogam( day, lu, jpngr, dnfix_cpc, dnfix_cgc )
c$$$      !******************************************************************************
c$$$      ! SUBROUTINE N_UPTAKE BY CRYPTOGAMIC COVERS
c$$$      !-------------------------------------------------------------------------
c$$$      ! Simulated to match pattern and global total fixed N after Elbert et al.
c$$$      ! (2012), Nature Geoscience. Basic assumption: N uptake is driven by energy
c$$$      ! available (solar radiation ~ photosynthetically active radiation) and not
c$$$      ! absorbed by leafs or stems. N fixation by cryptogamic ground cover (CGC)
c$$$      ! thus scales with (1-VPC), where VPC is analogous to FPC but takes into
c$$$      ! account the shading by branches and stems. N fixation by cryptogamic
c$$$      ! plant covers (CPC) scales with SPC. 
c$$$      !-------------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$      ! ARGUMENTS
c$$$      INTEGER day, lu, jpngr
c$$$      REAL*8 dnfix_cpc, dnfix_cgc
c$$$      
c$$$      ! LOCAL VARIABLES
c$$$      INTEGER
c$$$     $     pft,ppft
c$$$      
c$$$      REAL*8
c$$$     $     fpc_ind,               ! phenology-modulated (!) fractional plant cover
c$$$     $     local_fpc_grid,        ! FPC w.r.t. grid cell area (is not the same as the global variable fpc_grid)
c$$$     $     vpc_ind,               ! fractional vegetation cover including stems and branches
c$$$     $     vpc_grid,              ! VPC w.r.t. grid cell area
c$$$     $     spc_grid,              ! fractional stem/branches cover
c$$$     $     fpc_grid_total,        ! fpc_grid summed over all PFTs in the present LU
c$$$     $     vpc_grid_total,        ! vpc_grid summed over all PFTs in the present LU
c$$$     $     spc_grid_total,        ! spc_grid summed over all PFTs in the present LU
c$$$     $     lm_tot(npft),
c$$$     $     scale
c$$$
c$$$      ! Initialisations
c$$$      vpc_grid_total = 0.
c$$$      fpc_grid_total = 0.
c$$$      spc_grid_total = 0.
c$$$
c$$$      !  ! Calculate ftemp
c$$$      !  if (soiltemp.ge.-40.) then
c$$$      !    tshift = 46.02d0
c$$$      !    ftemp = exp(308.56d0*(1.0/(20.+tshift)-1.0/
c$$$      ! $       (soiltemp+tshift)))                             ! Eq.8, XP08 (canexch.cpp:1018)
c$$$      !  else
c$$$      !    ftemp = 0.
c$$$      !  endif
c$$$      !  ftemp = min(ftemp, 1.)                              ! (canexch.cpp:1023)
c$$$      !  ftemp = max(ftemp, 0.)                              ! (canexch.cpp:1024)      
c$$$
c$$$      do pft=1,npft
c$$$        if ( present(pft,jpngr) .and. lu_category(pft) .eq. lu ) then
c$$$
c$$$        ! LM_TOT
c$$$        !--------------------------------------------------------------------------
c$$$        ! Non-linearity of Beer-Law causes very high FPC values when 2 Grasses are present.
c$$$        ! (Beer Law does NOT make sense for grasses, anyway.)
c$$$        ! Thus, use sum of all grass/moss-leaf masses and calculate FPC based on the sum.
c$$$        ! Then compute each PFT's FPC as the product of total-grass FPC times each PFT's leaf mass.
c$$$        !-------------------------------------------------------------------------
c$$$          lm_tot(pft) = 0.
c$$$          ! Grass: C3, C4 on natural, croplands, pasture, peatlands
c$$$          if (grass(pft)) then
c$$$            do ppft=1,npft
c$$$              if (lu_category(ppft).eq.lu_category(pft)) then
c$$$                if (grass(ppft)) lm_tot(pft) =
c$$$     $               lm_tot(pft)+lm_ind(ppft,jpngr,1)
c$$$              endif
c$$$            enddo
c$$$          ! Moss: moss on peatlands
c$$$          elseif (moss(pft)) then
c$$$            do ppft=1,npft
c$$$              if (lu_category(ppft).eq.lu_category(pft)) then
c$$$                if (moss(ppft)) lm_tot(pft) =
c$$$     $               lm_tot(pft)+lm_ind(ppft,jpngr,1)
c$$$              endif
c$$$            enddo
c$$$          ! Tree: tree on natural lands, peatlands
c$$$          else
c$$$            lm_tot(pft) = lm_ind(pft,jpngr,1)
c$$$          endif
c$$$          
c$$$          ! LAI
c$$$          !--------------------------------------------------------------------------
c$$$          if (crownarea(pft,jpngr).gt.0.) then
c$$$            lai_ind(pft,jpngr)=(lm_tot(pft)*sla(pft))/
c$$$     $           crownarea(pft,jpngr)
c$$$          else
c$$$            lai_ind(pft,jpngr)=0.
c$$$          endif
c$$$          
c$$$          ! FPC and VPC
c$$$          !--------------------------------------------------------------------------
c$$$          ! Note that this is not identical to how it's calculated in SR update_fpc,
c$$$          ! where the phenology scaling factor is not included in the exponent.
c$$$          ! Fractional plant cover accounts for the fraction of the grid cell covered
c$$$          ! by the photosyntesic plant tissue. To be modulated by daily phenology!
c$$$          !--------------------------------------------------------------------------
c$$$          fpc_ind = 1.-dexp(
c$$$     $                        -1.*kbeer*lai_ind(pft,jpngr)*dphen(day,pft)
c$$$     $                        )
c$$$          vpc_ind = 1.-dexp(
c$$$     $                        -1.*kbeer*(
c$$$     $                                      lai_ind(pft,jpngr)*dphen(day,pft)
c$$$     $                                      + pftpar(pft,46)
c$$$     $                                      )
c$$$     $                        )
c$$$          
c$$$          local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
c$$$          vpc_grid       = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
c$$$          
c$$$          if (lm_tot(pft).gt.0.) then
c$$$            local_fpc_grid = local_fpc_grid*lm_ind(pft,jpngr,1)
c$$$     $           /lm_tot(pft)
c$$$            vpc_grid = vpc_grid*lm_ind(pft,jpngr,1)/lm_tot(pft)
c$$$          else
c$$$            local_fpc_grid = 0.
c$$$            vpc_grid       = 0. 
c$$$          endif
c$$$
c$$$          spc_grid = vpc_grid - local_fpc_grid
c$$$
c$$$          ! Sum over pfts
c$$$          !--------------------------------------------------------------------------
c$$$          fpc_grid_total = fpc_grid_total + local_fpc_grid
c$$$          vpc_grid_total = vpc_grid_total + vpc_grid
c$$$          spc_grid_total = spc_grid_total + spc_grid
c$$$
c$$$          ! print*,'spc_grid',spc_grid
c$$$          
c$$$          !!          call update_fpc(pft,jpngr)
c$$$          !          
c$$$          !      ! VAI is analogous to LAI but accounts for stems and branches in addition to
c$$$          !      ! leafs.
c$$$          !          vpc_ind = 1. - dexp(
c$$$          !     $                          - 1.*kbeer*(
c$$$          !     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
c$$$          !     $                                         + pftpar(pft,46)
c$$$          !     $                                         )
c$$$          !     $                          )
c$$$          !          vpc_grid = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
c$$$          !          vpc_grid_total = vpc_grid_total + vpc_grid
c$$$          !
c$$$          !      ! Calculate local FCP treating dphen analogously as for the calulation of VAI:
c$$$          !      ! FPC = 1-exp(-kbeer*LAI*dphen) instead of FPC = dphen*(1-exp(-kbeer*LAI))
c$$$          !!           fpc_ind = 1. - dexp(
c$$$          !!     $                           -1.*kbeer*(
c$$$          !!     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
c$$$          !!     $                                         )
c$$$          !!     $                           )
c$$$          !          fpc_ind = (1. - dexp(
c$$$          !     $                           -1.*kbeer*(
c$$$          !     $                                         lai_ind(pft,jpngr)
c$$$          !     $                                         )
c$$$          !     $                           ))!*dphen(day,pft)
c$$$          !          local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
c$$$          !          fpc_grid_total = fpc_grid_total + local_fpc_grid
c$$$          !
c$$$          !          print*,'pft',pft
c$$$          !          print*,'local_fpc_grid     ',local_fpc_grid
c$$$          !          print*,'fpc_grid(pft,jpngr)',fpc_grid(pft,jpngr)
c$$$          !          
c$$$          !      ! Calculate fractional stem/branch cover of grid cell as the difference
c$$$          !          spc_grid = vpc_grid - local_fpc_grid
c$$$          !          spc_grid_total = spc_grid_total + spc_grid
c$$$         
c$$$        endif
c$$$      enddo
c$$$
c$$$      
c$$$      if (vpc_grid_total.gt.1.) then
c$$$        !        print*,'-----------------scaling-------------------'
c$$$        !        print*,'fpc_grid_total',fpc_grid_total
c$$$        !        print*,'vpc_grid_total',vpc_grid_total
c$$$        scale = 1. / vpc_grid_total
c$$$        fpc_grid_total = fpc_grid_total * scale
c$$$        vpc_grid_total = vpc_grid_total * scale
c$$$        spc_grid_total = spc_grid_total * scale
c$$$        !        print*,'fpc_grid_total',fpc_grid_total
c$$$        !        print*,'vpc_grid_total',vpc_grid_total
c$$$      endif
c$$$
c$$$      if (fpc_grid_total.gt.1.) then
c$$$        !        print*,'fpc_grid_total',fpc_grid_total
c$$$        stop 
c$$$      endif
c$$$
c$$$      ! Daily N fixed by cryptogamic ground and plant covers (communicated to calling SR)
c$$$      !-------------------------------------------------------------------------
c$$$      ! Fixation scales with daily photosynthetically active radiation and the
c$$$      ! branch/stem surface for CPC and the bare ground surface for CGC.
c$$$      
c$$$      dnfix_cpc = par_day(day) * max( 0., spc_grid_total) / glob_CPC_scal
c$$$      dnfix_cgc = par_day(day) * max( 0., (1.0 - vpc_grid_total) ) / glob_CGC_scal
c$$$
c$$$      end subroutine n_fixation_cryptogam
      
 
