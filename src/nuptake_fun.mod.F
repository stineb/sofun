module nuptake_fun
  !////////////////////////////////////////////////////////////////
  ! FUN NITROGEN UPTAKE MODULE
  ! Contains the "main" subroutine 'nuptake' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'waterbal' must contain this list 
  ! of subroutines (names that way).
  !   - nuptake
  !   - getpar_modl_nuptake
  !   - initio_nuptake
  !   - initoutput_nuptake
  !   - getout_daily_nuptake
  !   - getout_monthly_nuptake
  !   - writeout_ascii_nuptake
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - daily NPP ('dnpp')
  !   - soil temperature ('xxx')
  !   - inorganic N pools ('no3', 'nh4')
  !   - xxx 
  !----------------------------------------------------------------
  use params_core, only: ndayyear, nmonth, nlu, npft, maxgrid

  implicit none

  ! FUN PARAMETERS
  real :: MINIMUMCOSTFIX    ! minimum cost of N fixation (at optimal temperature)
  real :: FIXOPTIMUM        ! optimum temperature for N fixation
  real :: FIXWIDTH          ! shape parameter for width of N fixation cost function
  real :: KN_ACTIVE_NO3     ! N-availability constant in cost function for active NO3 uptake
  real :: KN_ACTIVE_NH4     ! N-availability constant in cost function for active NH4 uptake
  real :: KC_ACTIVE_NO3     ! root exploration constant in cost function for active NO3 uptake
  real :: KC_ACTIVE_NH4     ! root exploration constant in cost function for active NH4 uptake
  real :: EPSILON_WTOT      ! minimum soil water content above which N becomes available  

  ! MODULE-SPECIFIC VARIABLES
  ! These are not required outside module STASH, but are used in different SRs of this module
  real, dimension(npft) :: dccost           ! daily mean C cost of N uptake [gC/gN] 
  real, dimension(npft) :: dnup_pas         ! daily passive N uptake [gN/m2/d]
  real, dimension(npft) :: dnup_act         ! daily active N uptake [gN/m2/d]  
  real, dimension(npft) :: dnup_fix         ! daily N uptake by plant symbiotic N fixation [gN/m2/d]
  real, dimension(npft) :: dnup_ret         ! daily N uptake [gN/m2/d]


  ! FUN OUTPUT VARIABLES
  real, dimension(npft,ndayyear,maxgrid) :: outdccost   ! daily mean C cost of N uptake (gC/gN) 
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_pas
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_act
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_fix
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_ret


contains

  subroutine initdaily_nuptake
    !////////////////////////////////////////////////////////////////
    ! Initialise daily variables with zero
    !----------------------------------------------------------------
    implicit none

    dnup_pas(:)    = 0.0
    dnup_act(:)    = 0.0
    dnup_fix(:)    = 0.0
    dnup_ret(:)    = 0.0

    return

  end subroutine initdaily_nuptake


  subroutine initio_nuptake( prefix )
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    character(len=*) :: prefix

    ! LOCAL VARIABLES
    character(len=256) :: filnam


    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------

    ! MEAN DAILY C COST OF N UPTAKE (gC/gN)
    filnam=trim(prefix)//'.d.ccost.out'
    open(400,file=filnam,err=888,status='unknown')

    ! PASSIVE N UPTAKE (gN)
    filnam=trim(prefix)//'.d.nup_pas.out'
    open(401,file=filnam,err=888,status='unknown')

    ! ACTIVE N UPTAKE (gN)
    filnam=trim(prefix)//'.d.nup_act.out'
    open(402,file=filnam,err=888,status='unknown')

    ! SYMBIOTIC BNF (gN)
    filnam=trim(prefix)//'.d.nup_fix.out'
    open(403,file=filnam,err=888,status='unknown')

    ! RETRANSLOCATED N FROM LABILE POOL TO SATISFY DEMAND (gN)
    filnam=trim(prefix)//'.d.nup_ret.out'
    open(404,file=filnam,err=888,status='unknown')

    return

    888  stop 'INITIO_NUPTAKE: error opening output files'

  end subroutine initio_nuptake



  subroutine getpar_modl_nuptake( paramfilnam )
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    character(len=*) paramfilnam

    ! LOCAL VARIABLES
    integer, parameter :: npar = 8
    real, dimension(npar) :: params_array


    open(unit=01,file=trim(paramfilnam),status='OLD')      
    read (01,*) params_array
    close (01)

    ! shape parameter of cost function of N fixation 
    ! Below parameters (MINIMUMCOSTFIX, FIXOPTIMUM, FIXWIDTH ) are based on 
    ! the assumption that the cost of symbiotic N fixation is the 
    ! inverse of nitrogenase activity. 
    ! After Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    MINIMUMCOSTFIX = params_array(1)

    ! shape parameter of cost function of N fixation 
    FIXOPTIMUM = params_array(2)

    ! shape parameter of cost function of N fixation 
    FIXWIDTH = params_array(3)

    ! N-availability constant in cost function for active NO3 uptake
    KN_ACTIVE_NO3 = params_array(4)
    
    ! N-availability constant in cost function for active NH4 uptake
    KN_ACTIVE_NH4 = params_array(5)

    ! root exploration constant in cost function for active NO3 uptake
    KC_ACTIVE_NO3 = params_array(6)
    
    ! root exploration constant in cost function for active NH4 uptake
    KC_ACTIVE_NH4 = params_array(7)

    ! minimum soil water content above which N becomes available 
    EPSILON_WTOT = params_array(8)

    return

  end subroutine getpar_modl_nuptake



  subroutine initoutput_nuptake
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    implicit none

    ! xxx remove their day-dimension
    outdccost(:,:,:) = 0.0

    return

  end subroutine initoutput_nuptake



  subroutine getout_daily_nuptake( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    ! Save the daily totals:
    ! xxx add lu-dimension and jpngr-dimension
    outdccost(:,doy,jpngr) = dccost(:)
    outdnup_pas(:,doy,jpngr) = dnup_pas(:)
    outdnup_act(:,doy,jpngr) = dnup_act(:)
    outdnup_fix(:,doy,jpngr) = dnup_fix(:)
    outdnup_ret(:,doy,jpngr) = dnup_ret(:)

    return  

  end subroutine getout_daily_nuptake


  subroutine writeout_ascii_nuptake( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use params_core, only: ndayyear, nmonth, npft
    use params_siml, only: firstyeartrend, spinupyears

    implicit none

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

        if (nlu>1) stop 'writeout_ascii_nuptake: write out lu-area weighted sum'
        if (npft>1) stop 'writeout_ascii_nuptake: think of something for ccost output'

        ! xxx lu-area weighted sum if nlu>0
        write(400,999) itime, sum(outdccost(:,day,jpngr))/npft
        write(401,999) itime, sum(outdnup_pas(:,day,jpngr))
        write(402,999) itime, sum(outdnup_act(:,day,jpngr))
        write(403,999) itime, sum(outdnup_fix(:,day,jpngr))
        write(404,999) itime, sum(outdnup_ret(:,day,jpngr))

      end do
    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_nuptake


  subroutine nuptake( jpngr, pft )
    !/////////////////////////////////////////////////////////////////
    ! SUBROUTINE N_UPTAKE FOR FUN APPROACH
    !-----------------------------------------------------------------
    ! This model calculates first the passive uptake of N via the transpiration stream
    ! then, if that uptake is insufficient to satisfy the demand from NPP, the remaining
    ! NPP is used to pay for N uptake by one or more of fixation and active uptake. The N
    ! extracted by this method must be the same as the N used
    ! for growth (n_grow=n_uptake) hence we calculate the optimum amount of carbon
    ! supplied to the roots/nodules to satisfy this constraint.
    ! Adopted from Fischer et al., 2010 by Beni Stocker, July 2012
    !-----------------------------------------------------------------
    use classdefs
    use params_core
    use params_modl, only: lu_category, cton_pro, nfixer
    use fluxes, only: dnpp, dnup, dcex, dwn, aet_d
    use treegeometry, only: nind, ispresent, fpc_grid
    use pools, only: nh4, no3, pexud, plabl, proot 

    implicit none
    
    ! ARGUMENTS
    integer, intent(in) :: jpngr, pft

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
    logical ranked         
    logical empty
    logical option(3)                           ! whether N-fixation is an option

    integer :: iz, istep, lu   ! counter variables
    integer :: iminCost                         ! info which source (-ID) is cheapest in each layer
    integer :: mloc(1)                          ! used to hold result of minloc
    integer :: cheapRank                        ! ranking of layers: e.g. cheapRank(1) contains layer number with cheapest N
    integer :: icost
    integer, parameter :: nstep = 100           ! chose this to be flexible again

    real :: dtransp                             ! daily transpiration for this PFT (mm)
    real :: avail_nh4                           ! available NH4 in soil layer (gN/m2)
    real :: avail_no3                           ! available NO3 in soil layer (gN/m2)
    real :: nh4_conc                            ! NH4 concentration (gN/gH2O)
    real :: no3_conc         
    real :: nh4_uptake                          ! (gN) 
    real :: no3_uptake                          ! (gN)
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
    real :: cost(3)                              ! cost of N for each layer and source (-ID) [gC/gN]
    real :: cost_of_n                            ! cheapest cost of N for each layer [gC/gN]
    real :: Cacq                                 ! C spent for N acquisition
    real :: Nacq                                 ! N acquired by spending Cacq
          

    ! xxx try
    soiltemp = 18.0    

    !-------------------------------------------------------------------------
    ! FUN treats different soil layers explicitly. In LPX N-pools represent a
    ! single soil N pool (100 cm). The code below is adopted for a
    ! formulation with different soil layers and soil N pools. To be consistent
    ! with (current) LPX, nsl_eff (effective number of soil layers) is
    ! set to 1. Thus, the model does not distinguish between N acquisition costs
    ! in different layers. An implementation of layer-specific N
    ! uptake, N-inorg pools, root turnover, N mineralization would be nice to have.
    !-------------------------------------------------------------------------


    
    ! Absolute water content corresponding to permanent wilting point [gH2O].
    !epsilon_wtot = Fpwp( nlayers-nsl+1, jpngr ) * sum( Dz_soil(:) ) * 1000.
    !-------------------------------------------------------------------------

    ! LU-PFT association
    lu = lu_category(pft)

    ! xxx think about this: order in which PFTs get access to Ninorg matters!
    
    ! xxx try:
    dtransp = aet_d(lu_category(pft))*fpc_grid(pft,jpngr)

    ! ! Change the order in which pfts get access to N stores daily.
    ! do ppft=1,npft
    !   pft = modulo(ppft+day-2,npft)+1

    if (dnpp(pft)%c12 < 0.0) then
      !-------------------------------------------------------------------------
      ! In case of negative daily NPP, no N uptake is necessary, no C is exuded
      ! npp has been calculated assuming no N limitation. On exit, npp has been reduced to account
      ! for expenditure on N uptake.
      !-------------------------------------------------------------------------
      dcex(pft) = 0.0
     
    else

      !//////////////////////////////////////////////////////////////////////////
      ! INITIALIZATION
      !-------------------------------------------------------------------------
      unmetDemand    = .true.

      Nacq_active    = 0.0                          ! active uptake, sum over sub-timesteps
      Nacq_fix       = 0.0                          ! N fixation, sum over sub-timesteps
      Cacq_active    = 0.0
      Cacq_fix       = 0.0
       
      no3_conc       = 9999.e27                     ! assign dummy value before variable is passed to SR
      nh4_conc       = 9999.e27                     ! assign dummy value before variable is passed to SR
      avail_no3      = 9999.e27                     ! assign dummy value before variable is passed to SR
      avail_nh4      = 9999.e27                     ! assign dummy value before variable is passed to SR
       
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
      n_demand = max( dnpp(pft)%c12/cton_pro(pft), 0.0 )  !in gN/m2/(timestep), timestetp=day or month 

      ! xxx debug
      !print*,'---------------------------------------'
      !print*,'NPP, N demand, Ninorg', dnpp(pft)%c12, n_demand, no3(lu,jpngr)%n14+nh4(lu,jpngr)%n14

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
      ! Just reduce the demand by amount retranslocated, not the labile N pool itself
      !--------------------------------------------------------------------------
      ! xxx debug
      !n_uptake_retrans = min( n_demand, plabl(pft,jpngr)%n%n14 )
      !n_demand_remaining = n_demand_remaining - n_uptake_retrans

      !print*,'N retranslocated, remaining N demand', n_uptake_retrans, n_demand_remaining


      !//////////////////////////////////////////////////////////////////////////
      ! NO3 AND NH4 CONCENTRATION IN SOIL WATER
      !-------------------------------------------------------------------------
      ! N is assumed immobile and unavailable in dry soils. Using the total soil
      ! water and soil N to calculate concentration. Should be equivalent to
      ! assuming concentration of N as being equal in liquid and ice and assuming
      ! only N in liquid is available. nh4_conc and no3_conc are in gN/gH2O.
      !-------------------------------------------------------------------------
      call update_Ninorg( lu, nh4(lu,jpngr)%n14, avail_nh4, nh4_conc )
      call update_Ninorg( lu, no3(lu,jpngr)%n14, avail_no3, no3_conc )

       
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
      ! No active control on passive uptake - always occurrs even if the unmet N
      ! demand is zero.
      !--------------------------------------------------------------------------         
      nh4_uptake = nh4_conc * dtransp
      no3_uptake = no3_conc * dtransp

      ! *** Associate cost of passive NO3 uptake due to pH equilibration by plant (Gutschik, 1981)? ***
       
      nh4_uptake = min( nh4_uptake, avail_nh4 )
      no3_uptake = min( no3_uptake, avail_no3 )
       
      ! xxx add more flexibility by allowing N uptake in excess of demand 
      !n_uptake_pass = nh4_uptake + no3_uptake

      ! Do not take up any N in excess of demand
      n_uptake_pass = min( (nh4_uptake + no3_uptake), n_demand_remaining )

      if ( n_demand_remaining <= n_uptake_pass ) then
        !--------------------------------------------------------------------------
        ! Passive uptake exceeds N-demand. Take up anyway.
        !--------------------------------------------------------------------------
        unmetDemand = .false.
        n_demand_remaining = 0.0

      else
        n_demand_remaining = n_demand_remaining - n_uptake_pass
      end if

       
      !--------------------------------------------------------------------------
      ! Update N taken up and available N in soil.
      !--------------------------------------------------------------------------
      nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - nh4_uptake
      no3(lu,jpngr)%n14 = no3(lu,jpngr)%n14 - no3_uptake
       
      call update_Ninorg( lu, nh4(lu,jpngr)%n14, avail_nh4, nh4_conc )
      call update_Ninorg( lu, no3(lu,jpngr)%n14, avail_no3, no3_conc )


      ! XXX May pools be emptied here?

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
          cost(icostFix) = fun_cost_fix( soiltemp )
        else
          cost(icostFix) = BIG_COST
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
        !nstep = 100
        
        !--------------------------------------------------------------------------
        ! Calculate NPP still to pay for after passive uptake and retranslocation
        ! Minimum function is used here to protect from numerical imprecision making
        ! 'npp_remaining' > 'dnpp' after division and multiplication with cton_pro.
        !--------------------------------------------------------------------------
        npp_remaining = min(n_demand_remaining * cton_pro(pft), dnpp(pft)%c12)

        !--------------------------------------------------------------------------
        ! Split NPP to pay for into steps
        !--------------------------------------------------------------------------
        npp_remaining_step = npp_remaining / real(nstep)

        !//////////////////////////////////////////////////////////////////////////
        ! SUB-TIMESTEP LOOP STARTS HERE
        !--------------------------------------------------------------------------
        substeps: do istep=1,nstep

          !--------------------------------------------------------------------------
          ! Initializations
          !--------------------------------------------------------------------------
          Nacq = 0.0
          Cacq = 0.0
          empty = .false.

          !//////////////////////////////////////////////////////////////////////////
          ! COST OPTIMIZATION
          !--------------------------------------------------------------------------
          ! Calculate the costs of the various possible sources and identify the
          ! cheapest. The cost of nitrogen fixation does not vary within a timestep,
          ! and has been calculated outside this loop.
          !--------------------------------------------------------------------------
             
          !--------------------------------------------------------------------------
          ! Cost of active uptake of NH4, NO3 and organic N
          !--------------------------------------------------------------------------
          cost(icostActiveNH4) = fun_cost_active_nh4( avail_nh4, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )
          cost(icostActiveNO3) = fun_cost_active_no3( avail_no3, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )
          ! cost(icostActiveOrg) = fun_cost_active_org(litterN)
           
          !--------------------------------------------------------------------------
          ! Find the lowest cost in each soil layer and each sub-timestep
          !--------------------------------------------------------------------------
          mloc(1:1) =  MINLOC( cost(:), mask=option(:) )
          iminCost  = mloc(1)
          cost_of_n = cost(iminCost)
                                 
          !--------------------------------------------------------------------------
          ! Spend C on extracting N
          !--------------------------------------------------------------------------
          if ( cost_of_n - BIG_COST > SMALLVALUE ) then
            ! If the cheapest N is very expensive, stores are empty
            ! xxx if cost exceeds maximum, don't take up N and put C to labl                 
            empty = .true.
             
          else
            ! Identify the process involved with the cheapest N
            icost = iminCost
             
            !--------------------------------------------------------------------------
            ! Calculate NPP decrement and N increment (for derivation, see bottom of this
            ! file). Cacq: C spent for N acquisition. Nacq: N acquired
            !--------------------------------------------------------------------------
            Cacq = npp_remaining_step/(cton_pro(pft)/cost(icost)+1.0)
            Nacq = Cacq/cost(icost)
             
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
              if ( Nacq > avail_nh4 ) then
                ! xxx better way to solve this?
                Nacq = avail_nh4
                Cacq = Nacq*cost(icost)
              endif

              nh4(lu,jpngr)%n14 = nh4(lu,jpngr)%n14 - Nacq

              call update_Ninorg( lu, nh4(lu,jpngr)%n14, avail_nh4, nh4_conc )

              Nacq_active = Nacq_active + Nacq
              Cacq_active = Cacq_active + Cacq
              nh4_uptake = nh4_uptake + Nacq
               
            elseif (icost==icostActiveNO3) then
              !//////////////////////////////////////////////////////////////////////////
              ! ACTIVE NO3 UPTAKE
              !--------------------------------------------------------------------------
              ! Decrement available store and update concentration.
              ! If store insufficient, reduce uptake and NPP increment.
              !--------------------------------------------------------------------------
              if ( Nacq > avail_no3 ) then
                ! xxx better way to solve this?
                Nacq = avail_no3
                Cacq = Nacq*cost(icost)
              endif

              no3(lu,jpngr)%n14 = no3(lu,jpngr)%n14 - Nacq

              call update_Ninorg( lu, no3(lu,jpngr)%n14, avail_no3, no3_conc )

              Nacq_active = Nacq_active + Nacq
              Cacq_active = Cacq_active + Cacq
              no3_uptake = no3_uptake + Nacq
               
              ! elseif ( icost == icostActiveOrg ) then
              !//////////////////////////////////////////////////////////////////////////
              ! ACTIVE UPTAKE OF ORGANIC N
              !--------------------------------------------------------------------------
              ! Does it deplete any explecit store?
              ! XXX To be implemented XXX

            else

              ! xxx instead of fixing, put C to labile pool and don't take up any N

              !//////////////////////////////////////////////////////////////////////////
              ! SYMBIOTIC N FIXATION
              !--------------------------------------------------------------------------
              if ( .not. nfixer(pft) ) then
                write(0,*) 'cost:', cost(:)
                write(0,*) 'Non-N-fixer tries to fix N! '//&
                  'BIG_COST too small?'
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

                cost(icostActiveNH4) = fun_cost_active_nh4( avail_nh4, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )

              elseif (icost==icostActiveNO3) then

                cost(icostActiveNO3) = fun_cost_active_no3( avail_no3, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )
                 
                ! If N-org uptake depletes stores of N-org, update costs here
                ! elseif (icost==icostActiveOrg) then

                ! cost(icostActiveOrg) = fun_cost_active_org(litterN)

              endif

              !--------------------------------------------------------------------------
              ! Decide which source of N is now cheapest in this layer. If cheapest layer
              ! is now not the same as the cheapest layer before, re-rank layers.
              !--------------------------------------------------------------------------
              mloc(1:1) = MINLOC( cost(:) )
              iminCost = mloc(1)
              cost_of_n = cost(iminCost)
                              
            endif                                    !  istep<nstep
             
          endif                                      !  cost>big_cost (empty)
           
          !//////////////////////////////////////////////////////////////////////////
          ! UPDATE FLUXES
          !--------------------------------------------------------------------------
          ! If pools got empty during integration loop, that is, before n_demand has
          ! fully been covered, 'empty' is true. In this case, all remaining npp is
          ! exuded. And the integration loop is exited.
          ! xxx think of a better way to do this.
          !--------------------------------------------------------------------------
          if ( empty ) then
             
            ! Do not reduce "real" npp here, but substract exu in npp-SR
            ! npp(pft,1) = npp(pft,1) - npp_remaining
            dcex(pft) = dcex(pft) + npp_remaining

            ! print*,'everything exuded'

            ! Leave integration loop (steps)
            exit substeps
             
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
           
        end do substeps


        !print*,'N by active uptake and fixation, remaining N demand', Nacq_active+Nacq_fix, n_demand_remaining

         
        !  print*,'fraction of npp paid for N uptake',
        !  $             (Cacq_fix+Cacq_active)/tmp
           
      endif                                          !unmetDemand
              

    endif                                              ! (npp(pft).lt.0.0)


    !--------------------------------------------------------------------------
    ! Update annual N-uptake of this PFT. N-retranslocation is not considered
    ! N-uptake.
    !--------------------------------------------------------------------------
    dnup(pft)%n14 = n_uptake_pass + Nacq_active + Nacq_fix  ! n_uptake_retrans is not considered uptake
    dnup_pas(pft) = n_uptake_pass
    dnup_act(pft) = Nacq_active                   
    dnup_fix(pft) = Nacq_fix  
    dnup_ret(pft) = n_uptake_retrans
    if (dnup(pft)%n14>0.0) then
      dccost(pft) = dcex(pft) / dnup(pft)%n14          
    else
      dccost(pft) = -9999
    endif


    !XXX debug
    !print*,'NPP, EXU          ', dnpp(pft)%c12, dcex(pft)
    !print*,'N uptake          ', dnup(pft)%n14
    !print*,'N retranslocated  ', n_uptake_retrans
    !if ((dnup(pft)%n14+n_uptake_retrans)>0.0) then
    !  print*,'C:N after N uptake', (dnpp(pft)%c12 - dcex(pft))/(dnup(pft)%n14+n_uptake_retrans)
    !else
    !  print*,'dnup + n_uptake_retrans',(dnup(pft)%n14+n_uptake_retrans)
    !end if
    !                 print*,'N uptake by BNF',Nacq_fix
    !                 print*,'N uptake by act',Nacq_active
    !print*,'actual CN ratio',(dnpp(pft)%c12 - dcex(pft))/(n_uptake_pass+Nacq_active+Nacq_fix+n_uptake_retrans)
    !                 if (Nacq_fix.gt.0.0) print*,'daily cost of BNF', Cacq_fix/Nacq_fix
    !                 if (Nacq_active.gt.0.0) print*,'daily cost of act', Cacq_active/Nacq_active


    return

  contains

    subroutine update_Ninorg( lu, Ninorg, avail_Ninorg, Ninorg_conc )
      !******************************************************************************
      ! Update N taken up and available N in soil.
      !--------------------------------------------------------------------------
      use fluxes, only: dwn
      implicit none

      integer, intent(in) :: lu
      real, intent(in)  :: Ninorg
      real, intent(out) :: avail_Ninorg
      real, intent(out) :: Ninorg_conc
                    
      Ninorg_conc = Ninorg/dwn(lu)
      
      if ( dwn(lu) > EPSILON_WTOT ) then 
        avail_Ninorg = Ninorg - EPSILON_WTOT*Ninorg_conc
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


  !$$$      subroutine n_fixation_cryptogam( day, lu, jpngr, dnfix_cpc, dnfix_cgc )
  !$$$      !******************************************************************************
  !$$$      ! SUBROUTINE N_UPTAKE BY CRYPTOGAMIC COVERS
  !$$$      !-------------------------------------------------------------------------
  !$$$      ! Simulated to match pattern and global total fixed N after Elbert et al.
  !$$$      ! (2012), Nature Geoscience. Basic assumption: N uptake is driven by energy
  !$$$      ! available (solar radiation ~ photosynthetically active radiation) and not
  !$$$      ! absorbed by leafs or stems. N fixation by cryptogamic ground cover (CGC)
  !$$$      ! thus scales with (1-VPC), where VPC is analogous to FPC but takes into
  !$$$      ! account the shading by branches and stems. N fixation by cryptogamic
  !$$$      ! plant covers (CPC) scales with SPC. 
  !$$$      !-------------------------------------------------------------------------
  !$$$
  !$$$      implicit none
  !$$$
  !$$$      ! ARGUMENTS
  !$$$      INTEGER day, lu, jpngr
  !$$$      REAL*8 dnfix_cpc, dnfix_cgc
  !$$$      
  !$$$      ! LOCAL VARIABLES
  !$$$      INTEGER
  !$$$     $     pft,ppft
  !$$$      
  !$$$      REAL*8
  !$$$     $     fpc_ind,               ! phenology-modulated (!) fractional plant cover
  !$$$     $     local_fpc_grid,        ! FPC w.r.t. grid cell area (is not the same as the global variable fpc_grid)
  !$$$     $     vpc_ind,               ! fractional vegetation cover including stems and branches
  !$$$     $     vpc_grid,              ! VPC w.r.t. grid cell area
  !$$$     $     spc_grid,              ! fractional stem/branches cover
  !$$$     $     fpc_grid_total,        ! fpc_grid summed over all PFTs in the present LU
  !$$$     $     vpc_grid_total,        ! vpc_grid summed over all PFTs in the present LU
  !$$$     $     spc_grid_total,        ! spc_grid summed over all PFTs in the present LU
  !$$$     $     lm_tot(npft),
  !$$$     $     scale
  !$$$
  !$$$      ! Initialisations
  !$$$      vpc_grid_total = 0.
  !$$$      fpc_grid_total = 0.
  !$$$      spc_grid_total = 0.
  !$$$
  !$$$      !  ! Calculate ftemp
  !$$$      !  if (soiltemp.ge.-40.) then
  !$$$      !    tshift = 46.02d0
  !$$$      !    ftemp = exp(308.56d0*(1.0/(20.+tshift)-1.0/
  !$$$      ! $       (soiltemp+tshift)))                             ! Eq.8, XP08 (canexch.cpp:1018)
  !$$$      !  else
  !$$$      !    ftemp = 0.
  !$$$      !  endif
  !$$$      !  ftemp = min(ftemp, 1.)                              ! (canexch.cpp:1023)
  !$$$      !  ftemp = max(ftemp, 0.)                              ! (canexch.cpp:1024)      
  !$$$
  !$$$      do pft=1,npft
  !$$$        if ( present(pft,jpngr) .and. lu_category(pft) .eq. lu ) then
  !$$$
  !$$$        ! LM_TOT
  !$$$        !--------------------------------------------------------------------------
  !$$$        ! Non-linearity of Beer-Law causes very high FPC values when 2 Grasses are present.
  !$$$        ! (Beer Law does NOT make sense for grasses, anyway.)
  !$$$        ! Thus, use sum of all grass/moss-leaf masses and calculate FPC based on the sum.
  !$$$        ! Then compute each PFT's FPC as the product of total-grass FPC times each PFT's leaf mass.
  !$$$        !-------------------------------------------------------------------------
  !$$$          lm_tot(pft) = 0.
  !$$$          ! Grass: C3, C4 on natural, croplands, pasture, peatlands
  !$$$          if (grass(pft)) then
  !$$$            do ppft=1,npft
  !$$$              if (lu_category(ppft).eq.lu_category(pft)) then
  !$$$                if (grass(ppft)) lm_tot(pft) =
  !$$$     $               lm_tot(pft)+lm_ind(ppft,jpngr,1)
  !$$$              endif
  !$$$            enddo
  !$$$          ! Moss: moss on peatlands
  !$$$          elseif (moss(pft)) then
  !$$$            do ppft=1,npft
  !$$$              if (lu_category(ppft).eq.lu_category(pft)) then
  !$$$                if (moss(ppft)) lm_tot(pft) =
  !$$$     $               lm_tot(pft)+lm_ind(ppft,jpngr,1)
  !$$$              endif
  !$$$            enddo
  !$$$          ! Tree: tree on natural lands, peatlands
  !$$$          else
  !$$$            lm_tot(pft) = lm_ind(pft,jpngr,1)
  !$$$          endif
  !$$$          
  !$$$          ! LAI
  !$$$          !--------------------------------------------------------------------------
  !$$$          if (crownarea(pft,jpngr).gt.0.) then
  !$$$            lai_ind(pft,jpngr)=(lm_tot(pft)*sla(pft))/
  !$$$     $           crownarea(pft,jpngr)
  !$$$          else
  !$$$            lai_ind(pft,jpngr)=0.
  !$$$          endif
  !$$$          
  !$$$          ! FPC and VPC
  !$$$          !--------------------------------------------------------------------------
  !$$$          ! Note that this is not identical to how it's calculated in SR update_fpc,
  !$$$          ! where the phenology scaling factor is not included in the exponent.
  !$$$          ! Fractional plant cover accounts for the fraction of the grid cell covered
  !$$$          ! by the photosyntesic plant tissue. To be modulated by daily phenology!
  !$$$          !--------------------------------------------------------------------------
  !$$$          fpc_ind = 1.-dexp(
  !$$$     $                        -1.*kbeer*lai_ind(pft,jpngr)*dphen(day,pft)
  !$$$     $                        )
  !$$$          vpc_ind = 1.-dexp(
  !$$$     $                        -1.*kbeer*(
  !$$$     $                                      lai_ind(pft,jpngr)*dphen(day,pft)
  !$$$     $                                      + pftpar(pft,46)
  !$$$     $                                      )
  !$$$     $                        )
  !$$$          
  !$$$          local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !$$$          vpc_grid       = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !$$$          
  !$$$          if (lm_tot(pft).gt.0.) then
  !$$$            local_fpc_grid = local_fpc_grid*lm_ind(pft,jpngr,1)
  !$$$     $           /lm_tot(pft)
  !$$$            vpc_grid = vpc_grid*lm_ind(pft,jpngr,1)/lm_tot(pft)
  !$$$          else
  !$$$            local_fpc_grid = 0.
  !$$$            vpc_grid       = 0. 
  !$$$          endif
  !$$$
  !$$$          spc_grid = vpc_grid - local_fpc_grid
  !$$$
  !$$$          ! Sum over pfts
  !$$$          !--------------------------------------------------------------------------
  !$$$          fpc_grid_total = fpc_grid_total + local_fpc_grid
  !$$$          vpc_grid_total = vpc_grid_total + vpc_grid
  !$$$          spc_grid_total = spc_grid_total + spc_grid
  !$$$
  !$$$          ! print*,'spc_grid',spc_grid
  !$$$          
  !$$$          !!          call update_fpc(pft,jpngr)
  !$$$          !          
  !$$$          !      ! VAI is analogous to LAI but accounts for stems and branches in addition to
  !$$$          !      ! leafs.
  !$$$          !          vpc_ind = 1. - dexp(
  !$$$          !     $                          - 1.*kbeer*(
  !$$$          !     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
  !$$$          !     $                                         + pftpar(pft,46)
  !$$$          !     $                                         )
  !$$$          !     $                          )
  !$$$          !          vpc_grid = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !$$$          !          vpc_grid_total = vpc_grid_total + vpc_grid
  !$$$          !
  !$$$          !      ! Calculate local FCP treating dphen analogously as for the calulation of VAI:
  !$$$          !      ! FPC = 1-exp(-kbeer*LAI*dphen) instead of FPC = dphen*(1-exp(-kbeer*LAI))
  !$$$          !!           fpc_ind = 1. - dexp(
  !$$$          !!     $                           -1.*kbeer*(
  !$$$          !!     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
  !$$$          !!     $                                         )
  !$$$          !!     $                           )
  !$$$          !          fpc_ind = (1. - dexp(
  !$$$          !     $                           -1.*kbeer*(
  !$$$          !     $                                         lai_ind(pft,jpngr)
  !$$$          !     $                                         )
  !$$$          !     $                           ))!*dphen(day,pft)
  !$$$          !          local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !$$$          !          fpc_grid_total = fpc_grid_total + local_fpc_grid
  !$$$          !
  !$$$          !          print*,'pft',pft
  !$$$          !          print*,'local_fpc_grid     ',local_fpc_grid
  !$$$          !          print*,'fpc_grid(pft,jpngr)',fpc_grid(pft,jpngr)
  !$$$          !          
  !$$$          !      ! Calculate fractional stem/branch cover of grid cell as the difference
  !$$$          !          spc_grid = vpc_grid - local_fpc_grid
  !$$$          !          spc_grid_total = spc_grid_total + spc_grid
  !$$$         
  !$$$        endif
  !$$$      enddo
  !$$$
  !$$$      
  !$$$      if (vpc_grid_total.gt.1.) then
  !$$$        !        print*,'-----------------scaling-------------------'
  !$$$        !        print*,'fpc_grid_total',fpc_grid_total
  !$$$        !        print*,'vpc_grid_total',vpc_grid_total
  !$$$        scale = 1. / vpc_grid_total
  !$$$        fpc_grid_total = fpc_grid_total * scale
  !$$$        vpc_grid_total = vpc_grid_total * scale
  !$$$        spc_grid_total = spc_grid_total * scale
  !$$$        !        print*,'fpc_grid_total',fpc_grid_total
  !$$$        !        print*,'vpc_grid_total',vpc_grid_total
  !$$$      endif
  !$$$
  !$$$      if (fpc_grid_total.gt.1.) then
  !$$$        !        print*,'fpc_grid_total',fpc_grid_total
  !$$$        stop 
  !$$$      endif
  !$$$
  !$$$      ! Daily N fixed by cryptogamic ground and plant covers (communicated to calling SR)
  !$$$      !-------------------------------------------------------------------------
  !$$$      ! Fixation scales with daily photosynthetically active radiation and the
  !$$$      ! branch/stem surface for CPC and the bare ground surface for CGC.
  !$$$      
  !$$$      dnfix_cpc = par_day(day) * max( 0., spc_grid_total) / glob_CPC_scal
  !$$$      dnfix_cgc = par_day(day) * max( 0., (1.0 - vpc_grid_total) ) / glob_CGC_scal
  !$$$
  !$$$      end subroutine n_fixation_cryptogam



end module nuptake_fun