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
  real :: KN_ACTIVE_NINORG  ! N-availability constant in cost function for active NH4 uptake
  real :: KC_ACTIVE_NINORG  ! root exploration constant in cost function for active NH4 uptake
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
    use params_siml, only: spinup
    use params_modl, only: lu_category, cton_pro, nfixer
    use fluxes, only: dnpp, dnup, dcex
    use waterbal_stash, only: aet_d  ! xxx make 'aet_d' a required global variable?
    use treegeometry, only: nind, ispresent, fpc_grid
    use pools, only: pninorg, pexud, plabl, proot, dwtot
    use statevars, only: dtemp_soil

    implicit none
    
    ! ARGUMENTS
    integer, intent(in) :: jpngr, pft
    
    ! LOCAL VARIABLES
    integer, parameter :: nsl_eff = 1           ! effective number of soil layers for vertical distribution of N inorg
    integer, parameter :: icostActiveNinorg = 1 ! process number for active uptake of inorganic N
    integer, parameter :: icostFix = 2          ! process number for symbiotic N-fixation
    integer, parameter :: MAXSTEP = 100         ! maximum number of sub-time steps in N uptake optimisation
    integer, parameter :: PARTS = 100           ! number of parts of split up of inorganic N for uptake

    real, parameter :: SMALLVALUE = 1.e-9       ! to prevent numerical instability
    real, parameter :: BIG_COST = 999999.       ! dummy value

    logical unmetDemand
    logical ranked         
    logical empty
    logical option(2)                           ! whether N-fixation is an option

    integer :: iz, istep, lu   ! counter variables
    integer :: iminCost                         ! info which source (-ID) is cheapest in each layer
    integer :: mloc(1)                          ! used to hold result of minloc
    integer :: cheapRank                        ! ranking of layers: e.g. cheapRank(1) contains layer number with cheapest N
    integer :: icost
    integer, parameter :: nstep = 100           ! chose this to be flexible again
    integer, save :: invocation = 0
    logical :: continue_uptake = .true.
    integer, parameter :: calib_year = 99999*ndayyear   !300*ndayyear xxx change calib_year to 300 again xxx
    real, dimension(npft,maxgrid), save :: aNacq_active
    real, dimension(npft,maxgrid), save :: aCacq_active
    real, dimension(npft,maxgrid), save :: amean_cost_pvy
    real :: dmean_cost
    character(len=1) :: sub
    real :: acc

    real :: dtransp                             ! daily transpiration for this PFT (mm)
    real :: avail_ninorg                        ! available inorganic N in soil layer (gN/m2)
    real :: ninorg_conc                         ! inorganic N concentration (gN/gH2O)
    real :: ninorg_uptake                       ! (gN) 
    real :: n_uptake_pass                       ! (gN)
    real :: npp_remaining
    real :: npp_remaining_step
    real :: dNacq_active
    real :: dNacq_fix
    real :: dCacq_active
    real :: dCacq_fix
    real :: n_demand
    real :: n_demand_remaining
    real :: n_uptake_retrans
    real :: cost(2)                              ! cost of N for each source (-ID) [gC/gN]
    real :: cost_of_n                            ! cheapest cost of N for each layer [gC/gN]
    real :: Cacq                                 ! C spent for N acquisition
    real :: Nacq                                 ! N acquired by spending Cacq
          

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

    if (pft==npft) invocation = invocation + 1

    if (spinup) then
      acc = max( 0.0, 2.0 - real(invocation/ndayyear)/100 ) 
    else
      acc = 0.0
    end if

    ! xxx try
    ! print*,'acc',acc
    ! pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + acc
    ! if (spinup.and.invocation/ndayyear<200) then
    !   pninorg(lu,jpngr)%n14 = 0.05
    ! end if
    pninorg(lu,jpngr)%n14 = 0.05


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

      dNacq_active    = 0.0                          ! active uptake, sum over sub-timesteps
      dNacq_fix       = 0.0                          ! N fixation, sum over sub-timesteps
      dCacq_active    = 0.0
      dCacq_fix       = 0.0
       
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
      ! INORGANIC N CONCENTRATION IN SOIL WATER
      !-------------------------------------------------------------------------
      ! N is assumed immobile and unavailable in dry soils. Using the total soil
      ! water and soil N to calculate concentration. Should be equivalent to
      ! assuming concentration of N as being equal in liquid and ice and assuming
      ! only N in liquid is available. nh4_conc and no3_conc are in gN/gH2O.
      !-------------------------------------------------------------------------
      call update_Ninorg( lu, pninorg(lu,jpngr)%n14, avail_ninorg, ninorg_conc, dwtot(lu,jpngr) )

      ! xxx debug
      !print*,'dwtot(lu,jpngr)  ',dwtot(lu,jpngr)

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
      ninorg_uptake = ninorg_conc * dtransp

      ! *** Associate cost of passive NO3 uptake due to pH equilibration by plant (Gutschik, 1981)? ***
       
      ninorg_uptake = min( ninorg_uptake, avail_ninorg )
       
      ! xxx add more flexibility by allowing N uptake in excess of demand 

      ! Do not take up any N in excess of demand
      n_uptake_pass = min( ninorg_uptake, n_demand_remaining )

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
      pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - ninorg_uptake
       
      call update_Ninorg( lu, pninorg(lu,jpngr)%n14, avail_ninorg, ninorg_conc, dwtot(lu,jpngr) )


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


      if (unmetDemand.and.invocation<=calib_year) then
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
          cost(icostFix) = fun_cost_fix( dtemp_soil(lu,jpngr) )
        else
          cost(icostFix) = BIG_COST
          option(icostFix) = .false.
        endif

        !//////////////////////////////////////////////////////////////////////////
        ! SUB-TIMESTEP LOOP STARTS HERE
        !--------------------------------------------------------------------------
        empty = .false.
        istep = 0
        continue_uptake = .true.

        do while (continue_uptake)

          istep = istep + 1

          cost(icostActiveNinorg) = fun_cost_active_ninorg( avail_ninorg, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )

          ! C:N ratio constraint 
          Cacq = npp_remaining_step/(cton_pro(pft)/cost(icostActiveNinorg)+1.0)

          Nacq = Cacq/cost(icostActiveNinorg)
          if ( Nacq > avail_ninorg ) then
            Nacq = avail_ninorg
            Cacq = Nacq*cost(icostActiveNinorg)
            empty = .true.
          endif

          dNacq_active = dNacq_active + Nacq
          dCacq_active = dCacq_active + Cacq
          ninorg_uptake = ninorg_uptake + Nacq
 
          pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - Nacq
          call update_Ninorg( lu, pninorg(lu,jpngr)%n14, avail_ninorg, ninorg_conc, dwtot(lu,jpngr) )

          npp_remaining = npp_remaining - Cacq

          dcex(pft) = dcex(pft) + Cacq

          dmean_cost = dCacq_active / dNacq_active

          if (istep>=nstep.or.empty) continue_uptake = .false.

        end do
           
      else if (unmetDemand.and.invocation>calib_year) then

        !//////////////////////////////////////////////////////////////////////////
        ! SUB-TIMESTEP LOOP STARTS HERE
        !--------------------------------------------------------------------------
        empty = .false.
        istep = 0
        continue_uptake = .true.

        do while (continue_uptake)

          istep = istep + 1

          if (dNacq_active>0.0) then
            dmean_cost = dCacq_active / dNacq_active
          else
            dmean_cost = 0.0
          end if

          if (dmean_cost>amean_cost_pvy(pft,jpngr)) then
            ! print*,'amean_cost_pvy',amean_cost_pvy(pft,jpngr)
            ! print*,'dmean_cost    ',dmean_cost
            continue_uptake = .false.
            Cacq = 0.0
          else
            Cacq = npp_remaining_step
            ! print*,'dmean_cost    ',dmean_cost
            ! print*,'amean_cost_pvy',amean_cost_pvy(pft,jpngr)
          end if


          cost(icostActiveNinorg) = fun_cost_active_ninorg( avail_ninorg, proot(pft,jpngr)%c%c12 * nind(pft,jpngr) )
          Nacq = Cacq/cost(icostActiveNinorg)

          if ( Nacq > avail_ninorg ) then
            Nacq = avail_ninorg
            Cacq = Nacq*cost(icostActiveNinorg)
            empty = .true.
          endif

          dNacq_active = dNacq_active + Nacq
          dCacq_active = dCacq_active + Cacq
          ninorg_uptake = ninorg_uptake + Nacq
 
          pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - Nacq
          call update_Ninorg( lu, pninorg(lu,jpngr)%n14, avail_ninorg, ninorg_conc, dwtot(lu,jpngr) )

          npp_remaining = npp_remaining - Cacq
          dcex(pft) = dcex(pft) + Cacq

          if (empty) continue_uptake = .false.

        end do

        ! print*,'invocation', invocation

      endif                                          !unmetDemand


    endif                                              ! (npp(pft).lt.0.0)


    !--------------------------------------------------------------------------
    ! Update N-uptake of this PFT. N-retranslocation is not considered
    ! N-uptake.
    !--------------------------------------------------------------------------
    ! daily
    dnup(pft)%n14 = n_uptake_pass + dNacq_active + dNacq_fix  ! n_uptake_retrans is not considered uptake
    dnup_pas(pft) = n_uptake_pass
    dnup_act(pft) = dNacq_active                   
    dnup_fix(pft) = dNacq_fix  
    dnup_ret(pft) = n_uptake_retrans
    if (dnup(pft)%n14>0.0) then
      dccost(pft) = dmean_cost       
    else
      dccost(pft) = -9999
    endif

    if (mod((invocation-1),ndayyear)==0) then
      !f first day of the year
      aNacq_active(pft,jpngr) = 0.0
      aCacq_active(pft,jpngr) = 0.0

    else if (mod(invocation,ndayyear)==0) then
      ! last day of the year
      amean_cost_pvy(pft,jpngr) = aCacq_active(pft,jpngr) / aNacq_active(pft,jpngr)
      ! print*,'amean_cost_pvy(pft,jpngr)',amean_cost_pvy(pft,jpngr)

    end if
    aNacq_active(pft,jpngr) = aNacq_active(pft,jpngr) + dNacq_active
    aCacq_active(pft,jpngr) = aCacq_active(pft,jpngr) + dCacq_active



    !XXX debug
    ! print*,'NPP, EXU          ', dnpp(pft)%c12, dcex(pft)
    ! print*,'N uptake          ', dnup(pft)%n14
    ! print*,'N retranslocated  ', n_uptake_retrans
    ! if ((dnup(pft)%n14+n_uptake_retrans)>0.0) then
    !  print*,'C:N after N uptake', (dnpp(pft)%c12 - dcex(pft))/(dnup(pft)%n14+n_uptake_retrans)
    ! end if
    !                 print*,'N uptake by BNF',dNacq_fix
    !                 print*,'N uptake by act',dNacq_active
    !print*,'actual CN ratio',(dnpp(pft)%c12 - dcex(pft))/(n_uptake_pass+dNacq_active+dNacq_fix+n_uptake_retrans)
    !                 if (dNacq_fix.gt.0.0) print*,'daily cost of BNF', dCacq_fix/dNacq_fix
    !                 if (dNacq_active.gt.0.0) print*,'daily cost of act', dCacq_active/dNacq_active


    return

  contains

    subroutine update_Ninorg( lu, Ninorg, avail_Ninorg, Ninorg_conc, wtot )
      !**************************************************************************
      ! Update N taken up and available N in soil.
      !--------------------------------------------------------------------------
      implicit none

      integer, intent(in) :: lu
      real, intent(in)  :: Ninorg
      real, intent(out) :: avail_Ninorg
      real, intent(out) :: Ninorg_conc
      real, intent(in)  :: wtot           ! total soil water content (mm)
                    
      Ninorg_conc = Ninorg/wtot

      if ( wtot > EPSILON_WTOT ) then 
        avail_Ninorg = Ninorg - EPSILON_WTOT*Ninorg_conc
      else
        avail_Ninorg = 0.0
        Ninorg_conc  = 0.0
      endif

      return

    end subroutine update_Ninorg


    function fun_cost_active_ninorg( avail_ninorg, croot )
      !******************************************************************************
      ! Cost of active inorganic N uptake by exudation of labile C. 
      ! From Fisher er al., 2010.
      !--------------------------------------------------------------------------      
      use params_modl
      implicit none

      ! arguments
      real, intent(in) :: avail_ninorg
      real, intent(in) :: croot

      ! function return variable
      real :: fun_cost_active_ninorg         ! gC/gN

      if ( croot > SMALLVALUE ) then
        if ( avail_ninorg > SMALLVALUE ) then
          fun_cost_active_ninorg = KN_ACTIVE_NINORG / avail_ninorg * KC_ACTIVE_NINORG / croot 
          ! fun_cost_active_ninorg = KN_ACTIVE_NH4 / avail_nh4 + KC_ACTIVE_NH4 / croot ! xxx in FUN 2.0 it's a sum
        else
          fun_cost_active_ninorg = BIG_COST
        endif

      else
        fun_cost_active_ninorg = BIG_COST
      endif

    end function fun_cost_active_ninorg


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
    integer, parameter :: npar = 6
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

    ! N-availability constant in cost function for active inorganic N uptake
    KN_ACTIVE_NINORG = params_array(4)
    
    ! root exploration constant in cost function for active inorganic N uptake
    KC_ACTIVE_NINORG = params_array(5)
    
    ! minimum soil water content above which N becomes available 
    EPSILON_WTOT = params_array(6)

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