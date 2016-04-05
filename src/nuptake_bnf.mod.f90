module md_nuptake
  !////////////////////////////////////////////////////////////////
  ! FUN NITROGEN UPTAKE MODULE
  ! Contains the "main" subroutine 'nuptake' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'nuptake' must contain this list 
  ! of subroutines (names that way).
  !   - nuptake
  !   - getpar_modl_nuptake
  !   - initio_nuptake
  !   - initoutput_nuptake
  !   - getout_daily_nuptake
  !   - getout_monthly_nuptake
  !   - writeout_ascii_nuptake
  ! Required module-independent model state variables (necessarily 
  ! updated by 'nuptake') are:
  !   - daily NPP ('dnpp')
  !   - soil temperature ('xxx')
  !   - inorganic N _pools ('no3', 'nh4')
  !   - xxx 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: ndayyear, nmonth, nlu, npft, maxgrid
  use md_classdefs

  implicit none

  private
  public nuptake, getpar_modl_nuptake, initdaily_nuptake, initio_nuptake, &
    initoutput_nuptake, getout_daily_nuptake, writeout_ascii_nuptake, &
    calc_dnup, outtype_calc_dnup

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type params_nuptake_type
    real :: eff_nup           ! uptake efficiency for equation
    real :: minimumcostfix    ! minimum cost of N fixation (at optimal temperature)
    real :: fixoptimum        ! optimum temperature for N fixation
    real :: fixwidth          ! shape parameter for width of N fixation cost function
  end type params_nuptake_type

  type( params_nuptake_type ) :: params_nuptake

  !----------------------------------------------------------------
  ! module-specific (private) variables
  !----------------------------------------------------------------
  real, dimension(npft) :: dccost           ! daily mean C cost of N uptake [gC/gN] 
  real, dimension(npft) :: dnup_pas         ! daily passive N uptake [gN/m2/d]
  real, dimension(npft) :: dnup_act         ! daily active N uptake [gN/m2/d]  
  real, dimension(npft) :: dnup_fix         ! daily N uptake by plant symbiotic N fixation [gN/m2/d]
  real, dimension(npft) :: dnup_ret         ! daily N uptake [gN/m2/d]

  type outtype_calc_dnup
    real :: act
    real :: fix
  end type outtype_calc_dnup

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdccost   ! daily mean C cost of N uptake (gC/gN) 
  real, allocatable, dimension(:,:,:) :: outdnup_pas
  real, allocatable, dimension(:,:,:) :: outdnup_act
  real, allocatable, dimension(:,:,:) :: outdnup_fix
  real, allocatable, dimension(:,:,:) :: outdnup_ret


contains


  subroutine nuptake( jpngr )
    !/////////////////////////////////////////////////////////////////
    ! SUBROUTINE NUPTAKE ASSUMING CONSTANT EXUDATION PER UNIT ROOT MASS
    !-----------------------------------------------------------------
    ! This model calculates first the passive uptake of N via the 
    ! transpiration stream.
    !-----------------------------------------------------------------
    use md_classdefs
    use md_plant, only: dcex, dnup, params_pft_plant, ispresent, plabl
    use md_ntransform, only: pninorg

    ! arguments
    integer, intent(in) :: jpngr

    ! local variables
    integer :: lu, pft
    real    :: avail_ninorg                        ! available inorganic N in soil layer (gN/m2)
    real    :: ninorg_conc                         ! inorganic N concentration (gN/gH2O)
    real    :: n_uptake_pass                       ! (gN)
    real    :: n_uptake_retrans
    real    :: dNacq_act
    real    :: dNacq_fix
    real    :: dmean_cost

    type( outtype_calc_dnup ) :: out_calc_dnup

    ! xxx debug
    real :: test_eff_bnf
    real :: test_eff_act

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft

      if ( ispresent(pft,jpngr) ) then

        lu = params_pft_plant(pft)%lu_category

        ! ! xxx try:
        ! dtransp = daet(lu_category(pft)) * fpc_grid(pft,jpngr)

        !//////////////////////////////////////////////////////////////////////////
        ! INITIALIZATION
        !-------------------------------------------------------------------------
        n_uptake_pass = 0.0
        
        dNacq_act = 0.0                          ! active uptake, sum over sub-timesteps
        dNacq_fix = 0.0                          ! N fixation, sum over sub-timesteps
        dCexu     = 0.0

        if ( dcex(pft)>0.0 ) then
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
          ! n_uptake_retrans = min( n_demand, plabl(pft,jpngr)%n%n14 )
          ! n_demand_remaining = n_demand_remaining - n_uptake_retrans


          ! !//////////////////////////////////////////////////////////////////////////
          ! ! PASSIVE UPTAKE
          ! ! No active control on passive uptake - always occurrs even if the unmet N
          ! ! demand is zero.
          ! !--------------------------------------------------------------------------
          ! n_uptake_pass = ninorg_conc * dtransp(pft) / 1000.0     ! [dtransp] = g H2O; [ninorg_conc] = g N / (mm H2O) = g N / (kg H2O)
          ! n_uptake_pass = min( n_uptake_pass, avail_ninorg )

          ! write(0,*) 'n_uptake_pass ',n_uptake_pass 

          ! Update
          pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - n_uptake_pass
          ! avail_ninorg          = calc_avail_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   
          ! ninorg_conc           = calc_conc_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   

          ! write(0,*) 'avail_ninorg',avail_ninorg
          ! write(0,*) 'ninorg_conc ',ninorg_conc 

          !//////////////////////////////////////////////////////////////////////////
          ! ACTIVE UPTAKE
          ! Active N uptake is a function of initial N available and C exuded
          !--------------------------------------------------------------------------
          out_calc_dnup = calc_dnup( dcex(pft), pninorg(lu,jpngr)%n14 )

          ! write(0,*) 'dcex(pft)      ', dcex(pft)      
          ! write(0,*) 'in SR nuptake: dcex(pft)            ',dcex(pft)      
          ! write(0,*) 'in SR nuptake: dNacq_act          ',dNacq_act 
          ! write(0,*) 'in SR nuptake: pninorg(lu,jpngr)%n14 ',pninorg(lu,jpngr)%n14 

          if ((out_calc_dnup%act+out_calc_dnup%fix)>0.0) then
            dmean_cost = dcex(pft)  / (out_calc_dnup%act+out_calc_dnup%fix)
          else
            dmean_cost = 9999.0
          end if

          ! Update
          pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - out_calc_dnup%act

        end if

        !--------------------------------------------------------------------------
        ! Update N-uptake of this PFT. N-retranslocation is not considered
        ! N-uptake.
        !--------------------------------------------------------------------------
        ! daily
        dnup(pft)%n14 = n_uptake_pass + out_calc_dnup%act + out_calc_dnup%fix  ! n_uptake_retrans is not considered uptake
        dnup_pas(pft) = n_uptake_pass
        dnup_act(pft) = out_calc_dnup%act                   
        dnup_fix(pft) = out_calc_dnup%fix  
        dnup_ret(pft) = n_uptake_retrans
        if (dnup(pft)%n14>0.0) then
          dccost(pft) = dmean_cost       
        else
          dccost(pft) = 0.0
        endif

        !--------------------------------------------------------------------------
        ! N acquisition to labile pool
        !--------------------------------------------------------------------------
        call ncp( dnup(pft), plabl(pft,jpngr)%n )

      end if 

    end do
    ! write(0,*) '---- finished nuptake'

  end subroutine nuptake


  function calc_dnup( cexu, n0, soiltemp ) result( out_dnup )
    !/////////////////////////////////////////////////////////////////
    ! With a FUN-like approach:
    ! dCexu/dNup = K / (N0 - Nup); K=1/eff_nup
    ! => Nup(Cexu) = N0 * ( 1.0 - exp( - eff_nup * cexu ) )
    !-----------------------------------------------------------------
    ! arguments
    real, intent(in) :: cexu      ! C exuded (gC/m2/d)
    real, intent(in) :: n0        ! initial available N (gN/m2)
    real, intent(in) :: soiltemp  ! soil temperature (deg C)

    ! function return variable
    type( outtype_calc_dnup ) :: out_calc_dnup

    ! local variables
    real :: mydnup_act
    real :: mydnup_fix
    real :: cost_bnf
    real :: eff_bnf
    real :: cexu_act
    real :: cexu_bnf

    !-----------------------------------------------------------------
    ! get cost of BNF at this soil temperature. Cost = dCex / dNfix
    !-----------------------------------------------------------------
    cost_bnf = fun_cost_fix( soiltemp )

    !-----------------------------------------------------------------
    ! get inverse of cost = efficiency: eff_bnf = dNfix / dCex
    !-----------------------------------------------------------------
    eff_bnf = 1.0 / cost_bnf

    !-----------------------------------------------------------------
    ! Find amount of active uptake (~Cex) for which eff_act = eff_bnf
    ! eff_act = dNup_act / dCex
    ! Nup_act = N0 * ( 1.0 - exp( -K * Cex ) )
    ! dNup_act / dCex = K * exp( -K * Cex)
    ! dNup_act / dCex = eff_bnf
    ! ==> Cex = - 1/K * ln( bnf_eff/K )
    !-----------------------------------------------------------------
    cexu_act = -1.0 / eff_nup * log( eff_bnf / ( n0 * eff_nup ) )

    if (cexu_act < cexu) then
      !-----------------------------------------------------------------
      ! Remaining Cex is consumed by N fixing processes
      !-----------------------------------------------------------------
      cexu_bnf = cexu - cexu_act

      !-----------------------------------------------------------------
      ! N uptake via BNF
      !-----------------------------------------------------------------
      out_calc_dnup%fix = cexu_bnf * eff_bnf

      !-----------------------------------------------------------------
      ! N uptake via active uptake
      !-----------------------------------------------------------------
      out_calc_dnup%act = n0 * ( 1.0 - exp( - eff_nup * cexu_act ) )

    else

      out_calc_dnup%fix = 0.0
      out_calc_dnup%act = n0 * ( 1.0 - exp( - eff_nup * cexu ) )

    end if


  end function calc_dnup


  function fun_cost_fix( soiltemp )
    !////////////////////////////////////////////////////////////////
    ! Cost of symbiotic N fixation is the inverse of nitrogenase activity
    ! after Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    !---------------------------------------------------------------- 
    use md_params_modl

    real, intent(in) :: soiltemp
    real :: fun_cost_fix                 ! function return variable

    fun_cost_fix = MINIMUMCOSTFIX + exp((soiltemp-FIXOPTIMUM)**2/(2*FIXWIDTH**2))    ! inverse gauss function  (take WARMEST layer)

  end function fucost_fix
  

  ! function calc_avail_ninorg( ninorg, wtot ) result( avail_ninorg )
  !   !//////////////////////////////////////////////////////////////////////////
  !   ! Returns N available for uptake accounting for constraint of mobility by
  !   ! soil moisture.
  !   !--------------------------------------------------------------------------
  !   ! arguments
  !   real, intent(in)  :: ninorg 
  !   real, intent(in)  :: wtot           ! total soil water content (mm)

  !   ! function return value
  !   real, intent(out) :: avail_ninorg

  !   if ( wtot > EPSILON_WTOT ) then 
  !     avail_ninorg = ninorg - EPSILON_WTOT * ninorg / wtot
  !   else
  !     avail_ninorg = 0.0
  !   endif

  ! end function calc_avail_ninorg


  ! function calc_conc_ninorg( ninorg, wtot ) result( conc_ninorg )
  !   !//////////////////////////////////////////////////////////////////////////
  !   ! Returns N available for uptake accounting for constraint of mobility by
  !   ! soil moisture.
  !   !--------------------------------------------------------------------------
  !   ! arguments
  !   real, intent(in)  :: ninorg 
  !   real, intent(in)  :: wtot           ! total soil water content (mm)

  !   ! function return value
  !   real, intent(out) :: conc_ninorg

  !   if ( wtot > EPSILON_WTOT ) then 
  !     conc_ninorg = ninorg / wtot
  !   else
  !     conc_ninorg = 0.0
  !   endif

  ! end function calc_conc_ninorg


  subroutine initdaily_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Initialise daily variables with zero
    !----------------------------------------------------------------
    dnup_pas(:)    = 0.0
    dnup_act(:)    = 0.0
    dnup_fix(:)    = 0.0
    dnup_ret(:)    = 0.0

  end subroutine initdaily_nuptake


  subroutine initio_nuptake()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

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

    ! C EXUDATION
    filnam=trim(prefix)//'.d.cex.out'
    open(105,file=filnam,err=888,status='unknown')


    !----------------------------------------------------------------
    ! ANNUAL OUTPUT
    !----------------------------------------------------------------

    ! ANNUAL TOTAL C EXUDATION
    filnam=trim(prefix)//'.a.cex.out'
    open(405,file=filnam,err=888,status='unknown')


    return

    888  stop 'INITIO_NUPTAKE: error opening output files'

  end subroutine initio_nuptake



  subroutine getpar_modl_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! uptake efficiency for equation
    ! dCexu/dNup = K / (N0 - Nup); K=1/eff_nup
    params_nuptake%eff_nup = getparreal( 'params/params_nuptake_constexu.dat', 'eff_nup' )

    ! shape parameter of cost function of N fixation 
    ! Below parameters (minimumcostfix, fixoptimum, fixwidth ) are based on 
    ! the assumption that the cost of symbiotic N fixation is the 
    ! inverse of nitrogenase activity. 
    ! After Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    params_nuptake%minimumcostfix = getparreal( 'params/params_nuptake_constexu.dat', 'minimumcostfix' )

    ! shape parameter of cost function of N fixation 
    params_nuptake%fixoptimum = getparreal( 'params/params_nuptake_constexu.dat', 'fixoptimum' )
 
    ! shape parameter of cost function of N fixation 
    params_nuptake%fixwidth = getparreal( 'params/params_nuptake_constexu.dat', 'fixwidth' )


  end subroutine getpar_modl_nuptake


  subroutine initdaily_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Initialise daily variables with zero
    !----------------------------------------------------------------
    dnup_pas(:) = 0.0
    dnup_act(:) = 0.0
    dnup_fix(:) = 0.0
    dnup_ret(:) = 0.0

  end subroutine initdaily_nuptake


  subroutine initio_nuptake()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_params_siml, only: runname, loutnuptake

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

    if (loutnuptake) then
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

    end if

    return

    888  stop 'INITIO_NUPTAKE: error opening output files'

  end subroutine initio_nuptake


  subroutine initoutput_nuptake
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    use md_params_siml, only: init, loutnuptake

    if (loutnuptake) then

      if (init) allocate( outdccost   (npft,ndayyear,maxgrid) ) ! daily mean C cost of N uptake (gC/gN) 
      if (init) allocate( outdnup_pas (npft,ndayyear,maxgrid) )
      if (init) allocate( outdnup_act (npft,ndayyear,maxgrid) )
      if (init) allocate( outdnup_fix (npft,ndayyear,maxgrid) )
      if (init) allocate( outdnup_ret (npft,ndayyear,maxgrid) )

      outdccost  (:,:,:) = 0.0 ! daily mean C cost of N uptake (gC/gN) 
      outdnup_pas(:,:,:) = 0.0
      outdnup_act(:,:,:) = 0.0
      outdnup_fix(:,:,:) = 0.0
      outdnup_ret(:,:,:) = 0.0

    end if

  end subroutine initoutput_nuptake



  subroutine getout_daily_nuptake( jmoy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use md_vars_core, only: dcex 

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    outdccost(:,doy,jpngr)   = dccost(:)
    outdnup_pas(:,doy,jpngr) = dnup_pas(:)
    outdnup_act(:,doy,jpngr) = dnup_act(:)
    outdnup_fix(:,doy,jpngr) = dnup_fix(:)
    outdnup_ret(:,doy,jpngr) = dnup_ret(:)
    outdcex(:,doy,jpngr)     = dcex(:)

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    outacex(:,jpngr) = outacex(:,jpngr) + dcex(:)


  end subroutine getout_daily_nuptake


  subroutine writeout_ascii_nuptake( year )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use md_params_core, only: ndayyear, npft, nlu
    use md_params_siml, only: firstyeartrend, spinupyears, daily_out_startyr, &
      daily_out_endyr, outyear

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real :: itime
    integer :: day, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii_nuptake: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (loutnuptake) then
      if ( .not. spinup .and. outyear>=daily_out_startyr .and. outyear<=daily_out_endyr ) then
        ! Write daily output only during transient simulation
        do day=1,ndayyear

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

          if (nlu>1) stop 'writeout_ascii_nuptake: write out lu-area weighted sum'
          if (npft>1) stop 'writeout_ascii_nuptake: think of something for ccost output'

          ! xxx lu-area weighted sum if npft>0
          write(400,999) itime, sum(outdccost(:,day,jpngr)) / real(npft) 
          write(401,999) itime, sum(outdnup_pas(:,day,jpngr))
          write(402,999) itime, sum(outdnup_act(:,day,jpngr))
          write(403,999) itime, sum(outdnup_fix(:,day,jpngr))
          write(404,999) itime, sum(outdnup_ret(:,day,jpngr))

        end do
      end if
    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_nuptake


  ! subroutine n_fixation_cryptogam( day, lu, jpngr, dnfix_cpc, dnfix_cgc )
  ! !******************************************************************************
  ! ! SUBROUTINE N_UPTAKE BY CRYPTOGAMIC COVERS
  ! !-------------------------------------------------------------------------
  ! ! Simulated to match pattern and global total fixed N after Elbert et al.
  ! ! (2012), Nature Geoscience. Basic assumption: N uptake is driven by energy
  ! ! available (solar radiation ~ photosynthetically active radiation) and not
  ! ! absorbed by leafs or stems. N fixation by cryptogamic ground cover (CGC)
  ! ! thus scales with (1-VPC), where VPC is analogous to FPC but takes into
  ! ! account the shading by branches and stems. N fixation by cryptogamic
  ! ! plant covers (CPC) scales with SPC. 
  ! !-------------------------------------------------------------------------
  !
  ! implicit none
  !
  ! ! ARGUMENTS
  ! INTEGER day, lu, jpngr
  ! REAL*8 dnfix_cpc, dnfix_cgc
  ! 
  ! ! LOCAL VARIABLES
  ! INTEGER
  !$     pft,ppft
  ! 
  ! REAL*8
  !$     fpc_ind,               ! phenology-modulated (!) fractional plant cover
  !$     local_fpc_grid,        ! FPC w.r.t. grid cell area (is not the same as the global variable fpc_grid)
  !$     vpc_ind,               ! fractional vegetation cover including stems and branches
  !$     vpc_grid,              ! VPC w.r.t. grid cell area
  !$     spc_grid,              ! fractional stem/branches cover
  !$     fpc_grid_total,        ! fpc_grid summed over all PFTs in the present LU
  !$     vpc_grid_total,        ! vpc_grid summed over all PFTs in the present LU
  !$     spc_grid_total,        ! spc_grid summed over all PFTs in the present LU
  !$     lm_tot(npft),
  !$     scale
  !
  ! ! Initialisations
  ! vpc_grid_total = 0.
  ! fpc_grid_total = 0.
  ! spc_grid_total = 0.
  !
  ! !  ! Calculate ftemp
  ! !  if (soiltemp.ge.-40.) then
  ! !    tshift = 46.02d0
  ! !    ftemp = exp(308.56d0*(1.0/(20.+tshift)-1.0/
  ! ! $       (soiltemp+tshift)))                             ! Eq.8, XP08 (canexch.cpp:1018)
  ! !  else
  ! !    ftemp = 0.
  ! !  endif
  ! !  ftemp = min(ftemp, 1.)                              ! (canexch.cpp:1023)
  ! !  ftemp = max(ftemp, 0.)                              ! (canexch.cpp:1024)      
  !
  ! do pft=1,npft
  !   if ( present(pft,jpngr) .and. lu_category(pft) .eq. lu ) then
  !
  !   ! LM_TOT
  !   !--------------------------------------------------------------------------
  !   ! Non-linearity of Beer-Law causes very high FPC values when 2 Grasses are present.
  !   ! (Beer Law does NOT make sense for grasses, anyway.)
  !   ! Thus, use sum of all grass/moss-leaf masses and calculate FPC based on the sum.
  !   ! Then compute each PFT's FPC as the product of total-grass FPC times each PFT's leaf mass.
  !   !-------------------------------------------------------------------------
  !     lm_tot(pft) = 0.
  !     ! Grass: C3, C4 on natural, croplands, pasture, peatlands
  !     if (grass(pft)) then
  !       do ppft=1,npft
  !         if (lu_category(ppft).eq.lu_category(pft)) then
  !           if (grass(ppft)) lm_tot(pft) =
  !$               lm_tot(pft)+lm_ind(ppft,jpngr,1)
  !         endif
  !       enddo
  !     ! Moss: moss on peatlands
  !     elseif (moss(pft)) then
  !       do ppft=1,npft
  !         if (lu_category(ppft).eq.lu_category(pft)) then
  !           if (moss(ppft)) lm_tot(pft) =
  !$               lm_tot(pft)+lm_ind(ppft,jpngr,1)
  !         endif
  !       enddo
  !     ! Tree: tree on natural lands, peatlands
  !     else
  !       lm_tot(pft) = lm_ind(pft,jpngr,1)
  !     endif
  !     
  !     ! LAI
  !     !--------------------------------------------------------------------------
  !     if (crownarea(pft,jpngr).gt.0.) then
  !       lai_ind(pft,jpngr)=(lm_tot(pft)*sla(pft))/
  !$           crownarea(pft,jpngr)
  !     else
  !       lai_ind(pft,jpngr)=0.
  !     endif
  !     
  !     ! FPC and VPC
  !     !--------------------------------------------------------------------------
  !     ! Note that this is not identical to how it's calculated in SR update_fpc,
  !     ! where the phenology scaling factor is not included in the exponent.
  !     ! Fractional plant cover accounts for the fraction of the grid cell covered
  !     ! by the photosyntesic plant tissue. To be modulated by daily phenology!
  !     !--------------------------------------------------------------------------
  !     fpc_ind = 1.-dexp(
  !$                        -1.*kbeer*lai_ind(pft,jpngr)*dphen(day,pft)
  !$                        )
  !     vpc_ind = 1.-dexp(
  !$                        -1.*kbeer*(
  !$                                      lai_ind(pft,jpngr)*dphen(day,pft)
  !$                                      + pftpar(pft,46)
  !$                                      )
  !$                        )
  !     
  !     local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !     vpc_grid       = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !     
  !     if (lm_tot(pft).gt.0.) then
  !       local_fpc_grid = local_fpc_grid*lm_ind(pft,jpngr,1)
  !$           /lm_tot(pft)
  !       vpc_grid = vpc_grid*lm_ind(pft,jpngr,1)/lm_tot(pft)
  !     else
  !       local_fpc_grid = 0.
  !       vpc_grid       = 0. 
  !     endif
  !
  !     spc_grid = vpc_grid - local_fpc_grid
  !
  !     ! Sum over pfts
  !     !--------------------------------------------------------------------------
  !     fpc_grid_total = fpc_grid_total + local_fpc_grid
  !     vpc_grid_total = vpc_grid_total + vpc_grid
  !     spc_grid_total = spc_grid_total + spc_grid
  !
  !     ! print*,'spc_grid',spc_grid
  !     
  !     !!          call update_fpc(pft,jpngr)
  !     !          
  !     !      ! VAI is analogous to LAI but accounts for stems and branches in addition to
  !     !      ! leafs.
  !     !          vpc_ind = 1. - dexp(
  !     !     $                          - 1.*kbeer*(
  !     !     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
  !     !     $                                         + pftpar(pft,46)
  !     !     $                                         )
  !     !     $                          )
  !     !          vpc_grid = vpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !     !          vpc_grid_total = vpc_grid_total + vpc_grid
  !     !
  !     !      ! Calculate local FCP treating dphen analogously as for the calulation of VAI:
  !     !      ! FPC = 1-exp(-kbeer*LAI*dphen) instead of FPC = dphen*(1-exp(-kbeer*LAI))
  !     !!           fpc_ind = 1. - dexp(
  !     !!     $                           -1.*kbeer*(
  !     !!     $                                         lai_ind(pft,jpngr)*dphen(day,pft)
  !     !!     $                                         )
  !     !!     $                           )
  !     !          fpc_ind = (1. - dexp(
  !     !     $                           -1.*kbeer*(
  !     !     $                                         lai_ind(pft,jpngr)
  !     !     $                                         )
  !     !     $                           ))!*dphen(day,pft)
  !     !          local_fpc_grid = fpc_ind * crownarea(pft,jpngr) * nind(pft,jpngr)
  !     !          fpc_grid_total = fpc_grid_total + local_fpc_grid
  !     !
  !     !          print*,'pft',pft
  !     !          print*,'local_fpc_grid     ',local_fpc_grid
  !     !          print*,'fpc_grid(pft,jpngr)',fpc_grid(pft,jpngr)
  !     !          
  !     !      ! Calculate fractional stem/branch cover of grid cell as the difference
  !     !          spc_grid = vpc_grid - local_fpc_grid
  !     !          spc_grid_total = spc_grid_total + spc_grid
  !    
  !   endif
  ! enddo
  !
  ! 
  ! if (vpc_grid_total.gt.1.) then
  !   !        print*,'-----------------scaling-------------------'
  !   !        print*,'fpc_grid_total',fpc_grid_total
  !   !        print*,'vpc_grid_total',vpc_grid_total
  !   scale = 1. / vpc_grid_total
  !   fpc_grid_total = fpc_grid_total * scale
  !   vpc_grid_total = vpc_grid_total * scale
  !   spc_grid_total = spc_grid_total * scale
  !   !        print*,'fpc_grid_total',fpc_grid_total
  !   !        print*,'vpc_grid_total',vpc_grid_total
  ! endif
  !
  ! if (fpc_grid_total.gt.1.) then
  !   !        print*,'fpc_grid_total',fpc_grid_total
  !   stop 
  ! endif
  !
  ! ! Daily N fixed by cryptogamic ground and plant covers (communicated to calling SR)
  ! !-------------------------------------------------------------------------
  ! ! Fixation scales with daily photosynthetically active radiation and the
  ! ! branch/stem surface for CPC and the bare ground surface for CGC.
  ! 
  ! dnfix_cpc = par_day(day) * max( 0., spc_grid_total) / glob_CPC_scal
  ! dnfix_cgc = par_day(day) * max( 0., (1.0 - vpc_grid_total) ) / glob_CGC_scal
  !
  ! end subroutine n_fixation_cryptogam


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


end module md_nuptake
