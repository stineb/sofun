module _nuptake
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
  !----------------------------------------------------------------
  use _params_core, only: ndayyear, nmonth, nlu, npft, maxgrid
  use _classdefs

  implicit none

  ! MODULE-SPECIFIC PARAMETERS
  real :: exurate           ! C exudation rate per unit root mass
  real :: eff_nup           ! uptake efficiency for equation
  real :: EPSILON_WTOT      ! minimum soil water content above which N becomes available  
  real :: MINIMUMCOSTFIX    ! minimum cost of N fixation (at optimal temperature)
  real :: FIXOPTIMUM        ! optimum temperature for N fixation
  real :: FIXWIDTH          ! shape parameter for width of N fixation cost function

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
    ! SUBROUTINE NUPTAKE ASSUMING CONSTANT EXUDATION PER UNIT ROOT MASS
    !-----------------------------------------------------------------
    ! This model calculates first the passive uptake of N via the 
    ! transpiration stream.
    !-----------------------------------------------------------------
    use _classdefs
    use _params_modl, only: lu_category
    use _vars_core, only: pninorg, plabl, proot, dwtot
    use _vars_core, only: dcex, dnup
    use _vars_core, only: aet_d
    use _vars_core, only: fpc_grid

    ! arguments
    integer, intent(in) :: jpngr, pft

    ! local variables
    integer :: lu
    real    :: dtransp                             ! daily transpiration for this PFT (mm)
    real    :: avail_ninorg                        ! available inorganic N in soil layer (gN/m2)
    real    :: ninorg_conc                         ! inorganic N concentration (gN/gH2O)
    real    :: n_uptake_pass                       ! (gN)
    real    :: n_uptake_retrans
    real    :: dNacq_active
    real    :: dNacq_fix
    real    :: dCacq_active
    real    :: dCacq_fix
    real    :: dmean_cost

    lu = lu_category(pft)

    ! xxx think about this: order in which PFTs get access to Ninorg matters!
    
    ! write(0,*) '---- in nuptake:'

    if ( proot(pft,jpngr)%c%c12 == 0.0 ) then

      ! write(0,*) 'no roots'

      ! no roots => no N uptake, no C exudation
      n_uptake_pass    = 0.0
      dNacq_active     = 0.0
      dNacq_fix        = 0.0
      n_uptake_retrans = 0.0

    else

      ! xxx try:
      dtransp = aet_d(lu_category(pft)) * fpc_grid(pft,jpngr)

      !//////////////////////////////////////////////////////////////////////////
      ! INITIALIZATION
      !-------------------------------------------------------------------------
      dNacq_active    = 0.0                          ! active uptake, sum over sub-timesteps
      dNacq_fix       = 0.0                          ! N fixation, sum over sub-timesteps
      dCacq_active    = 0.0
      dCacq_fix       = 0.0

      ! write(0,*) 'lu, jpngr',lu, jpngr
      ! write(0,*) 'pninorg(lu,jpngr)%n14',pninorg(lu,jpngr)%n14
      ! write(0,*) 'dwtot(lu,jpngr)      ',dwtot(lu,jpngr)

      ninorg_conc   = calc_conc_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   
      avail_ninorg  = calc_avail_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   

      ! write(0,*) 'avail_ninorg',avail_ninorg
      ! write(0,*) 'ninorg_conc ',ninorg_conc 
      ! stop

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


      !//////////////////////////////////////////////////////////////////////////
      ! PASSIVE UPTAKE
      ! No active control on passive uptake - always occurrs even if the unmet N
      ! demand is zero.
      !--------------------------------------------------------------------------
      n_uptake_pass = ninorg_conc * dtransp    
      n_uptake_pass = min( n_uptake_pass, avail_ninorg )

      ! write(0,*) 'n_uptake_pass ',n_uptake_pass 

      ! Update
      pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - n_uptake_pass
      avail_ninorg          = calc_avail_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   
      ninorg_conc           = calc_conc_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   

      ! write(0,*) 'avail_ninorg',avail_ninorg
      ! write(0,*) 'ninorg_conc ',ninorg_conc 

      !//////////////////////////////////////////////////////////////////////////
      ! ACTIVE UPTAKE
      ! Active N uptake is a function of initial N available and C exuded
      !--------------------------------------------------------------------------         
      dCacq_active = exurate * proot(pft,jpngr)%c%c12
      dNacq_active = calc_dnup( dCacq_active, avail_ninorg )

      ! write(0,*) 'dCacq_active ',dCacq_active 
      ! write(0,*) 'dNacq_active ',dNacq_active 

      dcex(pft) = dcex(pft) + dCacq_active
      dmean_cost = dCacq_active / dNacq_active

      ! Update
      pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 - dNacq_active

    end if

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

    ! write(0,*) '---- finished nuptake'

  end subroutine nuptake


  function calc_dnup( cexu, n0 ) result( dnup )
    !/////////////////////////////////////////////////////////////////
    ! With a FUN-like approach:
    ! dCexu/dNup = K / (N0 - Nup); K=1/eff_nup
    ! => Nup(Cexu) = N0 * ( 1.0 - exp( - eff_nup * cexu ) )
    !-----------------------------------------------------------------
    ! arguments
    real :: cexu   ! C exuded (gC/m2/d)
    real :: n0     ! initial available N (gN/m2)

    ! function return variable
    real :: dnup

    dnup = n0 * ( 1.0 - exp( - eff_nup * cexu ) )

  end function calc_dnup


  function calc_avail_ninorg( ninorg, wtot ) result( avail_ninorg )
    !//////////////////////////////////////////////////////////////////////////
    ! Returns N available for uptake accounting for constraint of mobility by
    ! soil moisture.
    !--------------------------------------------------------------------------
    ! arguments
    real, intent(in)  :: ninorg 
    real, intent(in)  :: wtot           ! total soil water content (mm)

    ! function return value
    real, intent(out) :: avail_ninorg

    if ( wtot > EPSILON_WTOT ) then 
      avail_ninorg = ninorg - EPSILON_WTOT*ninorg/wtot
    else
      avail_ninorg = 0.0
    endif

  end function calc_avail_ninorg


  function calc_conc_ninorg( ninorg, wtot ) result( conc_ninorg )
    !//////////////////////////////////////////////////////////////////////////
    ! Returns N available for uptake accounting for constraint of mobility by
    ! soil moisture.
    !--------------------------------------------------------------------------
    ! arguments
    real, intent(in)  :: ninorg 
    real, intent(in)  :: wtot           ! total soil water content (mm)

    ! function return value
    real, intent(out) :: conc_ninorg

    if ( wtot > EPSILON_WTOT ) then 
      conc_ninorg = ninorg/wtot
    else
      conc_ninorg = 0.0
    endif

  end function calc_conc_ninorg


  function fun_cost_fix( soiltemp )
    !******************************************************************************
    ! Cost of symbiotic N fixation is the inverse of nitrogenase activity
    ! after Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    !--------------------------------------------------------------------------      
    use _params_modl

    real, intent(in) :: soiltemp
    real :: fun_cost_fix                 ! function return variable

    fun_cost_fix = MINIMUMCOSTFIX + exp((soiltemp-FIXOPTIMUM)**2/(2*FIXWIDTH**2))    ! inverse gauss function  (take WARMEST layer)

  end function fun_cost_fix
  

  subroutine initdaily_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Initialise daily variables with zero
    !----------------------------------------------------------------
    dnup_pas(:)    = 0.0
    dnup_act(:)    = 0.0
    dnup_fix(:)    = 0.0
    dnup_ret(:)    = 0.0

  end subroutine initdaily_nuptake


  subroutine initio_nuptake( prefix )
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    ! arguments
    character(len=*) :: prefix

    ! local variables
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



  subroutine getpar_modl_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    ! local variables
    integer, parameter :: npar = 6
    real, dimension(npar) :: params_array
    character(len=50)     ::paramfilnam
    
    paramfilnam = 'params_nuptake_constexu.dat'

    open(unit=04,file=trim(paramfilnam),status='OLD')      
    read (04,*) params_array
    close (04)

    ! C exudation rate per unit root mass
    exurate = params_array(1)

    ! uptake efficiency for equation
    ! dCexu/dNup = K / (N0 - Nup); K=1/eff_nup
    eff_nup = params_array(2)

    ! minimum soil water content above which N becomes available 
    ! central value 0.05 mm
    EPSILON_WTOT = params_array(3)

    ! shape parameter of cost function of N fixation 
    ! Below parameters (MINIMUMCOSTFIX, FIXOPTIMUM, FIXWIDTH ) are based on 
    ! the assumption that the cost of symbiotic N fixation is the 
    ! inverse of nitrogenase activity. 
    ! After Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    MINIMUMCOSTFIX = params_array(4)

    ! shape parameter of cost function of N fixation 
    FIXOPTIMUM = params_array(5)

    ! shape parameter of cost function of N fixation 
    FIXWIDTH = params_array(6)


  end subroutine getpar_modl_nuptake



  subroutine initoutput_nuptake
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------

    ! xxx remove their day-dimension
    outdccost(:,:,:)   = 0.0
    outdnup_pas(:,:,:) = 0.0
    outdnup_act(:,:,:) = 0.0
    outdnup_fix(:,:,:) = 0.0
    outdnup_ret(:,:,:) = 0.0

    return

  end subroutine initoutput_nuptake



  subroutine getout_daily_nuptake( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    ! Save the daily totals:
    ! xxx add lu-dimension and jpngr-dimension
    outdccost(:,doy,jpngr)   = dccost(:)
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
    use _params_core, only: ndayyear, nmonth, npft
    use _params_siml, only: firstyeartrend, spinupyears


    ! arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! local variables
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



end module _nuptake
