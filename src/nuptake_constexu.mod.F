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
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: ndayyear, nmonth, nlu, npft, maxgrid
  use _vars_core, only: dcex, pninorg, plabl, proot, dcex, dnup, &
      dtransp, aet_d, fpc_grid

  use _classdefs

  implicit none

  ! MODULE-SPECIFIC PARAMETERS
  real :: exurate           ! C exudation rate per unit root mass
  real :: eff_nup           ! uptake efficiency for equation
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

  ! OUTPUT VARIABLES
  ! daily
  real, dimension(npft,ndayyear,maxgrid) :: outdccost   ! daily mean C cost of N uptake (gC/gN) 
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_pas
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_act
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_fix
  real, dimension(npft,ndayyear,maxgrid) :: outdnup_ret
  real, dimension(npft,ndayyear,maxgrid) :: outdcex

  ! annual
  real, dimension(npft,maxgrid) :: outacex

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

    ! arguments
    integer, intent(in) :: jpngr, pft

    ! local variables
    integer :: lu
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

    ! xxx try: 
    ! write(0,*) 'NUPTAKE: pninorg pool size set manually'
    ! pninorg%n14 = 0.01

    if ( proot(pft,jpngr)%c%c12 == 0.0 ) then

      ! write(0,*) 'no roots'

      ! no roots => no N uptake, no C exudation
      n_uptake_pass    = 0.0
      dNacq_active     = 0.0
      dNacq_fix        = 0.0
      n_uptake_retrans = 0.0

    else

      ! ! xxx try:
      ! dtransp = aet_d(lu_category(pft)) * fpc_grid(pft,jpngr)

      !//////////////////////////////////////////////////////////////////////////
      ! INITIALIZATION
      !-------------------------------------------------------------------------
      n_uptake_pass   = 0.0

      dNacq_active    = 0.0                          ! active uptake, sum over sub-timesteps
      dNacq_fix       = 0.0                          ! N fixation, sum over sub-timesteps
      dCacq_active    = 0.0
      dCacq_fix       = 0.0

      ! write(0,*) 'lu, jpngr',lu, jpngr
      ! write(0,*) 'pninorg(lu,jpngr)%n14',pninorg(lu,jpngr)%n14
      ! write(0,*) 'dwtot(lu,jpngr)      ',dwtot(lu,jpngr)

      ! ninorg_conc   = calc_conc_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   
      ! avail_ninorg  = calc_avail_ninorg( pninorg(lu,jpngr)%n14, dwtot(lu,jpngr) )   

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
      dCacq_active = exurate * proot(pft,jpngr)%c%c12
      ! dNacq_active = calc_dnup( dCacq_active, avail_ninorg )
      dNacq_active = calc_dnup( dCacq_active, pninorg(lu,jpngr)%n14 )

      ! write(0,*) 'dCacq_active ',dCacq_active 
      ! write(0,*) 'in SR nuptake: dCacq_active          ',dCacq_active 
      ! write(0,*) 'in SR nuptake: dNacq_active          ',dNacq_active 
      ! write(0,*) 'in SR nuptake: pninorg(lu,jpngr)%n14 ',pninorg(lu,jpngr)%n14 

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
      dccost(pft) = 0.0
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

    ! N uptake function of C exuded
    dnup = n0 * ( 1.0 - exp( - eff_nup * cexu ) )

  end function calc_dnup


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


  subroutine initio_nuptake()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use _params_siml, only: runname

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
    use _sofunutils, only: getparreal

    ! C exudation rate per unit root mass
    exurate = getparreal( 'params_nuptake_constexu.dat', 'exurate' )

    ! uptake efficiency for equation
    ! dCexu/dNup = K / (N0 - Nup); K=1/eff_nup
    eff_nup = getparreal( 'params_nuptake_constexu.dat', 'eff_nup' )

    ! shape parameter of cost function of N fixation 
    ! Below parameters (MINIMUMCOSTFIX, FIXOPTIMUM, FIXWIDTH ) are based on 
    ! the assumption that the cost of symbiotic N fixation is the 
    ! inverse of nitrogenase activity. 
    ! After Houlton et al., 2008. Minimum cost of N-fixation is 4.8 gC/gN
    ! (value from Gutschik 1981)
    MINIMUMCOSTFIX = getparreal( 'params_nuptake_constexu.dat', 'MINIMUMCOSTFIX' )

    ! shape parameter of cost function of N fixation 
    FIXOPTIMUM = getparreal( 'params_nuptake_constexu.dat', 'FIXOPTIMUM' )
 
    ! shape parameter of cost function of N fixation 
    FIXWIDTH = getparreal( 'params_nuptake_constexu.dat', 'FIXWIDTH' )


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
    outdcex(:,:,:)     = 0.0

    outacex(:,:)       = 0.0


  end subroutine initoutput_nuptake



  subroutine getout_daily_nuptake( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
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


  subroutine writeout_ascii_nuptake( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_core, only: ndayyear, npft, nlu
    use _params_siml, only: firstyeartrend, spinupyears

    ! arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! local variables
    real :: itime
    integer :: myday, myjpngr

    real :: myoutdccost
    real :: myoutdnup_pas
    real :: myoutdnup_act
    real :: myoutdnup_fix
    real :: myoutdnup_ret

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    myjpngr = 1


    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (.not.spinup) then
      ! Write daily output only during transient simulation
      do myday=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(myday-1)/real(ndayyear)

        if (nlu>1) stop 'writeout_ascii_nuptake: write out lu-area weighted sum'
        if (npft>1) stop 'writeout_ascii_nuptake: think of something for ccost output'

        ! xxx lu-area weighted sum if npft>0
        myoutdccost   = sum(outdccost(:,myday,myjpngr))
        myoutdccost   = myoutdccost / real(npft)
        myoutdnup_pas = sum(outdnup_pas(:,myday,myjpngr))
        myoutdnup_act = sum(outdnup_act(:,myday,myjpngr))
        myoutdnup_fix = sum(outdnup_fix(:,myday,myjpngr))
        myoutdnup_ret = sum(outdnup_ret(:,myday,myjpngr))

        write(105,999) itime, sum(outdcex(:,myday,myjpngr))
        write(400,999) itime, myoutdccost
        write(401,999) itime, myoutdnup_pas
        write(402,999) itime, myoutdnup_act
        write(403,999) itime, myoutdnup_fix
        write(404,999) itime, myoutdnup_ret

      end do
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(year) + real(firstyeartrend) - real(spinupyears)

    write(405,999) itime, sum(outacex(:,myjpngr))

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


end module _nuptake
