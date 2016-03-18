module _allocation
  !////////////////////////////////////////////////////////////////
  ! ALLOCATION MODULE
  ! Contains the "main" subroutine 'allocation_daily' and all 
  ! necessary subroutines for handling input/output, and auxiliary
  ! subroutines.
  ! Every module that implements 'allocation_daily' must contain 
  ! this list of subroutines (names that way).
  !   - allocation_daily
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _classdefs

  use _params_core, only: npft, nlu, maxgrid, ndaymonth, ndayyear, &
    c_molmass, n_molmass, nmonth

  implicit none

  private 
  public allocation_daily

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------


  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! output variables
  real, dimension(npft,maxgrid) :: outaCalclm
  real, dimension(npft,maxgrid) :: outaNalclm
  real, dimension(npft,maxgrid) :: outaCalcrm
  real, dimension(npft,maxgrid) :: outaNalcrm

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type paramstype_alloc
    ! Parameters determining the relationship of structural N and C to metabolic N
    ! From regressing Narea to metabolic Narea in Hikosaka data
    ! real, parameter    :: r_n_cw_v = 1.23223            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    ! real, parameter    :: ncw_min = 0.056               ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    real :: r_n_cw_v          = 0.1       ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    real :: ncw_min           = 0.1       ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    real :: r_ctostructn_leaf = 40        ! constant ratio of C to structural N
  end type paramstype_alloc

  type( paramstype_alloc ) :: params_alloc

  real, parameter :: frac_shoot = 0.5

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------
  type statetype_mustbe_zero_for_lai
    real :: cleaf
    real :: maxnv
  end type statetype_mustbe_zero_for_lai

  type leaftraits_type
    real :: narea
    real :: narea_metabolic
    real :: narea_structural
    real :: lma
    real :: nmass
    real :: r_cton_leaf
    real :: r_ntoc_leaf
  end type leaftraits_type

  ! states area global within module (instead of being passed on as arguments)
  type(statetype_mustbe_zero_for_lai) :: state_mustbe_zero_for_lai

  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dnroot


contains

  subroutine allocation_daily( jpngr, doy, moy, dom )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use _plant, only: params_plant, params_pft_plant, &
      pleaf, proot, plabl, r_cton_leaf, &
      r_ntoc_leaf, crownarea, drauto, dnpp, lai_ind, nind, fapar_ind, fpc_grid, &
      narea, narea_metabolic, narea_structural, lma, nmass, update_fpc_grid
    use _waterbal, only: solar
    use _gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use _findroot_fzeroin

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy     ! day of year
    integer, intent(in) :: moy     ! month of year
    integer, intent(in) :: dom     ! day of month

    ! local variables
    integer :: lu
    integer :: pft
    integer :: usemoy
    integer :: usedoy
    real    :: dclabl
    real    :: dnlabl
    logical :: cont          ! true if allocation to leaves (roots) is not 100% and not 0%
    real    :: max_dc
    real    :: min_dc
    real    :: abserr
    real    :: relerr
    real    :: nleaf0
    real    :: lai0, lai1
    integer, parameter :: nmax = 100

    type(outtype_zeroin)  :: out_zeroin
    type(leaftraits_type) :: traits

    integer, save      :: invocation = 0             ! internally counted simulation year
    integer, parameter :: spinupyr_phaseinit_2 = 1   ! this is unnecessary: might as well do flexible allocation right from the start.

    ! xxx debug
    real    :: test

    abserr=100.0*XMACHEPS !*10e5
    relerr=1000.0*XMACHEPS !*10e5

    !-------------------------------------------------------------------------
    ! Count number of calls (one for each simulation year) and allow flexible
    ! allocation only after year 'spinupyr_phaseinit_2'.
    !-------------------------------------------------------------------------
    if (doy==1) then
      invocation = invocation + 1
      ! write(0,*) 'WARNING: FIXED ALLOCATION'
    end if

    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      if (params_pft_plant(pft)%grass) then

        ! write(0,*) '--- allocation_daily, doy:',doy
        ! write(0,*) 'plabl(lu,jpngr)', plabl
        ! write(0,*) 'cleaf          ', pleaf(pft,jpngr)%c%c12
        ! write(0,*) 'croot          ', proot(pft,jpngr)%c%c12
        ! write(0,*) 'C:N in leaves  ', cton( pleaf(pft,jpngr), default=0.0 )

        if ( plabl(pft,jpngr)%c%c12>0.0 ) then
          !------------------------------------------------------------------
          ! Prescribe constant root:shoot ratio (in terms of C mass).
          ! Calculate maximum C allocatable based on current labile pool size.
          ! Maximum is the lower of all labile C and the C to be matched by 
          ! all labile N, discounted by the yield factor.
          !------------------------------------------------------------------
          if (pleaf(pft,jpngr)%c%c12==0.0) then
            ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
            r_cton_leaf(pft,jpngr) = get_rcton_init( solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            r_ntoc_leaf(pft,jpngr) = 1.0 / r_cton_leaf(pft,jpngr)
          end if
          max_dc = params_plant%growtheff * plabl(pft,jpngr)%c%c12

          !-------------------------------------------------------------------
          ! LEAF ALLOCATION
          !-------------------------------------------------------------------
          ! allocation to leaves is prescribed
          dcleaf(pft) = max_dc * frac_shoot

          call allocate_leaf( &
            dcleaf(pft), &
            pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
            plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
            solar%meanmppfd(:), mactnv_unitiabs(pft,:), &
            lai_ind(pft,jpngr), dnleaf(pft) &
            )

          !-------------------------------------------------------------------  
          ! Update leaf traits
          !-------------------------------------------------------------------  
          traits                      = get_leaftraits( lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          narea(pft)                  = traits%narea
          narea_metabolic(pft)        = traits%narea_metabolic
          narea_structural(pft)       = traits%narea_structural
          lma(pft,jpngr)              = traits%lma
          nmass(pft)                  = traits%nmass

          r_cton_leaf(pft,jpngr)      = traits%r_cton_leaf
          r_ntoc_leaf(pft,jpngr)      = traits%r_ntoc_leaf

          !-------------------------------------------------------------------  
          ! Update fpc_grid and fapar_ind (not lai_ind)
          !-------------------------------------------------------------------  
          call update_fpc_grid( pft, jpngr )

          !-------------------------------------------------------------------
          ! ROOT ALLOCATION
          !-------------------------------------------------------------------
          call allocate_root( &
            proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
            plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
            pft, dcroot(pft), dnroot(pft) &
            )
        
          !-------------------------------------------------------------------
          ! GROWTH RESPIRATION
          !-------------------------------------------------------------------
          ! add growth respiration to autotrophic respiration and substract from NPP
          ! (note that NPP is added to plabl in and growth resp. is implicitly removed
          ! from plabl above)
          drauto(pft)   = drauto(pft)     + ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) )
          dnpp(pft)%c12 = dnpp(pft)%c12   - ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) 

        else

          dcleaf(pft) = 0.0
          dcroot(pft) = 0.0
          dnleaf(pft) = 0.0
          dnroot(pft) = 0.0

          ! write(0,*) 'not growing ...'

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! if (doy==39) stop 'on day 39'

    ! test_calloc = test_calloc + dcleaf + dcroot
    ! write(0,*) 'test_calloc', test_calloc
    
    ! write(0,*) '--- END allocation_daily:'

  end subroutine allocation_daily


  subroutine allocate_leaf( mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd,    nv,    lai,   mydnleaf )
    !///////////////////////////////////////////////////////////////////
    ! LEAF ALLOCATION
    ! Sequence of steps:
    ! - increment foliage C pool
    ! - update LAI
    ! - calculate canopy-level foliage N as a function of LAI 
    ! - reduce labile pool by C and N increments
    !-------------------------------------------------------------------
    use _classdefs
    use _plant, only: params_plant

    ! arguments
    real, intent(in)                    :: mydcleaf
    real, intent(inout)                 :: cleaf, nleaf
    real, intent(inout)                 :: clabl, nlabl
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(out)                   :: lai
    real, optional, intent(out)         :: mydnleaf

    ! local variables
    real :: nleaf0
    real :: dclabl

    ! find LAI, given new leaf mass. This is necessary to get leaf-N as 
    ! a function of LAI.
    if (mydcleaf>0.0) then

      ! write(0,*) 'cleaf before ', cleaf 
      cleaf  = cleaf + mydcleaf
      ! write(0,*) 'mydcleaf     ', mydcleaf 

      ! write(0,*) 'cleaf = 0.5', get_lai( 0.5  , meanmppfd(:), nv(:) )
      ! write(0,*) 'cleaf = 100', get_lai( 100.0, meanmppfd(:), nv(:))
      ! stop

      ! Calculate LAI as a function of leaf C
      lai = get_lai( cleaf, meanmppfd(:), nv(:) )
      ! write(0,*) 'in allocate_leaf: cleaf, lai :', cleaf, lai
      ! write(0,*) 'mydcleaf ', mydcleaf 
      ! ! stop

      ! calculate canopy-level leaf N as a function of LAI
      nleaf0   = nleaf      
      nleaf    = get_canopy_leaf_n( lai, meanmppfd(:), nv(:) )
      mydnleaf = nleaf - nleaf0
      ! if (mydnleaf>0.0) then
      !   write(0,*) 'mydcleaf/dnleaf ', mydcleaf/mydnleaf 
      ! end if

      ! subtract from labile pool, making sure pool does not get negative
      dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcleaf )
      if ( (dclabl - clabl) > 1e-8 ) stop 'trying to remove too much from labile pool: leaf C'
      clabl  = clabl - dclabl

      ! write(0,*) 'params_pft_plant(pft)%r_ntoc_root  ', params_pft_plant(pft)%r_ntoc_root

    else

      lai      =  get_lai( cleaf, meanmppfd(:), nv(:) )
      mydnleaf = 0.0

    end if

  end subroutine allocate_leaf


  subroutine allocate_root( croot, nroot, clabl, nlabl, pft, mydcroot, mydnroot )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use _classdefs
    use _plant, only: params_plant, params_pft_plant

    ! arguments
    real, intent(inout)         :: croot, nroot
    real, intent(inout)         :: clabl, nlabl
    integer, intent(in)         :: pft
    real, optional, intent(out) :: mydcroot
    real, optional, intent(out) :: mydnroot

    ! local variables
    real :: dclabl

    if (clabl>0.0) then
      ! use remainder for allocation to roots
      mydcroot = params_plant%growtheff * clabl
      mydnroot = mydcroot * params_pft_plant(pft)%r_ntoc_root

      dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcroot )
      if ( (dclabl - clabl) > 1e-8 ) stop 'trying to remove too much from labile pool: root C'
      clabl  = clabl - dclabl

      if (mydcroot<0.0) stop 'root allocation neg.: C'
      if (mydnroot<0.0) stop 'root allocation neg.: N'

      croot = croot + mydcroot
      nroot = nroot + mydnroot

      ! if (present(verbose)) then
      !   write(0,*) 'after presumed allocation:'
      !   write(0,*) 'clabl', clabl ! OK
      !   write(0,*) 'nlabl', nlabl ! OK
      !   ! write(0,*) 'mydcleaf', mydcleaf
      !   ! write(0,*) 'mydcroot', mydcroot
      !   ! write(0,*) 'mydnleaf', mydnleaf
      !   ! write(0,*) 'mydnroot', mydnroot

      !   ! write(0,*) 'after presumed allocation:'
      !   ! write(0,*) 'B cleaf', cleaf
      !   ! write(0,*) 'B root', root
      ! end if
    end if

  end subroutine allocate_root


  function get_lai( cleaf, meanmppfd, nv ) result( lai )
    !////////////////////////////////////////////////////////////////
    ! Calculates LAI as a function of canopy-level leaf-C:
    ! Cleaf = Mc * c * ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) )
    ! Cannot be solved analytically for L = f(Cleaf). Therefore, 
    ! numerical root-searching algorithm is applied so that
    ! Cleaf / ( Mc * c ) - ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) ) = 0
    ! This is implemented in function 'mustbe_zero_for_lai()'.
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth
    use _findroot_fzeroin

    ! arguments 
    real, intent(in)                    :: cleaf
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv 

    ! local variables
    real                 :: abserr
    real                 :: relerr
    real                 :: lower
    real                 :: upper
    integer, parameter   :: nmax = 100
    type(outtype_zeroin) :: out_zeroin

    ! xxx debug
    real :: test

    ! function return value
    real, intent(out) :: lai

    ! local variables
    real :: maxnv

    if (cleaf>0.0) then
      ! Metabolic N is predicted and is optimised at a monthly time scale. 
      ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
      ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost.
      ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
      ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
      maxnv = maxval( meanmppfd(:) * nv(:) )
      ! write(0,*) 'maxnv ', maxnv
      ! stop

      ! write(0,*) 'what is LAI for Cleaf=', cleaf
      ! write(0,*) 'maxnv ', maxnv

      ! Update state. This derived-type variable is "global" within this module
      state_mustbe_zero_for_lai%cleaf = cleaf
      state_mustbe_zero_for_lai%maxnv = maxnv

      ! Calculate initial guess for LAI (always larger than actual LAI)
      lower = 0.0 ! uninformed lower bound
      upper = 1.0 / ( c_molmass * params_alloc%ncw_min * params_alloc%r_ctostructn_leaf ) * cleaf
      ! upper = 20.0
      ! write(0,*) 'upper ', upper

      ! write(0,*) 'lower =', lower
      ! test = mustbe_zero_for_lai( lower ) 
      ! write(0,*) '=> test =', test

      ! write(0,*) 'upper =', upper
      ! test = mustbe_zero_for_lai( upper ) 
      ! write(0,*) '=> test =', test

      ! call function zeroin to find root (value of LAI for which evaluation expression is zero)
      abserr=100.0*XMACHEPS !*10e5
      relerr=1000.0*XMACHEPS !*10e5

      ! write(0,*) 'abserr', abserr
      ! write(0,*) 'relerr', relerr
      ! stop 'here'

      ! write(0,*) '*** finding root of mustbe_zero_for_lai ***'
      out_zeroin = zeroin( mustbe_zero_for_lai, abserr, relerr, nmax, lower, upper )
      if ( out_zeroin%error /= 0 ) then
        lai = 0.0
        write(0,*) 'error code', out_zeroin%error
        stop 'zeroin for mustbe_zero_for_lai() failed'
      else
        lai = out_zeroin%root
      end if

      ! write(0,*) 'out_zeroin', out_zeroin
      ! write(0,*) 'cleaf', cleaf
      ! write(0,*) 'lai', lai
      ! stop
    
    else

      lai = 0.0

    end if

  end function get_lai


  function mustbe_zero_for_lai( mylai ) result( mustbe_zero )
    !/////////////////////////////////////////////////////////
    ! This function returns value of the expression 'mustbe_zero'. 
    ! If expression is zero, then 'mylai' is a root and is the 
    ! LAI for a given Cleaf (meanmppfd, cleaf, and nv are 
    ! passed on to this function as a derived type state.)
    !---------------------------------------------------------
    ! Cleaf = c_molmass * params_alloc%r_ctostructn_leaf * N_canop_cellwall
    ! N_canop_cellwall = LAI * params_alloc%ncw_min + nv * Iabs * params_alloc%r_n_cw_v
    ! Iabs = meanmppfd * (1-exp( -kbeer * LAI))
    ! ==> Cleaf = f(LAI) = c_molmass * params_alloc%r_ctostructn_leaf * [ meanmppfd * (1-exp( -kbeer * LAI)) * nv * params_alloc%r_n_cw_v + LAI * params_alloc%ncw_min ]
    ! ==> LAI = f(Cleaf) leads to inhomogenous equation. Therefore apply root finding algorithm so that:
    ! 0 = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - meanmppfd * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * nv * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    !---------------------------------------------------------
    use _plant, only: params_plant, get_fapar

    ! arguments
    real, intent(in) :: mylai

    ! function return value
    real, intent(out) :: mustbe_zero

    ! local variables
    real :: mycleaf
    real :: mymaxnv

    ! write(0,*) '--- in mustbe_zero_for_lai with mydcleaf=', mylai

    ! Read from updated state. This derived-type variable is "global" within this module
    mycleaf = state_mustbe_zero_for_lai%cleaf
    mymaxnv = state_mustbe_zero_for_lai%maxnv

    ! write(0,*) '----------'
    ! write(0,*) 'inside mustbe_zero_for_lai: '
    ! write(0,*) 'mylai', mylai
    ! write(0,*) 'mycleaf', mycleaf
    ! write(0,*) 'mymaxnv', mymaxnv

    ! mustbe_zero = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - meanmppfd * nv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    ! mustbe_zero = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - maxnv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    mustbe_zero = mycleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - mymaxnv * get_fapar( mylai ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min

    ! write(0,*) 'mustbe_zero                           ', mustbe_zero
    ! write(0,*) '-------------'


  end function mustbe_zero_for_lai


  function get_rcton_init( meanmppfd, nv ) result( rcton )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Cleaf = c_molmass * params_alloc%r_ctostructn_leaf * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * params_alloc%r_n_cw_v + LAI * params_alloc%ncw_min ]
    ! Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (params_alloc%r_n_cw_v + 1) + LAI * params_alloc%ncw_min ]
    ! linearization around LAI = 0 ==> (1-exp(-k*L)) ~= k*L
    ! ==> Cleaf ~= LAI * c_molmass * params_alloc%r_ctostructn_leaf * ( meanmppfd * kbeer * nv * params_alloc%r_n_cw_v + params_alloc%ncw_min )
    ! ==> Nleaf ~= LAI * n_molmass * ( meanmppfd * kbeer * nv * (params_alloc%r_n_cw_v + 1) + params_alloc%ncw_min )
    ! r_cton = Cleaf / Nleaf
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth
    use _plant, only: params_plant

    ! arguments
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: rcton

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )
    rcton = ( c_molmass * params_alloc%r_ctostructn_leaf * ( maxnv * params_plant%kbeer * params_alloc%r_n_cw_v + params_alloc%ncw_min )) / ( n_molmass * ( maxnv * params_plant%kbeer * (params_alloc%r_n_cw_v + 1.0) + params_alloc%ncw_min ) )

    ! rcton = ( c_molmass * params_alloc%r_ctostructn_leaf * ( meanmppfd * kbeer * nv * params_alloc%r_n_cw_v + params_alloc%ncw_min )) / ( n_molmass * ( meanmppfd * kbeer * nv * (params_alloc%r_n_cw_v + 1) + params_alloc%ncw_min ) )

  end function get_rcton_init


  function get_canopy_leaf_n_metabolic( mylai, meanmppfd, nv ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth
    use _plant, only: get_fapar

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: mynleaf_metabolic  ! mol N 

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    mynleaf_metabolic = maxnv * get_fapar( mylai )

  end function get_canopy_leaf_n_metabolic


  function get_canopy_leaf_n_structural( mylai, mynleaf_metabolic ) result( mynleaf_structural )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_structural = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth

    ! arguments
    real, intent(in) :: mylai
    real, intent(in) :: mynleaf_metabolic

    ! function return variable
    real, intent(out) :: mynleaf_structural  ! mol N 

    mynleaf_structural = mynleaf_metabolic * params_alloc%r_n_cw_v + mylai * params_alloc%ncw_min

    ! write(0,*) '--- in get_canopy_leaf_n_structural'
    ! write(0,*) 'mylai ', mylai
    ! write(0,*) 'mynleaf_metabolic ', mynleaf_metabolic
    ! write(0,*) 'params_alloc%r_n_cw_v ', params_alloc%r_n_cw_v
    ! write(0,*) 'params_alloc%ncw_min ', params_alloc%ncw_min
    ! write(0,*) 'mynleaf_structural ', mynleaf_structural
    ! write(0,*) '-------------------------------'

  end function get_canopy_leaf_n_structural


  function get_canopy_leaf_n( mylai, meanmppfd, nv ) result( mynleaf )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Nleaf = LAI * (n_metabolic + n_cellwall) * n_molmass
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    ! LAI * n_cellwall = LAI * (params_alloc%ncw_min + params_alloc%r_n_cw_v * n_metabolic)
    ! ==> Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (params_alloc%r_n_cw_v + 1) + LAI * params_alloc%ncw_min ]
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: mynleaf ! g N

    ! local variables
    real :: nleaf_metabolic   ! mol N m-2
    real :: nleaf_structural  ! mol N m-2

    nleaf_metabolic  = get_canopy_leaf_n_metabolic(  mylai, meanmppfd, nv )
    nleaf_structural = get_canopy_leaf_n_structural( mylai, nleaf_metabolic )
    mynleaf          = n_molmass * ( nleaf_metabolic + nleaf_structural )

    ! write(0,*) '--- in get_canopy_leaf_n'
    ! write(0,*) 'nleaf_metabolic ', nleaf_metabolic
    ! write(0,*) 'nleaf_structural ', nleaf_structural
    ! write(0,*) 'mynleaf ', mynleaf
    ! write(0,*) '-------------------------------'

    ! mynleaf = n_molmass * ( maxnv * get_fapar( mylai ) * ( 1.0 + params_alloc%r_n_cw_v ) + mylai * params_alloc%ncw_min )

  end function get_canopy_leaf_n


  function get_leaftraits( mylai, meanmppfd, nv ) result( traits )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = a + b * Narea_metabolic
    ! Carea            = c * Narea_structural
    !----------------------------------------------------------------
    use _params_core, only: c_content_of_biomass
    use _plant, only: get_fapar

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    type(leaftraits_type) :: traits

    ! local variables
    real :: mynarea_metabolic_canop   ! mol N m-2-ground
    real :: mynarea_structural_canop  ! mol N m-2-ground

    mynarea_metabolic_canop  = get_canopy_leaf_n_metabolic(  mylai, meanmppfd(:), nv(:) )     ! mol N m-2-ground    
    mynarea_structural_canop = get_canopy_leaf_n_structural( mylai, mynarea_metabolic_canop ) ! mol N m-2-ground
    
    traits%narea_metabolic  = n_molmass * mynarea_metabolic_canop / mylai   ! g N m-2-leaf
    traits%narea_structural = n_molmass * mynarea_structural_canop / mylai  ! g N m-2-leaf
    traits%narea            = n_molmass * ( mynarea_metabolic_canop + mynarea_structural_canop ) / mylai ! g N m-2-leaf
    traits%lma              = c_molmass * params_alloc%r_ctostructn_leaf * mynarea_structural_canop / mylai 
    traits%nmass            = traits%narea / ( traits%lma / c_content_of_biomass )
    traits%r_cton_leaf      = traits%lma / traits%narea
    traits%r_ntoc_leaf      = 1.0 / traits%r_cton_leaf

    ! write(0,*) '--- in get_leaftraits'
    ! write(0,*) 'mylai                  ', mylai
    ! write(0,*) 'traits%narea_metabolic ', traits%narea_metabolic 
    ! write(0,*) 'traits%narea_structural', traits%narea_structural
    ! write(0,*) 'traits%narea           ', traits%narea           
    ! write(0,*) 'traits%lma             ', traits%lma             
    ! write(0,*) 'traits%nmass           ', traits%nmass           
    ! write(0,*) 'traits%r_cton_leaf     ', traits%r_cton_leaf     
    ! write(0,*) 'traits%r_ntoc_leaf     ', traits%r_ntoc_leaf     
    ! write(0,*) '-------------------------------'
    ! stop 

  end function get_leaftraits


  subroutine initio_allocation()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use _params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------

    ! C ALLOCATED TO LEAF GROWTH 
    filnam=trim(prefix)//'.a.calclm.out'
    open(350,file=filnam,err=999,status='unknown')

    ! N ALLOCATED TO LEAF GROWTH 
    filnam=trim(prefix)//'.a.nalclm.out'
    open(351,file=filnam,err=999,status='unknown')

    ! C ALLOCATED TO ROOT GROWTH 
    filnam=trim(prefix)//'.a.calcrm.out'
    open(352,file=filnam,err=999,status='unknown')

    ! N ALLOCATED TO ROOT GROWTH 
    filnam=trim(prefix)//'.a.nalcrm.out'
    open(353,file=filnam,err=999,status='unknown')

    return

    999  stop 'INITIO_ALLOCATION: error opening output files'

  end subroutine initio_allocation


  subroutine initoutput_allocation()
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    ! xxx remove their day-dimension
    outaCalclm(:,:) = 0.0
    outaNalclm(:,:) = 0.0
    outaCalcrm(:,:) = 0.0
    outaNalcrm(:,:) = 0.0

    ! write(0,*) 'initialising outaCalloc',outaCalloc

  end subroutine initoutput_allocation


  subroutine getout_daily_allocation( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    outaCalclm(:,jpngr) = outaCalclm(:,jpngr) + dcleaf(:) 
    outaNalclm(:,jpngr) = outaNalclm(:,jpngr) + dnleaf(:)
    outaCalcrm(:,jpngr) = outaCalcrm(:,jpngr) + dcroot(:) 
    outaNalcrm(:,jpngr) = outaNalcrm(:,jpngr) + dnroot(:)

    ! write(0,*) 'collecting outaCalloc',outaCalloc

  end subroutine getout_daily_allocation


  subroutine writeout_ascii_allocation( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_siml, only: firstyeartrend, spinupyears

    ! arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! local variables
    real    :: itime
    integer :: jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(year) + real(firstyeartrend) - real(spinupyears)

    ! write(0,*) 'writing time, outaCalloc',itime, sum(outaCalloc(:,jpngr))

    write(350,999) itime, sum(outaCalclm(:,jpngr))
    write(351,999) itime, sum(outaNalclm(:,jpngr))
    write(352,999) itime, sum(outaCalcrm(:,jpngr))
    write(353,999) itime, sum(outaNalcrm(:,jpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module _allocation
