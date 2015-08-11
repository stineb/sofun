module _npp
  !////////////////////////////////////////////////////////////////
  ! NPP_LPJ MODULE
  ! Contains the "main" subroutine 'npp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'npp' must contain this list 
  ! of subroutines (names that way).
  !   - npp
  !   - getpar_modl_npp
  !   - initio_npp
  !   - initoutput_npp
  !   - getout_daily_npp
  !   - getout_monthly_npp
  !   - writeout_ascii_npp
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - daily NPP ('dnpp')
  !   - soil temperature ('xxx')
  !   - inorganic N _pools ('no3', 'nh4')
  !   - xxx 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  implicit none

  ! MODULE SPECIFIC PARAMETERS
  real :: r_root
  real :: r_sapw

contains

  subroutine npp( jpngr, dtemp, doy )
    !/////////////////////////////////////////////////////////////////////////
    ! NET PRIMARY PRODUCTIVITY
    ! Calculate maintenance and growth respiration and substract this from GPP 
    ! to get NPP before additional root respiration for nutrient uptake (see 
    ! SR nuptake). NPP is defined so as to include C allocated to growth as 
    ! well as to exudation. Thus, exudation is not part of autotrophic respir-
    ! ation, but diverted to the quickly decaying exudates pool. Exudates decay
    ! is calculated in SR 'littersom' and is kept track of as soil respiration 
    ! ('rsoil'). This implies that growth respiration is "paid" also on exu-
    ! dates. 
    !-------------------------------------------------------------------------
    use _classdefs
    use _params_core, only: npft, ndayyear
    use _params_modl, only: lu_category, tree, grass
    use _rates
    use _vars_core, only: dgpp, dnpp, drauto, drleaf, drroot, drsapw, dcex, dnup
    use _vars_core, only: nind, ispresent
    use _vars_core, only: proot, psapw, pexud, plabl
    use _nuptake, only: nuptake
    use _vars_core, only: dtemp_soil
    use _gpp, only: drd
    use _phenology, only: summergreen, shedleaves

    ! arguments
    integer, intent(in) :: jpngr
    real, intent(in)    :: dtemp      ! air temperature at this day

    ! xxx debug
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    real    :: ftemp_air
    real    :: ftemp_soil
    real    :: dresp_maint = 0.0    ! daily total maintenance respiration


    ! RA: MAINTENANCE RESPIRATION
    ! NPP = GPP - RA
    ! NUP: N-UPTAKE
    ! CEX: C EXUDATION
    ! ADD (NPP - CEX) AND NUP TO LABILE POOL

    ! write(0,*) '---- in npp:'

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft
      if (ispresent(pft,jpngr)) then

        ! xxx consider this: order of PFT within this loop is relevant for N uptake
        ! maybe permute order? maybe divide N availability up for each PFT?
        
        lu=lu_category(pft)
        
        ! reference temperature: 10°C
        ftemp_air  = ftemp( dtemp, "lloyd_and_taylor" )
        ftemp_soil = ftemp( dtemp_soil(lu,jpngr), "lloyd_and_taylor" )

        
        !/////////////////////////////////////////////////////////////////////////
        ! MAINTENANCE RESPIRATION
        ! use function 'resp_main'
        !-------------------------------------------------------------------------
        ! fine roots should have a higher repsiration coefficient than other tissues (Franklin et al., 2007).
        drleaf(pft) = drd(pft)  ! leaf respiration is given by dark respiration as calculated in P-model.       
        drroot(pft) = calc_resp_maint( proot(pft,jpngr)%c%c12 * nind(pft,jpngr), r_root )
        if (tree(pft)) then
          drsapw(pft) = calc_resp_maint( psapw(pft,jpngr)%c%c12 * nind(pft,jpngr), r_sapw )
        endif
        
        ! daily total maintenance repsiration
        dresp_maint = drleaf(pft) + drroot(pft)
        if (tree(pft)) then
          dresp_maint = dresp_maint + drsapw(pft)
        endif

        ! write(0,*) 'dresp_maint',dresp_maint

        !/////////////////////////////////////////////////////////////////////////
        ! DAILY NPP 
        ! NPP is the sum of C available for growth and for N uptake 
        ! This is where isotopic signatures are introduced because only 'dbminc'
        ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
        ! trophic respiration is immediate, it makes thus no sense to calculate 
        ! full isotopic effects of gross exchange _fluxes.
        !-------------------------------------------------------------------------
        drauto(pft) = dresp_maint
        dnpp(pft)   = carbon( dgpp(pft) - drauto(pft) )

        if ( dnpp(pft)%c12 < 0.0 ) then
          write(0,*) 'pft ',pft
          write(0,*) 'drd ',drd(pft)
          write(0,*) 'drroot',drroot(pft)
          write(0,*) 'dgpp ',dgpp(pft)
          write(0,*) 'dnpp ',dnpp(pft)
          write(0,*) 'NPP: dnpp negative'
        end if


        !/////////////////////////////////////////////////////////////////////////
        ! EXUDATION FOR N UPTAKE
        ! This calculates exudation 'dcex', N uptake 'dnup', ...
        ! Labile C exuded for N uptake in interaction with mycorrhiza.
        ! Calculate otpimal C expenditure for N uptake (FUN approach).
        ! PFT loop has to be closed above to get total NPP over all PFTs in each 
        ! LU. 
        !-------------------------------------------------------------------------
        ! SR nuptake calculates dcex and dnup (incl. dnup_act, dnup_pas, ...)
        call nuptake( jpngr, pft )

        ! Add exuded C to exudates pool (fast decay)
        call ccp( carbon( dcex(pft) ), pexud(pft,jpngr) )

        !/////////////////////////////////////////////////////////////////////////
        ! TO LABILE POOL
        ! NPP available for growth first enters the labile pool ('plabl ').
        ! XXX Allocation is called here without "paying"  growth respir.?
        !-------------------------------------------------------------------------
        !print*,'dnpp, dcex, dnup ', dnpp(pft), dcex(pft), dnup(pft)
        !print*,'C adding to plabl', ( dnpp(pft)%c12  )
        !print*,'N adding to plabl', ( dnup(pft) )
        !if (dnup(pft)%n14>0.0) then
        !  print*,'adding to plabl with C:N ratio ', ( dnpp(pft)%c12 - dcex(pft) ) /  dnup(pft)%n14
        !end if
        call orgcp( orgpool( cminus( dnpp(pft), carbon(dcex(pft)) ), dnup(pft) ), plabl(pft,jpngr) )

        ! write(0,*) '---------------in NPP'
        ! write(0,*) 'dclabl ',cminus( dnpp(pft), carbon(dcex(pft)))
        ! write(0,*) 'drd ',drd(pft)
        ! write(0,*) 'drroot',drroot(pft)
        ! write(0,*) 'dgpp(pft) ',dgpp(pft)
        ! write(0,*) 'dnpp(pft) ',dnpp(pft)
        ! write(0,*) 'dcex(pft) ',dcex(pft)
        ! write(0,*) 'dnup(pft) ',dnup(pft)
        ! write(0,*) 'plabl(pft,jpngr) ',plabl(pft,jpngr)

        !-------------------------------------------------------------------------
        ! Leaves are shed in (annual) grasses (=end of vegetation period) when 
        ! labile C pool gets negative.
        !-------------------------------------------------------------------------
        if ( dnpp(pft)%c12 < 0.0 ) then
          if (summergreen(pft)) then
            shedleaves(:,pft)   = .false.
            shedleaves(doy,pft) = .true.
          else
            stop 'labile C negative'
          end if
        end if

      endif
    end do

    ! write(0,*) '---- finished npp'

  end subroutine npp


  function calc_resp_maint( cmass, rresp, ftemp ) result( resp_maint )
    !////////////////////////////////////////////////////////////////
    ! Returns maintenance respiration
    !----------------------------------------------------------------
    ! arguments
    real, intent(in)           :: cmass   ! N mass per unit area [gN/m2]
    real, intent(in)           :: rresp   ! respiration coefficient [gC gC-1 d-1]
    real, intent(in), optional :: ftemp   ! temperature modifier

    ! function return variable
    real :: resp_maint                    ! return value: maintenance respiration [gC/m2]

    if (present(ftemp)) then
      resp_maint = cmass * rresp * ftemp
    else
      resp_maint = cmass * rresp
    end if

  end function calc_resp_maint


  subroutine getpar_modl_npp()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use _params_core, only: ndayyear
    use _sofunutils, only: getreal

    ! Fine-root specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.913 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    call getreal( 'params_npp_lpj.dat', 'r_root', r_root )
    r_root = r_root / ndayyear          ! conversion to rate per day

    ! xxx try:
    r_root = 0.005    

    ! Fine-root specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.044 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    ! (= 0.044 nmol mol-1 s-1; range: 0.5–10, 20 nmol mol-1 s-1 (Landsberg and Sands (2010))
    call getreal( 'params_npp_lpj.dat', 'r_sapw', r_sapw )
    r_sapw = r_sapw / ndayyear  ! conversion to rate per day


  end subroutine getpar_modl_npp

end module _npp
