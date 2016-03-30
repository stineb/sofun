module md_npp
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
  use md_classdefs
  use md_plant

  implicit none

  private
  public npp

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
    use md_params_core, only: npft, ndayyear
    use md_soiltemp, only: dtemp_soil
    use md_gpp, only: dgpp, drd
    use md_phenology, only: shedleaves

    ! arguments
    integer, intent(in) :: jpngr
    real, intent(in)    :: dtemp      ! air temperature at this day

    ! xxx debug
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu

    ! print*, '---- in npp:'

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft

      if ( ispresent(pft,jpngr) ) then

        ! print*, '---------------in NPP'
        ! print*, 'drd ',drd(pft)
        ! print*, 'dgpp(pft) ',dgpp(pft)
        ! print*, 'dnpp(pft) ',dnpp(pft)
        ! print*, 'dcex(pft) ',dcex(pft)
        ! print*, 'plabl(pft,jpngr) ',plabl(pft,jpngr)

        lu = params_pft_plant(pft)%lu_category
        
        !/////////////////////////////////////////////////////////////////////////
        ! MAINTENANCE RESPIRATION
        ! use function 'resp_main'
        !-------------------------------------------------------------------------
        ! fine roots should have a higher repsiration coefficient than other tissues (Franklin et al., 2007).
        drleaf(pft) = drd(pft)  ! leaf respiration is given by dark respiration as calculated in P-model.       
        drroot(pft) = calc_resp_maint( proot(pft,jpngr)%c%c12 * nind(pft,jpngr), params_plant%r_root, dtemp )
        if (params_pft_plant(pft)%tree) then
          drsapw(pft) = calc_resp_maint( psapw(pft,jpngr)%c%c12 * nind(pft,jpngr), params_plant%r_sapw, dtemp )
        endif
                
        !/////////////////////////////////////////////////////////////////////////
        ! DAILY NPP 
        ! NPP is the sum of C available for growth and for N uptake 
        ! This is where isotopic signatures are introduced because only 'dbminc'
        ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
        ! trophic respiration is immediate, it makes thus no sense to calculate 
        ! full isotopic effects of gross exchange _fluxes.
        !-------------------------------------------------------------------------
        dnpp(pft)   = carbon( dgpp(pft) - drgrow(pft) - drleaf(pft) - drroot(pft) )

        if ( dnpp(pft)%c12 < 0.0 ) then
          print*, 'pft    ',pft
          print*, 'drleaf ',drleaf(pft)
          print*, 'drroot ',drroot(pft)
          print*, 'dgpp   ',dgpp(pft)
          print*, 'dnpp   ',dnpp(pft)
          print*, 'NPP: dnpp negative'
          stop
        end if

        !/////////////////////////////////////////////////////////////////////////
        ! EXUDATION FOR N UPTAKE
        ! This calculates exudation 'dcex', N uptake 'dnup', ...
        ! Labile C exuded for N uptake in interaction with mycorrhiza.
        ! Calculate otpimal C expenditure for N uptake (FUN approach).
        ! PFT loop has to be closed above to get total NPP over all PFTs in each 
        ! LU. 
        !-------------------------------------------------------------------------
        ! dnup(pft) = nitrogen(0.0) ! XXX WILL BE DETERMINED IN ALLOCATION
        dcex(pft) = calc_cexu( proot(pft,jpngr)%c%c12 )      

        ! ! SR nuptake calculates dcex and dnup (incl. dnup_act, dnup_pas, ...)
        ! call nuptake( jpngr, pft )

        ! Add exuded C to exudates pool (fast decay)
        call ccp( carbon( dcex(pft) ), pexud(pft,jpngr) )

        !/////////////////////////////////////////////////////////////////////////
        ! TO LABILE POOL
        ! NPP available for growth first enters the labile pool ('plabl ').
        ! XXX Allocation is called here without "paying"  growth respir.?
        !-------------------------------------------------------------------------
        call orgcp( orgpool( cminus( dnpp(pft), carbon(dcex(pft)) ), dnup(pft) ), plabl(pft,jpngr) )

        ! print*, '---------------in NPP'
        ! print*, 'plabl  ', plabl(pft,jpngr)
        ! print*, 'drd    ', drd(pft)
        ! print*, 'drroot ',drroot(pft)
        ! print*, 'dgpp(pft) ',dgpp(pft)
        ! print*, 'dnpp(pft) ',dnpp(pft)
        ! print*, 'dnup(pft) ',dnup(pft)
        ! print*, 'dcex(pft) ',dcex(pft)
        ! if (doy==39) stop

        ! !-------------------------------------------------------------------------
        ! ! Leaves are shed in (annual) grasses (=end of vegetation period) when 
        ! ! labile C pool gets negative.
        ! !-------------------------------------------------------------------------
        ! if ( dnpp(pft)%c12 < 0.0 ) then
        !   if (summergreen(pft)) then
        !     shedleaves(:,pft)   = .false.
        !     shedleaves(doy,pft) = .true.
        !   else
        !     stop 'labile C negative'
        !   end if
        ! end if

      else

        dnpp(pft)   = carbon(0.0)
        drleaf(pft) = 0.0
        drroot(pft) = 0.0
        drsapw(pft) = 0.0
        dcex(pft)   = 0.0

      endif
    end do

    ! print*, '---- finished npp'

  end subroutine npp


  function calc_resp_maint( cmass, rresp, dtemp ) result( resp_maint )
    !////////////////////////////////////////////////////////////////
    ! Returns maintenance respiration
    !----------------------------------------------------------------
    use md_sofunutils, only: ftemp
    use md_gpp, only: ramp_gpp_lotemp     ! same ramp as for GPP 

    ! arguments
    real, intent(in)           :: cmass   ! N mass per unit area [gN/m2]
    real, intent(in)           :: rresp   ! respiration coefficient [gC gC-1 d-1]
    real, intent(in), optional :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: resp_maint                    ! return value: maintenance respiration [gC/m2]

    if (present(dtemp)) then
      resp_maint = cmass * rresp * ftemp( dtemp, "lloyd_and_taylor" ) * ramp_gpp_lotemp( dtemp )
    else
      resp_maint = cmass * rresp * ramp_gpp_lotemp( dtemp )
    end if

  end function calc_resp_maint


  function calc_cexu( croot ) result( cexu )
    !/////////////////////////////////////////////////////////////////
    ! Constant exudation rate
    !-----------------------------------------------------------------
    ! arguments
    real, intent(in)  :: croot

    ! function return variable
    real, intent(out) :: cexu

    cexu = params_plant%exurate * croot

  end function calc_cexu

end module md_npp
