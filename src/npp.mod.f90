module md_npp
  !////////////////////////////////////////////////////////////////
  ! NPP_LPJ MODULE
  ! Contains the "main" subroutine 'npp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'npp' must contain this list 
  ! of subroutines (names that way).
  !   - npp
  !   - ((interface%steering%init))io_npp
  !   - ((interface%steering%init))output_npp
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
  public npp, calc_cexu, calc_resp_maint, deactivate_root

contains

  subroutine npp( jpngr, dtemp )
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

    ! arguments
    integer, intent(in) :: jpngr
    real, intent(in)    :: dtemp      ! air temperature at this day

    ! local variables
    integer :: pft
    integer :: lu
    real :: avl

    ! print*, '---- in npp:'

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft

      if ( ispresent(pft,jpngr) ) then

        if (plabl(pft,jpngr)%c%c12<0.0) stop 'before npp labile C is neg.'
        if (plabl(pft,jpngr)%n%n14<0.0) stop 'before npp labile N is neg.'

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
        ! DAILY NPP AND C EXPORT
        ! NPP is the sum of C available for growth and for N uptake 
        ! This is where isotopic signatures are introduced because only 'dbminc'
        ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
        ! trophic respiration is immediate, it makes thus no sense to calculate 
        ! full isotopic effects of gross exchange _fluxes.
        ! Growth respiration ('drgrow') is deduced from 'dnpp' in allocation SR.
        !-------------------------------------------------------------------------
        dnpp(pft) = carbon( dgpp(pft) - drleaf(pft) - drroot(pft) )
        dcex(pft) = calc_cexu( proot(pft,jpngr)%c%c12 , dtemp )     
        avl       = plabl(pft,jpngr)%c%c12 + dnpp(pft)%c12 - dcex(pft)

        ! If C used for root respiration and export is not available, then reduce 
        ! root mass to match 
        if ( avl < 0.0 ) then
          call deactivate_root( dgpp(pft), drleaf(pft), plabl(pft,jpngr)%c%c12, proot(pft,jpngr), drroot(pft), dnpp(pft)%c12, dcex(pft), dtemp, plitt_bg(pft,jpngr) )
        end if

        !/////////////////////////////////////////////////////////////////////////
        ! TO LABILE POOL
        ! NPP available for growth first enters the labile pool ('plabl ').
        ! XXX Allocation is called here without "paying"  growth respir.?
        !-------------------------------------------------------------------------
        call ccp( carbon( dcex(pft) ), pexud(pft,jpngr) )
        call ccp( cminus( dnpp(pft), carbon(dcex(pft)) ), plabl(pft,jpngr)%c )

        if (plabl(pft,jpngr)%c%c12< -1.0e-13) stop 'after npp labile C is neg.'
        if (plabl(pft,jpngr)%n%n14< -1.0e-13) stop 'after npp labile N is neg.'

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


  subroutine deactivate_root( mygpp, mydrleaf, myplabl, myproot, rroot, npp, cexu, dtemp, myplitt )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates amount of root mass supportable by (GPP-Rd+Clabl='avl'), so that
    ! NPP is zero and doesn't get negative. Moves excess from pool 'myproot' to
    ! pool 'myplitt'.
    !-------------------------------------------------------------------------
    ! argument
    real, intent(in) :: mygpp
    real, intent(in) :: mydrleaf
    real, intent(in) :: myplabl
    type( orgpool ), intent(inout) :: myproot
    real, intent(out) :: rroot
    real, intent(out) :: npp
    real, intent(out) :: cexu
    real, intent(in) :: dtemp
    type( orgpool ), intent(inout), optional :: myplitt
    
    ! local variables
    real :: croot_trgt
    real :: droot
    type( orgpool ) :: rm_turn

    real, parameter :: safety = 0.9999

    ! calculate target root mass
    croot_trgt = safety * ( mygpp - mydrleaf + myplabl) / ( params_plant%r_root + params_plant%exurate )
    droot      = ( 1.0 - croot_trgt / myproot%c%c12 )
    rm_turn    = orgfrac( droot, myproot )
    if (present(myplitt)) then
      call orgmv( rm_turn, myproot, myplitt )
    else
      myproot = orgminus( myproot, rm_turn )
    end if

    ! update fluxes based on corrected root mass
    rroot = calc_resp_maint( myproot%c%c12, params_plant%r_root, dtemp )
    npp   = mygpp - mydrleaf - rroot
    cexu  = calc_cexu( myproot%c%c12 , dtemp )     

  end subroutine deactivate_root


  function calc_resp_maint( cmass, rresp, dtemp ) result( resp_maint )
    !////////////////////////////////////////////////////////////////
    ! Returns maintenance respiration
    !----------------------------------------------------------------
    use md_rates, only: ftemp
    use md_gpp, only: ramp_gpp_lotemp     ! same ramp as for GPP 

    ! arguments
    real, intent(in)           :: cmass   ! N mass per unit area [gN/m2]
    real, intent(in)           :: rresp   ! respiration coefficient [gC gC-1 d-1]
    real, intent(in), optional :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: resp_maint                    ! return value: maintenance respiration [gC/m2]

    resp_maint = cmass * rresp * ramp_gpp_lotemp( dtemp )

    ! LPX-like temperature dependeneo of respiration rates
    ! resp_maint = cmass * rresp * ftemp( dtemp, "lloyd_and_taylor" ) * ramp_gpp_lotemp( dtemp )

  end function calc_resp_maint


  function calc_cexu( croot, dtemp ) result( cexu )
    !/////////////////////////////////////////////////////////////////
    ! Constant exudation rate
    !-----------------------------------------------------------------
    use md_gpp, only: ramp_gpp_lotemp     ! same ramp as for GPP 

    ! arguments
    real, intent(in)           :: croot
    real, intent(in), optional :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: cexu

    cexu = params_plant%exurate * croot * ramp_gpp_lotemp( dtemp )

  end function calc_cexu

end module md_npp
