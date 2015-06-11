module _vegdynamics

  use _params_core

  ! logical, dimension(npft,maxgrid) :: ispresent   ! boolean whether PFT is present
  ! real, dimension(npft,maxgrid)    :: fpc_grid    ! area fraction within gridcell occupied by PFT
  ! real, dimension(npft,maxgrid)    :: nind        ! number of individuals [1/m2]

  ! real, dimension(npft,maxgrid)    :: height      ! tree height (m)
  ! real, dimension(npft,maxgrid)    :: crownarea   ! individual's tree crown area

  logical, dimension(npft)         :: shedleaves  ! true when leaf shedding is triggered

contains

  subroutine estab_daily( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level metabolic N content per unit leaf area as a
    ! function of Vcmax25.
    !------------------------------------------------------------------
    use _params_core, only: npft
    use _params_modl, only: r_ntoc_root
    use _phenology, only: dtphen, summergreen, sprout
    use _vars_core, only: r_ntoc_leaf, pleaf, proot, lai_ind, fpc_ind, &
      sla, fpc_grid, ispresent, crownarea, nind, initpft

    ! xxx debug
    use _classdefs

    ! arguments
    integer, intent(in)               :: jpngr
    integer, intent(in)               :: doy

    ! local variables
    ! logical, save :: firstcall = .true.
    integer :: pft


    do pft=1,npft

      if (summergreen(pft)) then
        !----------------------------------------------------------
        ! GRASSES, summergreen
        !----------------------------------------------------------

        if ( sprout(doy,pft) ) then
          !----------------------------------------------------------
          ! beginning of season
          !----------------------------------------------------------
          ! write(0,*) 'starting to grow on day ',doy
          call initpft( pft, jpngr )

          pleaf(pft,jpngr) = get_sapl( r_ntoc_leaf(pft,jpngr) )
          proot(pft,jpngr) = get_sapl( r_ntoc_root(pft) )
          shedleaves(pft) = .false.

        else if ( dtphen(doy,pft) > 0.0 ) then
          ! write(0,*) 'continued growth on day ',doy

        else
          ! write(0,*) 'no growth on day ',doy

        end if

        ! update 
        lai_ind(pft,jpngr)  = get_lai_ind( pleaf(pft,jpngr)%c%c12, sla(pft,jpngr), crownarea(pft,jpngr) )
        fpc_ind(pft,jpngr)  = get_fapar( lai_ind(pft,jpngr) )
        fpc_grid(pft,jpngr) = get_fpc_grid( crownarea(pft,jpngr), nind(pft,jpngr), fpc_ind(pft,jpngr) )    
 
        ! write(0,*) 'pft   ',pft
        ! write(0,*) 'jpngr   ',jpngr
        ! write(0,*) 'in vegdynamics: pleaf(pft,jpngr) ',pleaf(pft,jpngr)
        ! write(0,*) 'in vegdynamics: proot(pft,jpngr) ',proot(pft,jpngr)
        ! write(0,*) 'C:N of pleaf     ',cton(pleaf(pft,jpngr), 0.0 )
        ! write(0,*) 'C:N of proot     ',cton(proot(pft,jpngr), 0.0 )

        ! write(0,*) 'in vegdynamics: lai_ind(pft,jpngr) ',lai_ind(pft,jpngr)
        ! write(0,*) 'in vegdynamics: sla(pft,jpngr)   ',sla(pft,jpngr)
        ! write(0,*) 'in vegdynamics: crownarea(pft,jpngr)   ',crownarea(pft,jpngr)
        ! write(0,*) 'in vegdynamics: fpc_ind(pft,jpngr) ',fpc_ind(pft,jpngr)
        ! write(0,*) 'fpc_grid(pft,jpngr)',fpc_grid(pft,jpngr)
! 
      else

        stop 'estab_daily not implemented for trees'

      end if

    end do

  end subroutine estab_daily


  function get_sapl( r_ntoc ) result( sapl )
    !//////////////////////////////////////////////////////////////////
    ! Calculates the sapling C and N content in respective pool
    ! (variable of type 'orgpool'), given C:N ratio in that pool. 
    !------------------------------------------------------------------
    use _classdefs

    ! arguments
    real, intent(in) :: r_ntoc

    ! local variables
    real, parameter :: sapling_c_init = 0.1     ! xxx to parameters?

    ! function return variable
    type(orgpool) :: sapl    ! sapling pool mass

    sapl%c%c12 = sapling_c_init
    sapl%n%n14 = sapling_c_init * r_ntoc

  end function get_sapl


  function get_fpc_grid( crownarea, nind, fpc_ind ) result( fpc_grid )
    !////////////////////////////////////////////////////////////////
    ! Function returns total fractional plant cover of a PFT
    ! Eq. 8 in Sitch et al., 2003
    !----------------------------------------------------------------

    ! arguments
    real, intent(in) :: crownarea
    real, intent(in) :: nind
    real, intent(in) :: fpc_ind

    ! function return variable
    real, intent(out) :: fpc_grid

    fpc_grid = crownarea * nind * fpc_ind

  end function get_fpc_grid


  function get_fapar( lai ) result( fapar )
    !////////////////////////////////////////////////////////////////
    ! Function returns fractional plant cover an individual
    ! Eq. 7 in Sitch et al., 2003
    !----------------------------------------------------------------
    use _params_modl, only: kbeer

    ! arguments
    real, intent(in) :: lai

    ! function return variable
    real, intent(out) :: fapar

    fapar = ( 1.0 - exp( -1.0 * kbeer * lai) )

  end function get_fapar


  function get_lai_ind( cleaf, sla, crownarea ) result( lai_ind )
    !////////////////////////////////////////////////////////////////
    ! Function returns leaf area index of an individual
    ! Eq. 5 in Sitch et al., 2003
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: cleaf
    real, intent(in) :: sla
    real, intent(in) :: crownarea

    ! function return variable
    real, intent(out) :: lai_ind

#if _check_sanity
    if (crownarea<0.0) stop 'GET_LAI_IND: crownarea < 0'
#endif

    ! write(0,*) '----'
    ! write(0,*) 'in get_lai_ind:'
    ! write(0,*) 'cleaf ',cleaf
    ! write(0,*) 'sla ',sla
    ! write(0,*) 'crownarea ',crownarea
    ! write(0,*) '----'

    if (crownarea==0.0) then
      lai_ind = 0.0
    else
      lai_ind = cleaf * sla / crownarea
    endif

  end function get_lai_ind

end module _vegdynamics
