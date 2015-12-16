module _vegdynamics
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  use _params_core

  ! logical, dimension(npft,maxgrid) :: ispresent   ! boolean whether PFT is present
  ! real, dimension(npft,maxgrid)    :: fpc_grid    ! area fraction within gridcell occupied by PFT
  ! real, dimension(npft,maxgrid)    :: nind        ! number of individuals [1/m2]

  ! real, dimension(npft,maxgrid)    :: height      ! tree height (m)
  ! real, dimension(npft,maxgrid)    :: crownarea   ! individual's tree crown area

  implicit none

contains

  subroutine estab_daily( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level metabolic N content per unit leaf area as a
    ! function of Vcmax25.
    !------------------------------------------------------------------
    use _params_core, only: npft
    use _phenology, only: dtphen, summergreen, sprout
    use _vars_core, only: pleaf, proot, lai_ind, fpc_ind, &
      fpc_grid, ispresent, crownarea, nind, initpft

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
          call add_sapl( pft, jpngr )

        end if

        !----------------------------------------------------------
        ! Update fpc_grid (and fpc_ind)
        !----------------------------------------------------------
        call update_fpc_grid( pft, jpngr )

        ! write(0,*) 'pft   ',pft
        ! write(0,*) 'jpngr   ',jpngr
        ! write(0,*) 'in vegdynamics: pleaf(pft,jpngr) ',pleaf(pft,jpngr)
        ! write(0,*) 'in vegdynamics: proot(pft,jpngr) ',proot(pft,jpngr)
        ! write(0,*) 'C:N of pleaf     ',cton(pleaf(pft,jpngr), 0.0 )
        ! write(0,*) 'C:N of proot     ',cton(proot(pft,jpngr), 0.0 )

        ! write(0,*) 'in vegdynamics: lai_ind(pft,jpngr) ',lai_ind(pft,jpngr)
        ! write(0,*) 'in vegdynamics: crownarea(pft,jpngr)   ',crownarea(pft,jpngr)
        ! write(0,*) 'in vegdynamics: fpc_ind(pft,jpngr) ',fpc_ind(pft,jpngr)
        ! write(0,*) 'fpc_grid(pft,jpngr)',fpc_grid(pft,jpngr)

      else

        stop 'estab_daily not implemented for trees'

      end if

    end do

  end subroutine estab_daily


  subroutine update_foliage_vars( pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Updates PFT-specific state variables after change in LAI.
    !------------------------------------------------------------------
    use _vars_core, only: pleaf, lma, crownarea, lai_ind, nind, fpc_ind, fpc_grid

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    if (pleaf(pft,jpngr)%c%c12==0.0) then
      lai_ind(pft,jpngr)  = 0.0
      fpc_ind(pft,jpngr)  = 0.0
      fpc_grid(pft,jpngr) = 0.0
    else
      ! This assumes that leaf canopy-average traits (LMA) do not change upon changes in LAI.
      lai_ind(pft,jpngr) = pleaf(pft,jpngr)%c%c12 / ( lma(pft,jpngr) * crownarea(pft,jpngr) * nind(pft,jpngr) )
      call update_fpc_grid( pft, jpngr )
    end if

  end subroutine update_foliage_vars


  subroutine update_fpc_grid( pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Updates PFT-specific state variables after change in LAI.
    !------------------------------------------------------------------
    use _vars_core, only: lai_ind, fpc_ind, fpc_grid, crownarea, nind

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    fpc_ind(pft,jpngr)  = get_fapar( lai_ind(pft,jpngr) )
    fpc_grid(pft,jpngr) = get_fpc_grid( crownarea(pft,jpngr), nind(pft,jpngr), fpc_ind(pft,jpngr) ) 

  end subroutine update_fpc_grid


  subroutine add_sapl( pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! 
    !------------------------------------------------------------------
    use _classdefs
    use _vars_core, only: plabl

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    real, parameter :: clabl_sapl = 5.0
    real, parameter :: nlabl_sapl = 0.12
    ! real, parameter :: nlabl_sapl = 0.5
    ! real, parameter :: cleaf_sapl = 5.0
    ! real, parameter :: croot_sapl = 1.0

    plabl(pft,jpngr)%c%c12 = clabl_sapl
    plabl(pft,jpngr)%n%n14 = nlabl_sapl

    ! pleaf(pft,jpngr)%c%c12 = cleaf_sapl
    ! proot(pft,jpngr)%c%c12 = croot_sapl

    ! pleaf(pft,jpngr)%n%n14 = cleaf_sapl * r_ntoc_leaf(pft,jpngr)
    ! proot(pft,jpngr)%n%n14 = croot_sapl * r_ntoc_root(pft)

  end subroutine add_sapl


  function get_fpc_grid( crownarea, nind, fpc_ind ) result( fpc_grid )
    !////////////////////////////////////////////////////////////////
    ! FRACTIONAL PLANT COVERAGE
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
    ! FOLIAGE PROJECTIVE COVER 
    ! = Fraction of Absorbed Photosynthetically Active Radiation
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


!   function get_lai_ind( cleaf, sla, crownarea ) result( lai_ind )
!     !////////////////////////////////////////////////////////////////
!     ! Function returns leaf area index of an individual
!     ! Eq. 5 in Sitch et al., 2003
!     !----------------------------------------------------------------
!     ! arguments
!     real, intent(in) :: cleaf
!     real, intent(in) :: sla
!     real, intent(in) :: crownarea

!     ! function return variable
!     real, intent(out) :: lai_ind

! #if _check_sanity
!     if (crownarea<0.0) stop 'GET_LAI_IND: crownarea < 0'
! #endif

!     ! write(0,*) '----'
!     ! write(0,*) 'in get_lai_ind:'
!     ! write(0,*) 'cleaf ',cleaf
!     ! write(0,*) 'sla ',sla
!     ! write(0,*) 'crownarea ',crownarea
!     ! write(0,*) '----'

!     if (crownarea==0.0) then
!       lai_ind = 0.0
!     else
!       lai_ind = cleaf * sla / crownarea
!     endif

!   end function get_lai_ind

end module _vegdynamics
