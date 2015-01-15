subroutine establishment( jpngr )
  !//////////////////////////////////////////////////////////////////////////////
  ! ESTABLISHMENT SUBROUTINE
  ! Defined here for fixed and dynamic vegetation. Activate either using
  ! sofun_module_control.inc
  ! b.stocker@imperial.ac.uk
  !-----------------------------------------------------------------------------
#include "sofun_module_control.inc"
  use classdefs
  use params_site, only: lpft1, lpft2, lpft3, lpft4, lpft5, lpft6,&
    lpft7, lpft8, lpft9, in_fpc_grid_1, in_fpc_grid_2, in_fpc_grid_3, &
    in_fpc_grid_4, in_fpc_grid_5, in_fpc_grid_6, in_fpc_grid_7, &
    in_fpc_grid_8, in_fpc_grid_9
  use treegeometry
  use pools
  use params_modl, only: tree

#ifndef _fix_veg
  use params_core, only: nlu
  use fluxes, only: aestab
  use gridvars, only: lu_area
#endif

  implicit none

  ! ARGUMENTS
  integer, intent(in)      :: jpngr

  ! LOCAL VARIABLES
  integer                  :: pft
  logical, dimension(npft) :: estab
  type(orgpool)            :: psapw_temp

#ifndef _fix_veg
  integer                  :: pft
  integer                  :: lu
  real, dimension(nlu)     :: fpc_total
  real, dimension(nlu)     :: fpc_tree_total   !total grid FPC for tree PFTs
  real, dimension(nlu)     :: fpc_grass_total  !total grid FPC for tree PFTs
  integer                  :: ntrees
  integer                  :: ngrass
  real                     :: bare
  real                     :: estab_rate       ! sapling establishment rate over area available for establishment (indiv/m2)
  real, dimension(nlu)     :: nind_add         ! grid-level establishment rate (indiv/m2) (=estab_grid)
  real                     :: nind_old         ! previous (year's) nind
#endif


#ifdef _fix_veg
  !//////////////////////////////////////////////////////////////////////////////
  ! GET FIXED PRESCRIBED VEGETATION COVER
  ! Subroutine defines 'fpc_grid' and updates 'nind' following geo-
  ! metry constraints.
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! define which PFTs to establish
  !-----------------------------------------------------------------------------
  estab(:) = .false.
  if (lpft1) then
    estab(1)           = .true.
    ispresent(1,jpngr) = .true.
  end if
  if (lpft2) then
    estab(2)           = .true.
    ispresent(2,jpngr) = .true.
  end if
  if (lpft3) then
    estab(3)           = .true.
    ispresent(3,jpngr) = .true.
  end if
  if (lpft4) then
    estab(4)           = .true.
    ispresent(4,jpngr) = .true.
  end if
  if (lpft5) then
    estab(5)           = .true.
    ispresent(5,jpngr) = .true.
  end if
  if (lpft6) then
    estab(6)           = .true.
    ispresent(6,jpngr) = .true.
  end if
  if (lpft7) then
    estab(7)           = .true.
    ispresent(7,jpngr) = .true.
  end if
  if (lpft8) then
    estab(8)           = .true.
    ispresent(8,jpngr) = .true.
  end if
  if (lpft9) then
    estab(9)           = .true.
    ispresent(9,jpngr) = .true.
  end if

  !-----------------------------------------------------------------------------
  ! Get fractional plant cover from site input file
  !-----------------------------------------------------------------------------
  if (lpft1) fpc_grid(1,jpngr) = in_fpc_grid_1
  if (lpft2) fpc_grid(2,jpngr) = in_fpc_grid_2
  if (lpft3) fpc_grid(3,jpngr) = in_fpc_grid_3
  if (lpft4) fpc_grid(4,jpngr) = in_fpc_grid_4
  if (lpft5) fpc_grid(5,jpngr) = in_fpc_grid_5
  if (lpft6) fpc_grid(6,jpngr) = in_fpc_grid_6
  if (lpft7) fpc_grid(7,jpngr) = in_fpc_grid_7
  if (lpft8) fpc_grid(8,jpngr) = in_fpc_grid_8
  if (lpft9) fpc_grid(9,jpngr) = in_fpc_grid_9

  ! xxx debug
  !print*,'A fpc_grid  ',fpc_grid(:,1)
  !print*,'A crownarea ',crownarea(:,1)
  !print*,'A leaf mass ',pleaf(:,1)
  !print*,'A labl      ',plabl(:,1)

  do pft=1,npft

    !-----------------------------------------------------------------------------
    ! Re-calculate crownarea
    !-----------------------------------------------------------------------------
    psapw_temp = psapw(pft,jpngr)
    call update_tree_geometry( &
      pft, jpngr, psapw_temp, ntoc( psapw(pft,jpngr), default=0.0 ) &
      )

    !-----------------------------------------------------------------------------
    ! Update fpc_grid, fpc_ind, and lai_ind
    !-----------------------------------------------------------------------------
    call update_fpc( pft, jpngr )

    !-----------------------------------------------------------------------------
    ! Re-calculate 'nind' based on new leaf mass and crownarea
    !-----------------------------------------------------------------------------
    nind(pft,jpngr) = get_nind( tree(pft), fpc_grid(pft,jpngr), crownarea(pft,jpngr), fpc_ind(pft,jpngr) )

  end do

  ! xxx debug
  !print*,'B fpc_grid ',fpc_grid(:,1)
  !print*,'B crownarea ',crownarea(:,1)


  return


#else
  !////////////////////////////////////////////////////////////////
  ! DEFINE DYNAMIC VEGETATION COVER
  ! Subroutine adds sapling mass and updates 'nind'. Tree geometry
  ! is updated using updated Cleaf, Croot, and (Cwood+Csap).
  ! b.stocker@imperial.ac.uk
  !----------------------------------------------------------------

  ! was done here: Check if 'survive' is false or nind<nind_min. If so, kill PFT.

  ! was done here: if aprec(jpngr)<aprec_min_estab; then estab=.false.
  ! minimum precipitation for establishment to bioclim? (100 mm/yr used in LPX)

  ! xxx debug
  !print*,'doing establishment for dynamic vegetation'

  !-----------------------------------------------------------------------------
  ! DYNAMIC ESTABLISHMENT ON BARE GROUND
  !-----------------------------------------------------------------------------
  ! define which PFTs to establish
  estab(:) = .false.
  if (lpft1) estab(1) = .true.
  if (lpft2) estab(2) = .true.
  if (lpft3) estab(3) = .true.
  if (lpft4) estab(4) = .true.
  if (lpft5) estab(5) = .true.
  if (lpft6) estab(6) = .true.
  if (lpft7) estab(7) = .true.
  if (lpft8) estab(8) = .true.
  if (lpft9) estab(9) = .true.

  ! START LU-LOOP
  do lu=1,nlu
  
    !-----------------------------------------------------------------------------
    ! Get total area of trees and grassses (FPC = fractional plant cover)
    !-----------------------------------------------------------------------------
    fpc_tree_total(lu)  = sum( fpc_grid(:,jpngr), mask=tree(:).and.islu(:,lu) )
    fpc_grass_total(lu) = sum( fpc_grid(:,jpngr), mask=grass(:).and.islu(:,lu) )
    fpc_total(lu)       = sum( fpc_grid(:,jpngr), mask=islu(:,lu) )

    !------------------------------------------------------------------------------
    ! Set 'estab' to false if FPC exceeds 1
    ! xxx this allows FPC>1 when trees and grasses coexist.
    !------------------------------------------------------------------------------
    do pft=1,npft
      if (islu(pft,lu)) then
        if(lu_area(lu,jpngr)<1.0e-12) estab(pft)=.false.
        if(tree(pft)) then
          if (fpc_tree_total(lu)>=1.0) estab(pft)=.false.
        elseif (grass(pft)) then 
          if (fpc_grass_total(lu)>=1.0) estab(pft)=.false.
        end if 
      end if
    end do


    !------------------------------------------------------------------------------
    ! introduce new PFTs if any
    !------------------------------------------------------------------------------
    do pft=1,npft
      if (.not.ispresent(pft,jpngr).and.estab(pft)) then !note aprec condition is included in estab
        ispresent(pft,jpngr) = .true.
        ! initialise pools with 0. initialisation with sapling pool size is done below.
        call initpft(pft,jpngr)
      end if
    end do


    !------------------------------------------------------------------------------
    ! count number of establishing PFTs
    !------------------------------------------------------------------------------
    ntrees=0
    ngrass=0
    do pft=1,npft
      if (islu(pft,lu)) then
        if(ispresent(pft,jpngr)) then
          if(estab(pft)) then
            if(tree(pft)) then  !trees
              ntrees=ntrees+1 
            elseif (grass(pft)) then !grasses
              ngrass=ngrass+1         
            end if
          end if
        end if
      end if
    end do
    
    !------------------------------------------------------------------------------
    ! calculate establishment rate, for grasses ~ bare ground
    !------------------------------------------------------------------------------

    ! Calculate grid-level establishment rate per woody PFT
    ! Space available for woody PFT establishment is proportion of grid cell
    ! not currently occupied by woody PFTs.
    ! xxx parameters 0.24 (=estab_max) and 5.0 (no explicit name in LPX) are
    ! adopted from LPX. Maybe add to modelparameters? Really necessary?
    if (ntrees>0.0) then
      estab_rate = max(0.0,0.24*(1.0-exp(5.0* &
          (fpc_tree_total(lu)-1.0))))/real(ntrees)
      nind_add(lu) = estab_rate*(1.0-fpc_tree_total(lu))
    else
      nind_add(lu) = 0.0
    end if

    ! Grasses and mosses can establish in non-vegetated areas
    if (ngrass>0.0) then
      bare = (1.0-fpc_total(lu))/(real(ngrass))
    end if



    do pft=1,npft
      
      if (islu(pft,lu)) then

        if (ispresent(pft,jpngr)) then

          if (estab(pft)) then
            ! ----------------------------------------------------------------------                        
            ! initialise with sapling pool size
            ! ----------------------------------------------------------------------                        

            if (tree(pft)) then
              ! ----------------------------------------------------------------------                        
              ! TREE ESTABLISHMENT
              ! ----------------------------------------------------------------------                        
              ! Add new saplings to current population
              nind_old=nind(pft,jpngr)
              nind(pft,jpngr)=nind_old+nind_add(lu)
              

              pleaf(pft,jpngr) = distributed_to_nind( &
                  pleaf(pft,jpngr), nind_old          &
                , lm_sapl(pft),     nind_add(lu)      &
                , nind(pft,jpngr)                     &
                )

              pwood(pft,jpngr) = distributed_to_nind( &
                  pwood(pft,jpngr), nind_old          &
                , hm_sapl(pft),     nind_add(lu)      &
                , nind(pft,jpngr)                     &
                )

              proot(pft,jpngr) = distributed_to_nind( &
                  proot(pft,jpngr), nind_old          &
                , rm_sapl(pft),     nind_add(lu)      &
                , nind(pft,jpngr)                     &
                )

              psapw_temp = distributed_to_nind( &
                  psapw(pft,jpngr), nind_old          &
                , sm_sapl(pft),     nind_add(lu)      &
                , nind(pft,jpngr)                     &
                )

              ! record implicit C and N fixation
              call orgcp( orgfrac(nind_add(lu), lm_sapl(pft)), aestab(lu) )
              call orgcp( orgfrac(nind_add(lu), sm_sapl(pft)), aestab(lu) )
              call orgcp( orgfrac(nind_add(lu), hm_sapl(pft)), aestab(lu) )
              call orgcp( orgfrac(nind_add(lu), rm_sapl(pft)), aestab(lu) )

             
              ! Calculate height, diameter and crown area for new average
              ! individual such that the basic allometric relationships (A-C below)
              ! are satisfied.
              call update_tree_geometry( pft, jpngr, psapw_temp, ntoc(psapw(pft,jpngr)) )

            else if (grass(pft)) then
              !------------------------------------------------------------------------------
              ! GRASS ESTABLISHMENT
              !------------------------------------------------------------------------------
              call orgcp( orgfrac(bare,lm_sapl(pft)), pleaf(pft,jpngr) )            
              call orgcp( orgfrac(bare,rm_sapl(pft)), proot(pft,jpngr) )

              ! record implicit C and N fixation
              call orgcp( orgfrac(bare,lm_sapl(pft)), aestab(lu) )        
              call orgcp( orgfrac(bare,rm_sapl(pft)), aestab(lu) )        

              ! if PFT is not yet ispresent:
              ! add C13 signature for saplings: C4: co2(2)-3.6,    C3: co2(2)-17.8
              ! add C14 signature for saplings: C4: co2(2)-3.6*2., C3: co2(2)-17.8*2.

              ! if PFT is ispresent: 
              ! add C13 and C14 signature from annual NPP

            end if

          end if  ! estab

          !------------------------------------------------------------------------------
          ! update LAI and FPC
          !------------------------------------------------------------------------------
          call update_fpc( pft, jpngr )

        end if  ! ispresent

      end if  ! islu

    end do                     !pft

  ! END LU LOOP
  end do

  return

contains

  function distributed_to_nind( pool1, nind1, pool2, nind2, nind_new )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable of type 'orgpool' and size
    !  of a fraction 'frac' of source pool ('from')
    !----------------------------------------------------------------
    implicit none

   ! arguments
    real, intent(in)           :: nind1, nind2, nind_new
    type (orgpool), intent(in) :: pool1, pool2
    
    ! function return variable
    type (orgpool) :: distributed_to_nind

    distributed_to_nind%c%c12 = (pool1%c%c12*nind1 + pool2%c%c12*nind2) &
      / nind_new

    distributed_to_nind%n%n14 = (pool1%n%n14*nind1 + pool2%n%n14*nind2) &
      / nind_new

  end function distributed_to_nind

#endif

end subroutine establishment

