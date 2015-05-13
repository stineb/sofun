module _dynveg


contains

  subroutine establishment( jpngr )
    !//////////////////////////////////////////////////////////////////////////////
    ! ESTABLISHMENT SUBROUTINE
    ! Defined here for fixed and dynamic vegetation. Activate either using
    ! sofun_module_control.inc
    ! b.stocker@imperial.ac.uk
    !-----------------------------------------------------------------------------
    use _classdefs
    use _params_site, only: lpft1, lpft2, lpft3, lpft4, lpft5, lpft6,&
      lpft7, lpft8, lpft9, in_fpc_grid_1, in_fpc_grid_2, in_fpc_grid_3, &
      in_fpc_grid_4, in_fpc_grid_5, in_fpc_grid_6, in_fpc_grid_7, &
      in_fpc_grid_8, in_fpc_grid_9
    use _vegdynamics
    use _pools
    use _params_modl, only: tree

    use _params_core, only: nlu
    use _land, only: ndlu,aestab
    use gridvars, only: lu_area

    implicit none

    ! ARGUMENTS
    integer, intent(in)      :: jpngr

    ! LOCAL VARIABLES
    integer                  :: pft
    logical, dimension(npft) :: estab
    type(orgpool)            :: psapw_temp

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
          ! initialise _pools with 0. initialisation with sapling pool size is done below.
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

  end subroutine establishment


  subroutine light( jpngr )
    !////////////////////////////////////////////////////////////////
    ! LIGHT COMPETITION
    ! reduce FPC of trees to 'fpc_tree_max', killing according number of
    ! individuals (=reducing 'nind')
    ! only effective with dynamic vegetation
    ! b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use _classdefs
    use _params_core, only: nlu, pft
    use _params_modl, only: sla, islu, tree, grass, fpc_tree_max, kbeer
    use _vegdynamics, only: fpc_grid, ispresent, nind, update_fpc
    use _pools
    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr


    ! LOCAL VARIABLES
    integer :: pft,lu

    real :: fpc_tree_total  ! total grid FPC for tree PFTs
    real :: fpc_grass_total ! total grid FPC for grass PFTs
    real :: fpc_total       ! total grid FPC for all PFTs
    real :: excess          ! tree FPC or grass cover to be reduced
    real :: nind_kill       ! reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
    real :: cleaf_target    ! required grass leaf mass so that FPC constraints are satisfied
    real :: frac_kill       ! fraction of grass PFT _pools removed so that FPC constraints are satisfied


    do lu=1,nlu

      !-----------------------------------------------------------------------------
      ! Get total area of trees and grassses (FPC = fractional plant cover)
      !-----------------------------------------------------------------------------
      fpc_tree_total  = sum( fpc_grid(:,jpngr), mask=tree(:).and.islu(:,lu) )
      fpc_grass_total = sum( fpc_grid(:,jpngr), mask=grass(:).and.islu(:,lu) )
      fpc_total       = sum( fpc_grid(:,jpngr), mask=islu(:,lu) )

      do pft=1,npft
        if (islu(pft,lu)) then
          if (ispresent(pft,jpngr)) then
            if (tree(pft)) then
              !-----------------------------------------------------------------------------
              ! TREES
              !-----------------------------------------------------------------------------
              if (fpc_tree_total>fpc_tree_max) then
                ! test kuno 2005-01-7: split excess according to total FPC rather than inc
                ! note that fpc_inc is already included in fpc_grid, fpc_tree_total, since
                ! it is used to calculate the excess fpc (check origin of these vars).
                ! kuno 2005-09-23:  mixed scheme is abolished, use only area-based scheme
                excess = (fpc_tree_total-fpc_tree_max)*fpc_grid(pft,jpngr)/fpc_tree_total
                                
                ! Reduce individual density (and thereby gridcell-level biomass)
                ! so that total tree FPC reduced to 'fpc_tree_max'
                nind_kill = nind(pft,jpngr) * (excess/fpc_grid(pft,jpngr))
                nind(pft,jpngr) = nind(pft,jpngr) - nind_kill

                ! Transfer lost biomass to litter
                call orgcpRec( orgfrac( nind_kill, pleaf(pft,jpngr) ), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )
                call orgcpRec( orgfrac( nind_kill, psapw(pft,jpngr) ), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )
                call orgcpRec( orgfrac( nind_kill, pwood(pft,jpngr) ), plitt_as(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )
                call orgcpRec( orgfrac( nind_kill, proot(pft,jpngr) ), plitt_bg(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )

              end if               ! fpc_tree_total > fpc_tree_max

            else
              !-----------------------------------------------------------------------------
              ! GRASSES
              !-----------------------------------------------------------------------------
              ! If grasses only, fpc_grass_total ~= 1; fpc_tree_total = 0
              if ( fpc_grass_total > ( 1.0 - min( fpc_tree_max, fpc_tree_total ) ) ) then           

                excess = ( fpc_grass_total - 1.0 + min( fpc_tree_max, fpc_tree_total ) ) * ( fpc_grid(pft,jpngr) / fpc_grass_total ) + 1.0d-15 

                cleaf_target = -1.0 / kbeer * dlog( 1.0d0 - ( fpc_grid(pft,jpngr) - excess ) ) / sla(pft)
                frac_kill = (pleaf(pft,jpngr)%c%c12 - cleaf_target) / pleaf(pft,jpngr)%c%c12

                ! Transfer lost biomass to litter
                call orgmvRec( orgfrac( frac_kill, pleaf(pft,jpngr)), pleaf(pft,jpngr), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )
                call orgmvRec( orgfrac( frac_kill, proot(pft,jpngr)), proot(pft,jpngr), plitt_bg(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )

              end if
            
            end if ! tree/grass
       
            !-----------------------------------------------------------------------------
            ! Update fractional plant coverage. Should now comply with fpc_tree_max and 
            ! not exceed 1.
            !-----------------------------------------------------------------------------
            call update_fpc(pft,jpngr)

          end if  ! ispresent

        end if  ! islu

      end do  !pft

    end do

    return

  end subroutine light

end module _dynveg
