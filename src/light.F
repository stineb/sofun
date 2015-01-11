      subroutine light( jpngr, year )
        !////////////////////////////////////////////////////////////////
        ! LIGHT COMPETITION
        ! only effective with dynamic vegetation
        ! b.stocker@imperial.ac.uk
        !----------------------------------------------------------------

        implicit none

        ! ARGUMENTS
        integer, intent(in) :: jpngr, year


        ! LOCAL VARIABLES
        integer :: pft,lu
      
        real, dimension(nlu) :: fpc_tree_total  ! total grid FPC for tree PFTs
        real, dimension(nlu) :: fpc_grass_total ! total grid FPC for grass PFTs
        real, dimension(nlu) :: fpc_moss_total  ! total grid FPC for grass PFTs
        real, dimension(nlu) :: fpc_inc_tree    ! this years total FPC increment for tree PFTs
        real, dimension(nlu) :: fpc_grow_tree   ! this years total FPC increment for GROWING tree PFTs
        real :: excess                          ! tree FPC or grass cover to be reduced
        real :: nind_kill                       ! reduction in individual density to reduce tree FPC
                     !  to permitted maximum (indiv/m2)
        real :: lm_kill                         ! reduction in grass PFT leaf mass to reduce grass
                     !  cover to permitted maximum (gC)
        real :: rm_kill                         ! reduction in grass PFT root mass to reduce grass
        ! cover to permitted maximum (gC)
        real :: fpc_ind
        real :: lm_old

        ! Calculate total woody FPC, FPC increment and grass cover (= crown area)
        fpc_inc_tree(:) = 0.d0
        fpc_grow_tree(:) = 0.d0
        lm_old = 0.0d0      


        do lu=1,nlu

          !-----------------------------------------------------------------------------
          ! Get total area of trees and grassses (FPC = fractional plant cover)
          !-----------------------------------------------------------------------------
          fpc_tree_total  = sum( fpc_grid(:,jpngr), mask=tree(:).and.islu(:,lu) )
          fpc_grass_total = sum( fpc_grid(:,jpngr), mask=grass(:).and.islu(:,lu) )
          fpc_total       = sum( fpc_grid(:,jpngr), mask=islu(:,lu) )

          ! LIGHT COMPETITION
          ! reduce FPC of trees to 'fpc_tree_max', killing according number of
          ! individuals (=reducing 'nind')
          ! -------------------------------------------------------------------------        
          do pft=1,npft
            if (islu(pft,lu)) then
              if (ispresent(pft,jpngr)) then
                if (tree(pft)) then

                  if (fpc_tree_total>fpc_tree_max) then
                    ! test kuno 2005-01-7: split excess according to total FPC rather than inc
                    ! note that fpc_inc is already included in fpc_grid, fpc_tree_total, since
                    ! it is used to calculate the excess fpc (check origin of these vars).
                    ! kuno 2005-09-23:  mixed scheme is abolished, use only area-based scheme
                    excess = (fpc_tree_total-fpc_tree_max)*fpc_grid(pft,jpngr)/fpc_tree_total
                                    
                      
                    ! Reduce individual density (and thereby gridcell-level biomass)
                    ! so that total tree FPC reduced to 'fpc_tree_max'
                    nind_kill = nind(pft,jpngr) * (excess/fpc_grid(pft,jpngr))

                    nind(pft,jpngr) = nind(pft,jpngr)-nind_kill

                    ! Transfer lost biomass to litter
                    call addc(jpngr,litter_ag_fast(pft,jpngr,:),lm_ind(pft,jpngr,1)*nind_kill,lm_ind(pft,jpngr,2:ncvar))
                    call addc(jpngr, litter_ag_slow(pft,jpngr,:), sm_ind(pft,jpngr,1)*nind_kill, sm_ind(pft,jpngr,2:ncvar) )
                    call addc(jpngr, litter_ag_slow(pft,jpngr,:), hm_ind(pft,jpngr,1)*nind_kill, hm_ind(pft,jpngr,2:ncvar) )
                    call addc(jpngr, litter_bg(pft,jpngr,:), rm_ind(pft,jpngr,1)*nind_kill, rm_ind(pft,jpngr,2:ncvar) )
                    call addn ... 

                  end if               ! fpc_tree_total > fpc_tree_max

                else                  ! grass and moss

                  ! If grasses only, fpc_grass_total ~= 1; fpc_tree_total = 0
                  if ( fpc_grass_total > ( 1.0 - min( fpc_tree_max, fpc_tree_total ) ) ) then           
                    ! note: fpc_tree_max scales with size,
                    ! but fpc_tree_total is absolute in units of 1 cell area
                            
                    ! total grass (except LU), spahni 2009-09-16: including moss
                    excess = ( fpc_grass_total - 1.0 + min( fpc_tree_max, fpc_tree_total ) ) * ( fpc_grid(pft,jpngr) / fpc_grass_total ) + 1.0d-15 
                    lm_old = lm_ind(pft,jpngr,1)
                    lm_ind(pft,jpngr,1) = -1.0 / kbeer * dlog( 1.0d0 - ( fpc_grid(pft,jpngr) - excess ) ) / sla(pft)

                    ! kuno 2005-06-17: lm_ind can get VERY large (virtually infinite),
                    ! if whole area covered by grass...?  this may have the potential to
                    ! cause a FPE
                            
                    ! note - if only one grass type is allowed to grow, its fpc never
                    ! reaches one, thus the land area has no effect other than the
                    ! limitation of light coming in (fpar < 1). Only with several species
                    ! is a light competition code needed. THIS MAY BE TRUE FOR GRASSES ONLY!

                    lm_kill = lm_old - lm_ind(pft,jpngr,1)
                    rm_kill = rm_ind(pft,jpngr,1) * (lm_kill/lm_old)
                    rm_ind(pft,jpngr,1) = rm_ind(pft,jpngr,1) - rm_kill
                    
                    ! Transfer lost biomass to litter
                    call addc(jpngr, litter_ag_fast(pft,jpngr,:), lm_kill*nind(pft,jpngr), lm_ind(pft,jpngr,2:ncvar) )

      #if _DyN
                    ! Transfer the according amount of N to litter.
                    if (afntoc_lm(pft,jpngr).ne.MISSING_VALUE) then
                      call addn(jpngr,
           $               litter_ag_fast(pft,jpngr,:),
           $               lm_kill*afntoc_lm(pft,jpngr)*nind(pft,jpngr)
           $               )
                      lm_ind(pft,jpngr,n_id) = lm_ind(pft,jpngr,n_id)
           $               -lm_kill*afntoc_lm(pft,jpngr)
                    end if
      #end if
                    if (grass(pft)) then
                      call addc(jpngr,
           $               litter_bg(pft,jpngr,:),
           $               rm_kill*nind(pft,jpngr),
           $               rm_ind(pft,jpngr,2:ncvar)
           $               )
      #if _DyN
                      ! Transfer the according amount of N to litter.
                      if (afntoc_rm(pft,jpngr).ne.MISSING_VALUE) then
                        call addn(jpngr,
           $                 litter_bg(pft,jpngr,:),
           $                 rm_kill*afntoc_rm(pft,jpngr)*nind(pft,jpngr)
           $                 )
                      rm_ind(pft,jpngr,n_id) = rm_ind(pft,jpngr,n_id)
           $               -rm_kill*afntoc_rm(pft,jpngr)
                      end if
      #end if
                    else
                      ! spahni 2009-09-16: make sure mosses have no litter below ground
                      litter_bg(pft,jpngr,1) = 0.0d0
                      litter_bg(pft,jpngr,2:ncvar) = MISSING_VALUE
                      litter_bg_save(pft,jpngr,1) = 0.0d0
      #if _DyN
                      litter_bg(pft,jpngr,n_id) = 0.0d0
                      litter_bg_save(pft,jpngr,2) = 0.0d0
      #end if        
                    end if
                  end if
                
                end if
           
              end if                   ! present

              end if

          end do                     !pft

        ! Update FPC Grid after light competition
        do pft=1,npft
          if (present(pft,jpngr)) then
        ! update fpc (for establishment routine)
            call update_fpc(pft,jpngr)
  #if _competition_scheme_sg==0
        ! kuno 2005-01-31:  debug see allocation
            fpc_grid_old(pft,jpngr)=fpc_grid(pft,jpngr)
  #end if
          end if
        end do

        return

      end subroutine light

