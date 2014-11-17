      subroutine getout_daily( jpngr, mo, day )
      !////////////////////////////////////////////////////////////////
      !  SR called daily to sum up output variables.
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl, only: npft, lu_category
      use statevars
      use outvars
      use pools

      ! ARGUMENTS
      integer, intent(in) :: jpngr
      integer, intent(in) :: mo
      integer, intent(in) :: day

      ! LOCAL VARIABLES
      integer :: pft

      !print*,'netmin_litt ',netmin_litt(1)
      !print*,'netmin_soil ',netmin_soil(1)

      !----------------------------------------------------------------
      ! Collect daily output variables
      ! so far not implemented for isotopes
      !----------------------------------------------------------------
      do pft=1,npft
        lu=lu_category(pft)

        outdgpp(pft,day,jpngr)   = dgpp(pft)
        outdnpp(pft,day,jpngr)   = dnpp(pft)%c12
        outdnup(pft,day,jpngr)   = dnup(pft)%n14
        outdcleaf(pft,day,jpngr) = pleaf(pft,jpngr)%c%c12
        outdcex(pft,day,jpngr)   = dcex(pft)
        outdnetmin(lu,day,jpngr) = outdnetmin(lu,day,jpngr) + netmin_litt(pft)
      
      end do

      do lu=1,nlu
      
        outdnetmin(lu,day,jpngr) = outdnetmin(lu,day,jpngr) + netmin_soil(lu)
        outdninorg(lu,day,jpngr) = nh4(lu,jpngr)%n14 + no3(lu,jpngr)%n14

      end do

      !----------------------------------------------------------------
      ! Collect annual output variables
      !----------------------------------------------------------------
      do pft=1,npft

        outagpp(pft,jpngr) = outagpp(pft,jpngr) + dgpp(pft)
        outanpp(pft,jpngr) = outanpp(pft,jpngr) + dnpp(pft)%c12
        outanup(pft,jpngr) = outanup(pft,jpngr) + dnup(pft)%n14

      end do


      end subroutine getout_daily
