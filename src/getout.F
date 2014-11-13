      subroutine getout_daily( jpngr )
      !////////////////////////////////////////////////////////////////
      !  SR called daily to sum up output variables.
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl, only: npft
      use statevars
      use outvars

      ! ARGUMENTS
      integer, intent(in) :: jpngr

      ! LOCAL VARIABLES
      integer :: pft

      do pft=1,npft

        ! Collect for annual output variables
        outagpp(pft,jpngr) = outagpp(pft,jpngr) + dgpp(pft)
        call ccp( dnpp(pft), outanpp(pft,jpngr) )
        call ncp( dnup(pft), outanup(pft,jpngr) )

      end do

      end subroutine getout_daily
