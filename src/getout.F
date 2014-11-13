      subroutine getout_daily
      !////////////////////////////////////////////////////////////////
      !  SR called daily to sum up output variables.
      !  June 2014
      !  b.stocker@imperial.ac.uk
      !----------------------------------------------------------------
      use params_modl, only: npft
      use dailyvars
      use annvars

      ! LOCAL VARIABLES
      integer :: pft

      do pft=1,npft

        ! Collect for annual output variables
        agpp(pft) = agpp(pft) + dgpp(pft)
        call ccp( dnpp(pft), anpp(pft) )
        call ncp( dnup(pft), anup(pft) )

      end do

      end subroutine getout_daily