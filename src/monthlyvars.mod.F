       module monthlyvars
       !////////////////////////////////////////////////////////////////
       !  Module contains all monthly variables (with explicit dimension
       !  for month). For output.
        ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
        ! contact: b.stocker@imperial.ac.uk
       !----------------------------------------------------------------
         use _params_core

         implicit none
         type(carbon), dimension(nlu,maxgrid) :: mrh

       end module monthlyvars

