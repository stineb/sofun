       module monthlyvars
       !////////////////////////////////////////////////////////////////
       !  Module contains all monthly variables (with explicit dimension
       !  for month). For output.
       !----------------------------------------------------------------
         use params_core

         implicit none
         type(carbon), dimension(nlu,maxgrid) :: mrh

       end module monthlyvars

