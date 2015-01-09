       module params_core
       !////////////////////////////////////////////////////////////////
       !  This module contains parameters that are not modified, but needed
       !  to define variables, dimension lengths, etc.
       !----------------------------------------------------------------
         implicit none
         
         integer, parameter :: recycle = 31             ! length of standard recycling period
         integer, parameter :: ndayyear = 365           ! number of days in a year
         integer, parameter :: nmonth = 12              ! number of months in a year
         integer, parameter :: maxlen_runname = 30      ! maximum length of runname (arbitrary)
         integer, parameter :: maxgrid = 1              ! number of spatial gridcells (dummy dimension for later code extension)
         integer, parameter :: nbucket = 2              ! number of buckets for soil water model
         integer, parameter :: npft = 1                 ! number of PFTs
         integer, parameter :: nlu = 1                  ! number of land units (tiles)
         integer, parameter :: lunat = 1                ! ID of natural land unit
         integer, parameter :: lucrop = 2               ! ID of crop land unit

         integer, parameter, dimension(nmonth) :: ndaymonth = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per month
         
         real, parameter    :: pi = 3.14159265359       ! pi - what else?
         real, parameter    :: dummy = -9999.           ! arbitrary dummy value

       end module params_core

