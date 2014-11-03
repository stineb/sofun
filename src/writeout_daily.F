      subroutine writeout_daily(year,day)
      !/////////////////////////////////////////////////////////////////////////
      ! Write daily ASCII output
      !-------------------------------------------------------------------------
      use params_siml, only: firstyeartrend, spinupyears
      use params_core, only: ndayyear
      use dailyvars

      ! Arguments
      integer, intent(in) :: year
      integer, intent(in) :: day

      ! Local variables
      real :: itime

      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

      ! Write daily value, summed over all PFTs
      write(101,999) itime, dgpp

      return

999   format (F20.8,F20.8)

      end subroutine writeout_daily