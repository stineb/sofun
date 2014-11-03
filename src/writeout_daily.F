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

      !/////////////////////////////////////////////////////////////////////////
      ! Collect variables to output variables
      !-------------------------------------------------------------------------
      do lu=1,nlu
        dreco(lu) = dreco(lu) + drhet(lu)
        do pft=1,npft
          if (lu==lu_category(pft)) then
            dreco(lu) = dreco(lu) + drsoil(pft) + drauto(pft) + dcex(pft)
          end if
        end do
      end do

      !/////////////////////////////////////////////////////////////////////////
      ! Write daily value, summed over all PFTs
      ! xxx implement taking sum over PFTs in this land use category
      !-------------------------------------------------------------------------
      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

      write(101,999) itime, dgpp

      return

999   format (F20.8,F20.8)

      end subroutine writeout_daily