      program main

      implicit none

      interface
        subroutine print_to_screen(a,b,c,d,e)
        implicit none
        integer, intent(in) :: a, b, c
        integer, intent(in), optional :: d, e
        end subroutine print_to_screen
      end interface

      call print_to_screen(1,2,3)
      call print_to_screen(1,2,3,d=4)
      call print_to_screen(1,2,3,e=5)

      end program main


      subroutine print_to_screen(a,b,c,d,e)

      implicit none
      integer, intent(in) :: a, b, c
      integer, intent(in), optional :: d, e

      print*,'a,b,c',a,b,c
      if (present(d)) print*,'d',d
      if (present(e)) print*,'e',e

      end subroutine print_to_screen
