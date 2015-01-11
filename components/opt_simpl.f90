      ! To complile, run 
      !     gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan opt_simpl.f90 -o runopt
      !     pgf95 -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform opt_simpl.f90 -o runopt
      !     pgf90 opt_simpl.f90 -o runopt (works!)

      module hi_precision

        implicit none
        integer, parameter :: hi = selected_real_kind( 8 )

      end module hi_precision


      program opt_simpl
        !******************************************************************
        ! Fortran90 version of fmin pulled from                           !
        ! http://people.scs.fsu.edu/~burkardt/f_src/test_min/test_min.f90 !
        ! and modified slightly by Harry Valentine on 28 October 2008     !
        !******************************************************************
        use hi_precision
        !use numerical_recipe
        implicit none
        external funcmin

        !interface
        !  function funcmin(x)
        !    real, intent(in) :: x
        !    real :: funcmin
        !  end function funcmin
        !end interface

        !interface
        !  function fmin( a, b, funcmin, tol )
        !    real :: a, b
        !    real :: funcmin
        !    real, intent(in) :: tol
        !  end function fmin
        !end interface

        real( kind = hi ) :: fmin

        real( kind = hi ) :: tol
        real( kind = hi ) :: minimum

        ! tolerance level
        tol = 1.0D-12

        ! get minimum using 'fmin'
        minimum = fmin( -1.0d+00, 1.0d+00, funcmin, tol )

        print*,'minimum: ', minimum * 2.0d+00

      end program opt_simpl


      function funcmin( x )
        !******************************************************************
        ! Fortran90 version of fmin pulled from                           !
        ! http://people.scs.fsu.edu/~burkardt/f_src/test_min/test_min.f90 !
        ! and modified slightly by Harry Valentine on 28 October 2008     !
        !******************************************************************
        use hi_precision
        implicit none

        real( kind = hi ), intent(in) :: x
        real( kind = hi )             :: y
        real( kind = hi )             :: funcmin

        y = (-1.0d+00) * (x - 2.0d+00) ** 2 + 1.0d+03

        funcmin = - y

      end function funcmin


      function fmin( a, b, funcmin, tol )
        !******************************************************************
        ! Fortran90 version of fmin pulled from                           !
        ! http://people.scs.fsu.edu/~burkardt/f_src/test_min/test_min.f90 !
        ! and modified slightly by Harry Valentine on 28 October 2008     !
        !******************************************************************
        use hi_precision
        implicit none
        !external funcmin

        real( kind = hi ) :: a
        real( kind = hi ) :: b
        real( kind = hi ) :: c
        real( kind = hi ) :: d
        real( kind = hi ) :: e
        real( kind = hi ) :: eps
        real( kind = hi ) :: fu
        real( kind = hi ) :: fv
        real( kind = hi ) :: fw
        real( kind = hi ) :: fx
        real( kind = hi ) :: midpoint
        real( kind = hi ) :: funcmin
        real( kind = hi ) :: p
        real( kind = hi ) :: fmin
        real( kind = hi ) :: q
        real( kind = hi ) :: r
        real( kind = hi ) :: tol
        real( kind = hi ) :: tol1
        real( kind = hi ) :: tol2
        real( kind = hi ) :: u
        real( kind = hi ) :: v
        real( kind = hi ) :: w
        real( kind = hi ) :: x

        !  C is the squared inverse of the golden ratio.
        c = 0.5d+00 * ( 3.0d+00 - sqrt ( 5.0d+00 ) )

        !  EPS is the square root of the relative machine precision.
        eps = sqrt ( epsilon ( eps ) )

        !  Initialization.
        v = a + c * ( b - a )
        w = v
        x = v
        e = 0.0d+00
        fx = funcmin( x )
        fv = fx
        fw = fx

        !  The main loop starts here.
        do

          midpoint = 0.5d+00 * ( a + b )
          tol1 = eps * abs ( x ) + tol / 3.0d+00
          tol2 = 2.0d+00 * tol1

          !  Check the stopping criterion.
          if ( abs ( x - midpoint ) <= ( tol2 - 0.5d+00 * ( b - a ) ) ) then
            exit
          end if

          !  Is golden-section necessary?
          if ( abs ( e ) <= tol1 ) then
            if ( midpoint <= x ) then
              e = a - x
            else
              e = b - x
            end if

            d = c * e

          !  Consider fitting a parabola.
          else

            r = ( x - w ) * ( fx - fv )
            q = ( x - v ) * ( fx - fw )
            p = ( x - v ) * q - ( x - w ) * r
            q = 2.0d+00 * ( q - r )
            if ( 0.0d+00 < q ) then
              p = -p
            end if
            q = abs ( q )
            r = e
            e = d

            !  Choose a golden-section step if the parabola is not advised.
            if ( &
              ( abs ( 0.5d+00 * q * r ) <= abs ( p ) ) .or. &
              ( p <= q * ( a - x ) ) .or. &
              ( q * ( b - x ) <= p ) ) then

              if ( midpoint <= x ) then
                e = a - x
              else
                e = b - x
              end if

              d = c * e

            !  Choose a parabolic interpolation step.
            else

              d = p / q
              u = x + d

              if ( ( u - a ) < tol2 ) then
                d = sign ( tol1, midpoint - x )
              end if

              if ( ( b - u ) < tol2 ) then
                d = sign ( tol1, midpoint - x )
              end if

          end if

        end if

        !  F must not be evaluated too close to X.
        if ( tol1 <= abs ( d ) ) then
          u = x + d
        end if

        if ( abs ( d ) < tol1 ) then
          u = x + sign ( tol1, d )
        end if

        fu = funcmin( u )

        !  Update the data.
        if ( fu <= fx ) then

          if ( x <= u ) then
            a = x
          else
            b = x
          end if

          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
          cycle

        end if

        if ( u < x ) then
          a = u
        else
          b = u
        end if

        if ( fu <= fw .or. w == x ) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if ( fu <= fv .or. v == x .or. v == w ) then
          v = u
          fv = fu
        end if

        end do

        fmin = x

        return

      end function fmin


