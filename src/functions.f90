! ===========================================================================
! functions.f90
! ===========================================================================
!> Weight functions of the particle-grid coupling: the B-splines W_n of order
!! n = 1 (linear, cloud in cell), 2 (quadratic) and 3 (cubic), in units of the
!! grid spacing. The same W_n deposits the mass on the grid and interpolates
!! the field back to the particles, which avoids a self-force. W_n has unit
!! area, vanishes for |y| >= (n+1)/2 and its weights on a uniform grid add up
!! to one.

module functions

  contains


  function Wn(n,y)

    implicit none

    integer n
    real(8) Wn,y


    if (n==1) then
      if (abs(y)<1.d0) then
        Wn = (1.d0-abs(y))
      else
        Wn = 0.d0
      end if

    else if (n==2) then

      if (abs(y)<0.5d0) then

        Wn = (0.75d0-y**2)

      else if (0.5<= abs(y) .and. abs(y)<1.5d0) then

        Wn = 0.125d0*(3.d0-2.d0*abs(y))**2

      else

        Wn = 0.d0

      end if

    else if (n==3) then

      if (abs(y)<1.d0) then

        Wn = 2.d0/3.d0 - y**2+abs(y)**3*0.5d0

      else if ((1.d0<=abs(y)) .and. (abs(y)<2.d0)) then

        Wn = 1.d0/6.d0*(2.d0-abs(y))**3

      else

        Wn = 0.d0

      end if


    else
       print *, "Weight function of order",n,"not implemented (valid orders are 1-3)"
       print *, "Aborting ..."
       stop 1
    end if
    return 
  end function Wn









  
end module functions
