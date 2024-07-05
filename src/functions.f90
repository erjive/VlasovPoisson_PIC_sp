! This module contains the shape and weight functions that represents 
! the macro-particles in the distribution function.
! 
! Based on the b-splines, the shape function of PIC methods
! is chosen as:
!
! S (x-xp) =   1  B  ((x-xp)/dx)
!  n          --   n
!             dx   
! Where B_n is the b-spline function of order n.
! The mathematical expressions were taken from SHARP code.

module functions

  contains

  function Sn(n,y,dx)

    implicit none
    integer n
    real(8) Sn,y,dx

    if (n==1) then

      if (abs(y) <= 0.5d0) then
        Sn = 1.d0
      else
        Sn = 0.d0
      end if

    else if (n==2) then

      if (abs(y)<1.0d0) then
        Sn = (1.d0-abs(y))
      else
        Sn = 0.d0
      end if

    else if (n==3) then

      if (abs(y)<=0.5d0) then

        Sn = (0.75d0-y**2)

      else if ((0.5<= abs(y)) .and. (abs(y)<1.5d0)) then

        Sn = 0.125d0*(3.d0-2.d0*abs(y))**2

      else

        Sn = 0.d0

      end if

    else if (n==4) then

      if (abs(y)<1.d0) then

        Sn = 2.d0/3.d0 - y**2+abs(y)**3*0.5d0

      else if ((1.d0<=abs(y)) .and. (abs(y)<2.d0)) then

        Sn = 1.d0/6.d0*(2.d0-abs(y))**3

      else

        Sn = 0.d0

      end if

    else if (n>4) then

       print *, "B-spline of order greater than 4 not implemented"
       print *, "Aborting ..."
       stop

    end if

    Sn = Sn/dx
    return 
  end function Sn


! The weight function W_n are just the S_(n+1) function multiplied by dx

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


    else if (n>3) then
       print *, "Weight function of order greater than 4 not implemented"
       print *, "Aborting ..."
       stop
    end if
    return 
  end function Wn


  function s1(y,dx)

    implicit none
    
    real(8) s1,y,dx  

    if (abs(y) <= 0.5d0) then
      s1 = 1.d0
    else
      s1 = 0.d0
    end if

      s1 = s1/dx

    return
  end function s1

  function s2(y,dx)

    implicit none
    
    real(8) s2,y,dx  

    if (abs(y)<1.d0) then
      s2 = (1.d0-abs(y))
    else
      s2 = 0.d0
    end if

    s2 = s2/dx

    return

  end function s2

  function s3(y,dx)

    implicit none
    
    real(8) s3,y,dx  

    if (abs(y)<=0.5d0) then

      s3 = (0.75d0-y**2)

    else if (0.5<= abs(y) .and. abs(y)<1.5d0) then

      s3 = 0.125d0*(3.d0-2.d0*abs(y))**2

    else

      s3 = 0.d0

    end if

    s3 = s3/dx

    return


  end function s3

  function s4(y,dx)

    implicit none
    
    real(8) s4,y,dx  

    if (abs(y)<1.d0) then

      s4 = 2.d0/3.d0 - y**2+abs(y)**3*0.5d0

    else if ((1.d0<=abs(y)) .and. (abs(y)<2.d0)) then

      s4 = 1.d0/6.d0*(2.d0-abs(y))**3

    else

      s4 = 0.d0

    end if

    s4 = s4/dx

    return


  end function s4


  function w1(y,dx)

    implicit none
    
    real(8) w1,y,dx  

    
    if (abs(y)<1.d0) then
      w1 = (1.d0-abs(y))
    else
      w1 = 0.d0
    end if

    return

  end function w1

  function w2(y,dx)

    implicit none
    
    real(8) w2,y,dx  

    if (abs(y)<=0.5d0) then

      w2 = (0.75d0-y**2)

    else if (0.5<= abs(y) .and. abs(y)<1.5d0) then

      w2 = 0.125d0*(3.d0-2.d0*abs(y))**2

    else

      w2 = 0.d0

    end if

    return

  end function w2

  function w3(y,dx)

    implicit none
    
    real(8) w3,y,dx  

    if (abs(y)<1.d0) then

      w3 = 2.d0/3.d0 - y**2+abs(y)**3*0.5d0

    else if ((1.d0<=abs(y)) .and. (abs(y)<2.d0)) then

      w3 = 1.d0/6.d0*(2.d0-abs(y))**3

    else

      w3 = 0.d0

    end if

    return

  end function w3
  
end module functions
