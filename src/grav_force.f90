
! **************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE   ***
! **************************************************

! Background field case.  Here we assume that the background
! gravitational field corresponds to the case of a constant
! density star of total mass 1 and radius 1, which has a
! gravitational potential "pot" given by:
!
! pot  =  1/2 ( r**2 - 3 )     r <  1
!
! pot  =  - 1 / r              r >= 1  
!
! for which the force is (force = - dpot/dr):
!
! force  = - r                 r <  1
! 
! force  = - 1 / r**2          r >= 1

subroutine grav_force



  use parameters
  use arrays
  use utils
  implicit none
  integer :: i
  real(8) :: smallpi
  real(8) :: sq,den      ! per-particle sqrt(1+r**2) and r**2+eps**2
  character(100) :: filename
  smallpi = acos(-1.0d0)

! Self-gravitating case.  In this case we need to
! solve the Poisson equation.

  if (autointeraction) then

     if (rmin/=0.0d0) then
        print *
        print *, 'For the self-gravitating case you should have rmin=0.'
        print *, 'Aborting ...'
        print *
        stop 1
     else
        call poisson_rk
!       Keep the self-gravity potential apart before the background and the
!       centrifugal barrier are added on top: the energy needs it with a
!       factor 1/2 (see energy.f90).
        potself_part = pot_part
     end if

  end if 

! Background, and the centrifugal term L**2/(2 r**2) of every particle.
! The isochrone and the point mass, the backgrounds used in production, are
! done in one parallel pass over the particles together with the centrifugal
! term, with sqrt(1+r**2) and r**2+eps**2 formed once per particle. As whole
! array expressions these were four serial passes and, once the force
! evaluation dominated the run, a large part of its cost. The expressions
! are the same, so the result is identical to the last bit.

  if (BGtype == "Isochrone") then

     if (autointeraction) then

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sq,den)
       do i=1,Npart
         sq  = sqrt(1.0D0+r_part(i)**2)
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = pot_part(i) + (-1.0D0/(1.0D0+sq))
         force_part(i) = force_part(i) + (-r_part(i)/(sq*(1.D0+sq)**2))
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
       end do
       !$OMP END PARALLEL DO

       pot = pot + (-1.0D0/(1.0D0+sqrt(1.0D0+r**2)))
       force = force + (-r/(sqrt(1.D0+r**2)*(1.D0+sqrt(1.D0+r**2))**2))

     else

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sq,den)
       do i=1,Npart
         sq  = sqrt(1.D0+r_part(i)**2)
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = (-1.0D0/(1.0D0+sq))
         force_part(i) = (-r_part(i)/sq*pot_part(i)**2)
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
       end do
       !$OMP END PARALLEL DO

     end if

  else if (BGtype == "Central") then

     if (autointeraction) then

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
       do i=1,Npart
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = pot_part(i) + (-1.0D0/r_part(i))
         force_part(i) = force_part(i) + (-1.0D0/r_part(i)**2)
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
       end do
       !$OMP END PARALLEL DO

       pot = pot + (-1.0D0/r)
       force = force + (-1.0D0/r**2)

     else

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
       do i=1,Npart
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = (-1.0D0/r_part(i))
         force_part(i) = (-1.0D0/r_part(i)**2)
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
       end do
       !$OMP END PARALLEL DO

     end if

  else

     if (BGtype == "sphere" .or. BGtype == "iso" .or. BGtype == "isotrun" .or. &
         BGtype == "nfw" .or. BGtype == "burkert") then
       call add_background
     else if (BGtype /= "null") then
!      "null": no background, only self-gravity, if any, and the centrifugal term.
       print *
       print *, 'Unknown type of gravitational force'
       print *, 'Aborting ...'
       print *
       stop 1
     end if

     !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
     do i=1,Npart
       den = r_part(i)**2 + eps*eps
       pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
       force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
     end do
     !$OMP END PARALLEL DO

  end if

        !pot   = pot + 0.5d0*Lfix**2/(r**2 + eps*eps)
        !force = force + Lfix**2*r/(r**2 + eps*eps)**2


  !filename = 'vlasov_potpart'
  !call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,pot_part)

  !filename = 'vlasov_comparepotential'
  !call save1Ddata(directory,filename,Nr,t,r,pot)

contains

! Backgrounds given by closed formulas for the potential and the force
! (force = -dpot/dr, checked symbolically for every one). Particles feel
! them always; with self-gravity they are added to the self potential and
! force, on the particles and on the grid (which only exists then).
! Before, "sphere" replaced the self-gravity instead of adding to it, and
! iso, isotrun, nfw and burkert only filled the grid arrays: particles
! never felt them, and without self-gravity those arrays were not even
! allocated.

  subroutine add_background

    if (autointeraction) then
       pot_part   = pot_part   + bgpot(r_part)
       force_part = force_part + bgforce(r_part)
       pot(1:Nr)   = pot(1:Nr)   + bgpot(r(1:Nr))
       force(1:Nr) = force(1:Nr) + bgforce(r(1:Nr))
    else
       pot_part   = bgpot(r_part)
       force_part = bgforce(r_part)
    end if

  end subroutine add_background

  elemental real(8) function bgpot(x)

    real(8), intent(in) :: x

    select case (BGtype)
    case ("sphere")
!      Constant density star of mass 1 and radius 1.
       if (x<1.d0) then
          bgpot = 0.5d0*(x**2 - 3.d0)
       else
          bgpot = - 1.d0/x
       end if
    case ("iso")
       bgpot = 3.0d0*log(x)
    case ("isotrun")
       bgpot = (10.0d0/6.0d0)*( 2.0d0*atan(x)/x + log(x**2+1) )
    case ("nfw")
       bgpot = -16.0d0*log( 1.0d0+x )/x
    case ("burkert")
       bgpot = ( 10.0d0/(3.0d0*x) )*( 2.0d0*(1.0d0+x)*atan(x) -2.0d0*(1.0d0+x)*log(1.0d0+x) &
               -(1.0d0-x)*log(1.0d0+x**2) )
    case default
       bgpot = 0.0d0
    end select

  end function bgpot

  elemental real(8) function bgforce(x)

    real(8), intent(in) :: x

    select case (BGtype)
    case ("sphere")
       if (x<1.d0) then
          bgforce = - x
       else
          bgforce = - 1.d0/x**2
       end if
    case ("iso")
       bgforce = -3.0d0/x
    case ("isotrun")
       bgforce = -(10.0d0/3.0d0)*( x-atan(x) )/x**2
    case ("nfw")
       bgforce = -16.0d0*( log(1.0d0+x)-x/(1.0d0+x) )/x**2
    case ("burkert")
       bgforce = -( 10.0d0/(3.0d0*x*x) )*( log( (1.0d0+x**2)*(1.0d0+x)**2 ) - 2.0d0*atan(x) )
    case default
       bgforce = 0.0d0
    end select

  end function bgforce

end subroutine grav_force
