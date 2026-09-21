! ===========================================================================
! grav_force.f90
! ===========================================================================
!> Potential and radial force on every particle: self-gravity (if
!! autointeraction), plus the background BGtype, plus the centrifugal term
!! L**2/(2 r**2) of the particle's own angular momentum.
!!
!! Backgrounds (unit mass and scale; force = -dpot/dr):
!!   Isochrone  pot = -1/(1+sqrt(1+r**2))
!!   Central    pot = -1/r
!!   sphere     constant density star of radius 1: pot = (r**2-3)/2 inside,
!!              -1/r outside
!!   iso, isotrun, nfw, burkert   see bgpot below
!!   null       none

subroutine grav_force



  use parameters
  use arrays
  use utils
  implicit none
  integer :: i
  real(8) :: sq,den      ! per-particle sqrt(1+r**2) and r**2+eps**2

! Self-gravity: solve the Poisson equation for the current particles.

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
         pot_part(i)   = pot_part(i) + (-1.0D0/abs(r_part(i)))
         force_part(i) = force_part(i) + (-sign(1.0D0,r_part(i))/r_part(i)**2)
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
       end do
       !$OMP END PARALLEL DO

       pot = pot + (-1.0D0/abs(r))
       force = force + (-sign(1.0D0,r)/r**2)

     else

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
       do i=1,Npart
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = (-1.0D0/abs(r_part(i)))
         force_part(i) = (-sign(1.0D0,r_part(i))/r_part(i)**2)
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
       print *
       print *, 'Unknown type of gravitational force'
       print *, 'Aborting ...'
       print *
       stop 1
     else if (.not. autointeraction) then
!      "null": no background, only the centrifugal term (plus self-gravity, if
!      any). With self-gravity poisson_rk has just set both arrays; without it
!      nothing has, and the centrifugal term below would be added on top of the
!      values left by the previous call, growing without bound
!      (AUDITORIA_2026-09-20.md, point 3).
       pot_part   = 0.0D0
       force_part = 0.0D0
     end if

     !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(den)
     do i=1,Npart
       den = r_part(i)**2 + eps*eps
       pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/den
       force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2
     end do
     !$OMP END PARALLEL DO

  end if

contains

! Backgrounds given by closed formulas for the potential and the force
! (force = -dpot/dr, checked symbolically for every one). Particles feel
! them always; with self-gravity they are added to the self potential and
! force, on the particles and on the grid, which only exists then.

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

! The closed forms are evaluated at |r|, and the force carries the sign of r.
! A particle may step to r < 0 inside a time step (main.f90 reflects it only
! afterwards) and the grid has ghost points at negative radii, so the
! background has to be even in r and its force odd. Written in terms of r
! itself, "nfw" and "burkert" are neither (nfw even reverses the sign of the
! force, pushing a particle that crosses the origin further out), "sphere"
! takes its inner branch for every r < 1 including r < -1, and "iso" is not
! defined at all for r < 0 (AUDITORIA_2026-09-20.md, point 4).

  elemental real(8) function bgpot(x)

    real(8), intent(in) :: x
    real(8) :: a

    a = abs(x)

    select case (BGtype)
    case ("sphere")
!      Constant density star of mass 1 and radius 1.
       if (a<1.d0) then
          bgpot = 0.5d0*(a**2 - 3.d0)
       else
          bgpot = - 1.d0/a
       end if
    case ("iso")
       bgpot = 3.0d0*log(a)
    case ("isotrun")
       bgpot = (10.0d0/6.0d0)*( 2.0d0*atan(a)/a + log(a**2+1) )
    case ("nfw")
       bgpot = -16.0d0*log( 1.0d0+a )/a
    case ("burkert")
       bgpot = ( 10.0d0/(3.0d0*a) )*( 2.0d0*(1.0d0+a)*atan(a) -2.0d0*(1.0d0+a)*log(1.0d0+a) &
               -(1.0d0-a)*log(1.0d0+a**2) )
    case default
       bgpot = 0.0d0
    end select

  end function bgpot

  elemental real(8) function bgforce(x)

    real(8), intent(in) :: x
    real(8) :: a

    a = abs(x)

    select case (BGtype)
    case ("sphere")
       if (a<1.d0) then
          bgforce = - a
       else
          bgforce = - 1.d0/a**2
       end if
    case ("iso")
       bgforce = -3.0d0/a
    case ("isotrun")
       bgforce = -(10.0d0/3.0d0)*( a-atan(a) )/a**2
    case ("nfw")
       bgforce = -16.0d0*( log(1.0d0+a)-a/(1.0d0+a) )/a**2
    case ("burkert")
       bgforce = -( 10.0d0/(3.0d0*a*a) )*( log( (1.0d0+a**2)*(1.0d0+a)**2 ) - 2.0d0*atan(a) )
    case default
       bgforce = 0.0d0
    end select

!   Odd extension to r < 0.

    if (x < 0.0d0) bgforce = - bgforce

  end function bgforce

end subroutine grav_force
