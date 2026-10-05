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
!!   iso, isotrun, nfw, burkert   see bgpot in utils.f90
!!   null       none

subroutine grav_force



  use parameters
  use arrays
  use utils
  implicit none
  integer :: i
  real(8) :: sq,den      ! per-particle sqrt(1+r**2) and r**2+eps**2
  real(8), parameter :: dmin = tiny(1.0d0)

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
!
! The two isochrone branches also write the force the same way. They used to
! differ: the branch without self-gravity reused the potential it had just
! computed, -r/sq*pot**2, instead of -r/(sq*(1+sq)**2). The two agree in
! exact arithmetic but not in floating point -- they differ in the last bit
! for 222113 of 400002 radii tested -- so a run with negligible mass and
! autointeraction=.true. did not reproduce the run without it, which is
! precisely the control this audit leans on (AUDITORIA_2026-09-20.md,
! point 15).
!
! The denominators of the centrifugal term are bounded below by the smallest
! normal number. The term vanishes for a particle with L = 0, but written as
! L**2/(2 r**2) and L**2 r/r**4 it evaluated 0/0 when such a particle sat
! exactly at r = 0, where the isochrone and the sphere are regular
! (pot = -1/2 and -3/2, no force); the NaN then reached the total energy and
! every h_k of the run (AUDITORIA_FISICA_2026-09-21.md, N2; test U7). With
! the bound it is 0/dmin = 0, and nothing changes for r**4 above 2.2e-308.
! A particle with L > 0 exactly at r = 0, which the dynamics cannot reach,
! gets a huge potential and no centrifugal force instead of NaN.
! A test on L inside these loops, with "if" or with "merge", was tried
! first: the loops stop being vectorized and a run without self-gravity
! takes 1.7 times as long.

  if (BGtype == "Isochrone") then

     if (autointeraction) then

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sq,den)
       do i=1,Npart
         sq  = sqrt(1.0D0+r_part(i)**2)
         den = r_part(i)**2 + eps*eps
         pot_part(i)   = pot_part(i) + (-1.0D0/(1.0D0+sq))
         force_part(i) = force_part(i) + (-r_part(i)/(sq*(1.D0+sq)**2))
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/max(den,dmin)
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/max(den**2,dmin)
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
         force_part(i) = (-r_part(i)/(sq*(1.D0+sq)**2))
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/max(den,dmin)
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/max(den**2,dmin)
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
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/max(den,dmin)
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/max(den**2,dmin)
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
         pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/max(den,dmin)
         force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/max(den**2,dmin)
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
       pot_part(i)   = pot_part(i) + 0.5d0*l_part(i)**2/max(den,dmin)
       force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/max(den**2,dmin)
     end do
     !$OMP END PARALLEL DO

  end if

contains

! Backgrounds given by closed formulas for the potential and the force,
! bgpot and bgforce of utils.f90. Particles feel them always; with
! self-gravity they are added to the self potential and force, on the
! particles and on the grid, which only exists then.

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

end subroutine grav_force
