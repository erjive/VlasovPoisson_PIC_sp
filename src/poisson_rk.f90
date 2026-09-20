  subroutine poisson_rk

! Originally written by Miguel Alcubierre; adapted to this code by Erik.
! *******************
! ***   POISSON   ***
! *******************

! This subroutine solves the Poisson equation
! for the self-gravitating case:
!
! Laplacian (pot)  =  4 pi rho
!
! with pot the gravitational field and rho the mass
! density.  Once we have pot, the gravitational
! force is calculated as:
!
! force = - grad (pot)
!
! In spherical symmetry the above equations reduce to:
!
!  2
! d pot  +  (2/r) d pot  =  4 pi rho
!  r               r
!
! force  =  - d pot /dr

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j

  real(8) ra,rb,slope,c0,c3,c4
  real(8) rho0,pi
  real(8) cutoff_interp
  integer :: Wgrid,jc,m
  real(8) :: w,rj,pj,fj,rm,sgn

! *******************
! ***   NUMBERS   ***
! *******************

  pi = acos(-1.d0)


! ****************************************
! ***   FIND GRAVITATIONAL POTENTIAL   ***
! ****************************************

! Initialize arrays to zero.

  pot = 0.0d0
  dev_pot = 0.0d0

! The equation is integrated in the enclosed mass rather than in dPhi/dr,
!
!   dM/dr = 4 pi r**2 rho,      dPhi/dr = M/r**2,
!
! which removes the 2/r term and with it the only source of stiffness at the
! origin. Between two grid points the density is the straight line through
! its two values,
!
!   rho(x) = rho(i-1) + slope (x - ra),   slope = (rho(i)-rho(i-1))/dr,
!
! so M is a quartic in x with no linear or quadratic term,
!
!   M(x) = c0 + c3 x**3 + c4 x**4,   c3 = 4 pi (rho(i-1) - slope ra)/3,
!                                    c4 = pi slope,
!
! and both integrals are done in closed form, with no quadrature error:
!
!   Int M/x**2 dx = -c0/x + c3 x**2/2 + c4 x**3/3.
!
! For the linear weight (bsplineorder = 1) the straight line between grid
! points reproduces the deposited profile exactly, so every particle
! contributes its whole mass to the field whatever its radius, including the
! first cell. The Runge-Kutta this file is named after integrated dPhi/dr
! directly and lost that mass near the origin: a particle sitting on r(1)
! produced no field at all (see AUDITORIA_2026-09-20.md, point 2).

! First calculate the density
  call avg_density

! Between the origin and r(1) the density is even, so the ghost point carries
! the same value as r(1) and rho is constant there:
!
!   M(r1) = (4/3) pi rho(1) r1**3,   Phi(r1) - Phi(0) = (2/3) pi rho(1) r1**2.
!
! Phi(0) = 0 is arbitrary; the additive constant is fixed at the end. Note 2
! of the original code still holds: with rmin = 0 there are two ghost points
! to the left of the origin and the first point with positive r is i=1.

  rho0 = avg_rho(1)

! dev_pot carries M while the integration runs, and is divided by r**2 at the
! end. It is a local role: outside this loop dev_pot is always dPhi/dr.

  dev_pot(1) = 4.d0/3.d0*pi*rho0*r(1)**3
  pot(1)     = 2.d0/3.d0*pi*rho0*r(1)**2

  do i=2,Nr

     ra = r(i-1)
     rb = r(i)

     slope = (avg_rho(i) - avg_rho(i-1))/dr

     c3 = 4.d0*pi*(avg_rho(i-1) - slope*ra)/3.d0
     c4 = pi*slope
     c0 = dev_pot(i-1) - c3*ra**3 - c4*ra**4

     dev_pot(i) = c0 + c3*rb**3 + c4*rb**4

     pot(i) = pot(i-1) + (- c0/rb + c3*rb**2/2.d0 + c4*rb**3/3.d0) &
                       - (- c0/ra + c3*ra**2/2.d0 + c4*ra**3/3.d0)

  end do

! From the enclosed mass to dPhi/dr.

  dev_pot(1:Nr) = dev_pot(1:Nr)/r(1:Nr)**2

! Ghost points using symmetries.

  pot(-1) = pot(2)
  pot( 0) = pot(1)

  dev_pot(-1) = - dev_pot(2)
  dev_pot( 0) = - dev_pot(1)

! Calculate the force.

  force = - dev_pot

! We now need to substract a constant to the solution
! to make sure that far away the potential behaves as
! pot ~ 1/r.  This condition implies that we should
! have pot + r*dev_pot = 0 far away.  But we won't since
! we arbitrarily fixed pot(r=0)=0 above. So now we
! find the value of pot + r*dev_pot at the outer boundary
! and just substract it from the whole solution.

  pot = pot - (pot(Nr) - force(Nr)*r(Nr))

! The centrifugal term depends on each particle's L, so it is added per
! particle in grav_force, not on the grid.

! So far we have solved for the potential and the force 
! felt on the mesh. In order to calculate the force that 
! particles felt, we need to interpolate the potential 
! and force on each of them.

  pot_part   = 0.0D0
  force_part = 0.0D0

! Interpolation to the particles with the same weight W_n used in the
! deposit, over every grid point j within its support. The grid is uniform,
! r_j = r(1) + (j-1) dr, and the field is known beyond the stored points:
!
!   j <= 0   the mirror of point 1-j (r_j = -r_(1-j)): the potential is even
!            and the force odd, the symmetry f(r,p) = f(-r,-p) at the origin;
!   beyond   a point past r(Nr), stored or mirrored, takes the exterior
!            solution Phi = Phi(r_Nr) r_Nr/r and F = F(r_Nr) (r_Nr/r)**2 of
!            the mass on the grid.
!
! So the weights of every particle add up to one, at any radius, including
! particles closer to the origin than r(1) or beyond the last point. The loop
! runs in parallel over particles.

  Wgrid = (bsplineorder+2)/2
  cutoff_interp = 0.5d0*dble(bsplineorder+1)*dr

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,jc,rj,w,pj,fj,m,rm,sgn)

  do i=1,Npart

    jc  = nint((r_part(i)-r(1))/dr) + 1

    do j=jc-Wgrid,jc+Wgrid

      rj = r(1) + dble(j-1)*dr
      if (abs(r_part(i)-rj) >= cutoff_interp) cycle

!     Point j <= 0 is the mirror of point m = 1-j (potential even, force
!     odd); a point beyond the grid, stored or mirrored, takes the exterior
!     solution.
      if (j <= 0) then
        m = 1-j
        sgn = -1.0d0
      else
        m = j
        sgn = 1.0d0
      end if
      if (m > Nr) then
        rm = r(1) + dble(m-1)*dr
        pj = pot(Nr)*r(Nr)/rm
        fj = sgn*force(Nr)*(r(Nr)/rm)**2
      else
        pj = pot(m)
        fj = sgn*force(m)
      end if

      w = Wn(bsplineorder,(r_part(i)-rj)/dr)

      pot_part(i)   = pot_part(i) + pj*w
      force_part(i) = force_part(i) + fj*w

    end do
  end do
  !$OMP END PARALLEL DO

  end subroutine poisson_rk


