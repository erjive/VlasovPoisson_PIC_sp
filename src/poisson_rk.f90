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

  real(8) spot,sdev_pot
  real(8) poth,dev_poth
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

! Solve equation using second order Runge-Kutta.

! First calculate the density
  call avg_density


! Note 1: The first point must be treated differently
! since we have a division by zero at the origin.
! If we take  dev_pot~a*r  for r<<1, with "a" constant,
! one can show that a=(4/3)*pi*rho(r=0)
! (but notice that this is only second order).

! Note 2: Remember that in order to apply correctly
! the boundaty conditions, we have two ghost points
! to the left of the origin, so the first point with
! positive r is i=1.

! First point with positive r (i=1).

  rho0 = avg_rho(1)

  pot(1) = 2.d0/3.d0*pi*rho0*r(1)**2
  dev_pot(1) = 4.d0/3.d0*pi*rho0*r(1)

! All other points.

  do i=2,Nr

!    Calculate sources at left point.

     spot = dev_pot(i-1)
     sdev_pot = - 2.d0*dev_pot(i-1)/r(i-1) + 4.d0*pi*avg_rho(i-1)

!    Advance half a step.

     poth = pot(i-1) + 0.5d0*dr*spot
     dev_poth = dev_pot(i-1) + 0.5d0*dr*sdev_pot

!    Calculate sources at intermediate point.

     spot = dev_poth
     sdev_pot = - 4.d0*dev_poth/(r(i-1) + r(i)) &
          + 2.d0*pi*(avg_rho(i-1) + avg_rho(i))

!    Advanced full step.

     pot(i) = pot(i-1) + dr*spot
     dev_pot(i) = dev_pot(i-1) + dr*sdev_pot

  end do

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


