  subroutine poisson_rk

! Erik: This subroutine was originally written by Miguel Alcubierre
! I adapted it for the vlasov_PIC code. 
! WARNING: This subroutine does not execute in parallel
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

  integer i,j,l

  real(8) spot,sdev_pot
  real(8) poth,dev_poth
  real(8) rho0,pi
  character(100) :: filename

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

  rho0 = 0.5d0*(avg_rho(0) + avg_rho(1))

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

! *******************************
! ***   ADD ANGULAR MOMENTUM  ***
! *******************************
!  if(Lfix /= 0.0d0) then
!     pot   = pot + 0.5d0*Lfix*Lfix/(r*r + eps*eps)
!     force = force + Lfix*Lfix*r/(r*r + eps*eps)**2
!  end if

! So far we have solved for the potential and the force 
! felt on the mesh. In order to calculate the force that 
! particles felt, we need to interpolate the potential 
! and force on each of them.

  pot_part   = 0.0D0
  force_part = 0.0D0

  !$OMP PARALLEL DO SCHEDULE(GUIDED) collapse(2)

  do i=1,Npart
    do j=1,Nr
      if (abs(r_part(i)-r(j))<=(bsplineorder+1)*dr) then

        pot_part(i)   = pot_part(i) + pot(j)*Wn(bsplineorder,(r_part(i)-r(j))/dr)

        force_part(i) = force_part(i) + force(j)*Wn(bsplineorder,(r_part(i)-r(j))/dr)

      end if
    end do
  end do
  !$OMP END PARALLEL DO

!  filename = 'vlasov_potpart'
!  call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,pot_part)

!     if (mod(l,spatial_output).eq.0) then

!  filename = 'vlasov_potpoisson'
!  call save1Ddata(directory,filename,Nr,t,r,pot)

!  filename = 'vlasov_potpoisson_r0'
!  call save0Ddata(directory,filename,t,pot(1))

!     end if



! ***************
! ***   END   ***
! ***************

  end subroutine poisson_rk


