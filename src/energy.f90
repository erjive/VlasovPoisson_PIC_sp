

subroutine energy


! ******************************************************
! ***   FIND DENSITY AND CURRENT IN PHYSICAL SPACE   ***
! ******************************************************
!
! This subroutine integrates over all phase space
! to the energy.  This quantity is defined as:
!
!            /
! Energy  =  | (p**2/2m + phi(r) ) f dr dp
!            /

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor

  smallpi = acos(-1.0d0)


! Zero Angular Momentum
!  if (Lfix == 0.0d0) then

!     factor = 1.0d0

! Include Angular Momentum
!  else

     factor = 8.D0*smallpi**2*drc*dpc*dlc

!  endif  

  total_energy = 0.0D0
  kinetic = 0.D0
  potential = 0.D0

  !$OMP PARALLEL DO SCHEDULE(GUIDED) REDUCTION (+:kinetic,potential,total_energy)
  do i=1,Npart
    kinetic   = kinetic + 0.5D0*p_part(i)**2*f(i)*l_part(i)
    potential = potential + pot_part(i)*f(i)*l_part(i)
    total_energy = total_energy + (0.5D0*p_part(i)**2+pot_part(i))*f(i)*l_part(i)
  end do
  !$OMP END PARALLEL DO

  kinetic   = factor * kinetic
  potential = factor * potential
  total_energy    = factor * total_energy

end subroutine energy
