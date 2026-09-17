! ===========================================================================
! density.f90
! ===========================================================================
!> Mass density and radial current on the radial grid, deposited from the
!! particles with the weight function W_n of order n (bsplineorder) and width
!! dr, the same function that interpolates the field back to the particles.
!!
!! A particle is a thin shell of mass m_j = 8 pi**2 drc dpc dlc f_j L_j at
!! radius r_j. Grid point i receives m_j W_n((r_i-r_j)/dr), and its density is
!! that mass over the volume the weight covers,
!!
!!   V_i = Int W_n((r_i-r)/dr) 4 pi r**2 dr = 4 pi dr (r_i**2 + dr**2 (n+1)/12),
!!
!! since W_n has unit area and second moment (n+1)/12 in units of dr**2.
!! With this volume a uniform density is reproduced exactly and
!! Sum_i rho_i V_i is the deposited mass.
!!
!! With rmin=0 the distribution obeys f(r,p) = f(-r,-p): every particle has
!! an image at (-r_j,-p_j), and near the origin its weight on grid point i,
!! W_n((r_i+r_j)/dr), is added (for the current, with the opposite sign
!! because p changes sign). This is the same symmetry that fills the ghost
!! points, and it makes the formula for V_i hold down to the first point.

subroutine density

  use parameters
  use arrays
  use utils

  implicit none

  integer j
  real(8) :: smallpi,average_rho,outside
  logical, save :: warned = .false.
  character(20) filename ! Name of output file.

  smallpi = acos(-1.0d0)

  call deposit(.true.)

  rho = avg_rho

! With self-gravity, Poisson only sees the mass on the grid. Report, once,
! when part of it has left.

  if (autointeraction .and. .not. warned) then
    outside = 0.0D0
    do j=1,Npart
      if (r_part(j) > r(Nr)) outside = outside + f(j)*l_part(j)
    end do
    if (outside > 0.0D0) then
      print *
      print *, 'WARNING: at t = ',t,' a fraction ',outside/sum(f*l_part), &
               ' of the mass is beyond r = ',r(Nr)
      print *, '         and does not enter the Poisson equation.'
      warned = .true.
    end if
  end if

! Mean density in the shell r1 <= r <= r2, written to vlasov_rhomix.tl: the
! mass of the particles inside, 8 pi**2 drc dpc dlc Sum f L, over the volume
! 4 pi (r2**3 - r1**3)/3.

  average_rho = 0.D0

  do j=1,Npart
    if (r_part(j)>=r1 .and. r_part(j)<= r2) then
      average_rho = average_rho + f(j)*l_part(j)
    end if
  end do

  average_rho = 8.0D0*smallpi**2*drc*dpc*dlc*average_rho &
              / (4.0D0*smallpi*(r2**3-r1**3)/3.0D0)

  filename = 'vlasov_rhomix'
  call save0Ddata(directory,filename,t,average_rho)

end subroutine density


!> avg_rho alone, for Poisson; called on every force evaluation with
!! self-gravity, so it is the most expensive loop of those runs.

subroutine avg_density

  implicit none

  call deposit(.false.)

end subroutine avg_density


!> Deposit the particles on the grid: avg_rho always, curr if want_curr.

subroutine deposit(want_curr)

  use parameters
  use arrays
  use functions
  use utils

  implicit none

  logical, intent(in) :: want_curr

  integer i,j
  real(8) :: smallpi,factor,cutoff_w,vol,wd,wi,m
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)
  logical :: images

  smallpi = acos(-1.0d0)

! Mass of a particle is 4 pi factor f L; density = mass W / V_i.
  factor = 2.0d0*smallpi*drc*dpc*dlc

  avg_rho = 0.D0
  if (want_curr) curr = 0.D0

  images = (ghost > 0)

  call build_cell_list(cell_start,particle_order)

! W_n vanishes for |y| >= (n+1)/2. A particle is filed in the cell of its
! nearest grid point, so the particles that give grid point i a nonzero
! weight lie in cells i-Wcell..i+Wcell with Wcell = floor((n+2)/2). An image
! at -r_j only reaches the first points, whose cell range already includes
! the cells near the origin where such a particle is filed.

  cutoff_w = 0.5d0*dble(bsplineorder+1)*dr
  Wcell = (bsplineorder+2)/2

! Parallel over grid points only: avg_rho(i) and curr(i) belong to one thread
! for the whole inner loop over particles.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,vol,wd,wi,m)

  do i = 1, Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)
    vol = r(i)**2*dr + dr**3*dble(bsplineorder+1)/12.0d0

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        wd = 0.0d0
        wi = 0.0d0
        if (abs(r(i)-r_part(j)) < cutoff_w) wd = Wn(bsplineorder,(r(i)-r_part(j))/dr)
        if (images .and. abs(r(i)+r_part(j)) < cutoff_w) wi = Wn(bsplineorder,(r(i)+r_part(j))/dr)

        if (wd /= 0.0d0 .or. wi /= 0.0d0) then
          m = f(j)*l_part(j)/vol
          avg_rho(i) = avg_rho(i) + m*(wd + wi)
          if (want_curr) curr(i) = curr(i) + m*p_part(j)*(wd - wi)
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)

! Ghost points: the density is even in r, the radial current odd.

  do i=1,ghost
      avg_rho(1-i) = avg_rho(i)
      if (want_curr) curr(1-i) = -curr(i)
  end do

  avg_rho = factor*avg_rho
  if (want_curr) curr = factor*curr

end subroutine deposit
