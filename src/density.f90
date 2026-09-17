
! ===========================================================================
! density.f90
! ===========================================================================
!> Mass density and radial current on the radial grid, deposited from the
!! particles:
!!
!!   rho(r)  = (1/r**2) 2 pi Int f L dp dL,    curr(r) = (1/r**2) 2 pi Int p f L dp dL,
!!
!! so that 4 pi Int r**2 rho dr = 8 pi**2 Int f L dr dp dL is the mass.
!! "density" deposits rho and curr with the shape Sn of width drc (the
!! particle cell) and avg_rho with the weight Wn of width dr (the grid);
!! "avg_density" computes only avg_rho, the density Poisson uses. Only avg_rho
!! conserves the mass on the grid when drc /= dr: rho is an output diagnostic.

subroutine density

  use parameters
  use arrays
  use functions
  use utils

  implicit none

  integer i,j
  real(8) :: smallpi,factor,average_rho,outside
  logical, save :: warned = .false.
  real(8) :: cutoff_rho,cutoff_avg
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)

  character(20) filename ! Name of output file.

  smallpi = acos(-1.0d0)

  factor = 2.0*smallpi*drc*dpc*dlc

  rho = 0.D0
  avg_rho = 0.D0
  curr = 0.D0

! The shape and weight functions have compact support, so a cell list
! (build_cell_list in utils.f90) restricts the deposit on grid point i to
! the particles filed in nearby cells.

  call build_cell_list(cell_start,particle_order)

  cutoff_rho = dble(bsplineorder)*drc
  cutoff_avg = dble(bsplineorder)*dr
  Wcell = ceiling(max(cutoff_rho,cutoff_avg)/dr) + 1

! Parallel over grid points only: rho(i), curr(i) and avg_rho(i) belong to
! one thread for the whole inner loop over particles.
  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j)

  do i=1,Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        if (abs(r(i)-r_part(j))<=cutoff_rho) then

          rho(i) = rho(i) + f(j)*l_part(j)*Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
          curr(i) = curr(i)+f(j)*l_part(j)*p_part(j)*Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
        end if

        if (abs(r(i)-r_part(j))<=cutoff_avg) then

          avg_rho(i) = avg_rho(i) + f(j)/(dr+dr**3/(12.d0*r(i)**2))*l_part(j)*Wn(bsplineorder,(r(i)-r_part(j))/dr)

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)

! Ghost points from the reflection symmetry f(r,p) = f(-r,-p), which makes
! rho even in r. The grid is staggered, r(i) = (i-1/2) dr, so the ghost point
! 1-k sits at -r(k) and mirrors the physical point k.

  do i=1,ghost
      rho(1-i) = rho(i)
      avg_rho(1-i) = avg_rho(i)
  end do


  rho = factor*m0*rho/r**2
  avg_rho = factor*m0*avg_rho/r**2

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

  use parameters
  use arrays
  use functions
  use utils

  implicit none

  integer i,j
  real(8) :: smallpi,factor,diff,contribution,shell
  real(8) :: cutoff_avg
  integer :: Wcell,c,clo,chi,pp
  integer, allocatable :: cell_start(:),particle_order(:)

  smallpi = acos(-1.0d0)

  factor = 2.0*smallpi*drc*dpc*dlc

  avg_rho = 0.D0

  call build_cell_list(cell_start,particle_order)

  cutoff_avg = (dble(bsplineorder) + 1.0d0)*dr

! Wn of order n vanishes for |y| >= (n+1)/2, and a particle is filed in the
! cell of its nearest grid point, so the particles that give grid point i a
! nonzero weight lie in cells i-Wcell..i+Wcell with Wcell = floor((n+2)/2)
! (1, 2, 2 for n = 1, 2, 3). Further cells only added exact zeros; leaving
! them out, and forming the shell denominator once per grid point, does not
! change the result.

  Wcell = (bsplineorder+2)/2

! Parallel over grid points only, as in density.

  !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(c,clo,chi,pp,j,diff,contribution,shell)

  do i = 1, Nr

    clo = max(1,i-Wcell)
    chi = min(Nr,i+Wcell)
    shell = r(i)**2*dr+dr**3/12.d0

    do c=clo,chi
      do pp=cell_start(c),cell_start(c+1)-1
        j = particle_order(pp)

        diff = abs(r(i) - r_part(j))
        if (diff <= cutoff_avg) then
            contribution = f(j) / shell *l_part(j)* Wn(bsplineorder, diff / dr)
            avg_rho(i) = avg_rho(i) + contribution
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  deallocate(cell_start,particle_order)

! Ghost points by reflection, as in density.

  do i=1,ghost
      avg_rho(1-i) = avg_rho(i)
  end do

  avg_rho = factor*m0*avg_rho

end subroutine avg_density
