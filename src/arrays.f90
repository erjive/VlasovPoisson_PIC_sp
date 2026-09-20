! ===========================================================================
! arrays.f90
! ===========================================================================
!> Arrays of the radial grid and of the particles.

  module arrays

  implicit none

! Radial grid (with ghost points at the origin when rmin=0).

  real(8), allocatable, dimension(:) :: r        ! Radial coordinate of the grid points.
  real(8), allocatable, dimension(:) :: force    ! Self-gravity (plus background) force on the grid, -dPhi/dr.
  real(8), allocatable, dimension(:) :: pot      ! Self-gravity (plus background) potential on the grid.
  real(8), allocatable, dimension(:) :: dev_pot  ! dPhi_self/dr on the grid.

! Weights of the mass quadrature of poisson_rk: the enclosed mass advances as
! M(r_i) = M(r_(i-1)) + mcoefA(i) rho(i-1) + mcoefB(i) rho(i), and the same
! weights give each particle its own contribution to M (the self-force).

  real(8), allocatable, dimension(:) :: mcoefA
  real(8), allocatable, dimension(:) :: mcoefB

  real(8), allocatable, dimension(:) :: rho      ! Mass density on the grid (output).
  real(8), allocatable, dimension(:) :: avg_rho  ! Mass density on the grid used by Poisson.
  real(8), allocatable, dimension(:) :: curr     ! Radial mass current on the grid (output).

! Particles.

  real(8), allocatable, dimension(:) :: r_part       ! Radius.
  real(8), allocatable, dimension(:) :: p_part       ! Radial momentum.
  real(8), allocatable, dimension(:) :: l_part       ! Angular momentum (conserved).
  real(8), allocatable, dimension(:) :: f            ! Phase-space weight.
  real(8), allocatable, dimension(:) :: r_part_p     ! Radius at the start of the step (euler, leapfrog).
  real(8), allocatable, dimension(:) :: p_part_p     ! Momentum at the start of the step (euler, leapfrog).
  real(8), allocatable, dimension(:) :: p_part_h     ! Momentum at the half step (leapfrog).
  real(8), allocatable, dimension(:) :: pot_part     ! Total potential at the particle, centrifugal term included.
  real(8), allocatable, dimension(:) :: potself_part ! Self-gravity part of pot_part (zero without autointeraction).
  real(8), allocatable, dimension(:) :: force_part   ! Total radial force on the particle.
  real(8), allocatable, dimension(:) :: q0_part, j0_part ! Initial angle and radial action (integrator analytic only).

  end module arrays
