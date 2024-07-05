! *************************
! ***   MODULE ARRAYS   ***
! *************************

  module arrays

  implicit none

! The coordinates (r,p), plus the force and potential are 1D arrays.

  real(8), allocatable, dimension(:) :: r        ! Radial coordinate r.
  real(8), allocatable, dimension(:) :: p        ! Momentum coordinate p.
  real(8), allocatable, dimension(:) :: force    ! Gravitational force (F=-dphi/dr).
  real(8), allocatable, dimension(:) :: pot      ! Gravitational potential.
  real(8), allocatable, dimension(:) :: dev_pot  ! Derivative of gravitational potential (dphi/dr)

! Position, momentum, phase space density, potential and force at the particles

  real(8), allocatable, dimension(:) :: r_part     ! Position of the particles.
  real(8), allocatable, dimension(:) :: r_part_p   ! Old position of the particles.
  real(8), allocatable, dimension(:) :: p_part     ! Momentum of the particles.
  real(8), allocatable, dimension(:) :: p_part_h   ! Momentum of the particles at half time step.
  real(8), allocatable, dimension(:) :: p_part_p   ! Old momentum of the particles.
  real(8), allocatable, dimension(:) :: p_part_hp  ! Old momentum of the particles at half time step.
  real(8), allocatable, dimension(:) :: l_part     ! Momentum of the particles.
  real(8), allocatable, dimension(:) :: f          ! Phase space density.
!  real(8), allocatable, dimension(:) :: f_p          ! Phase space density.
  real(8), allocatable, dimension(:) :: pot_part   ! Potential valuated at the position particle
  real(8), allocatable, dimension(:) :: force_part ! Force applied to the particle.


! The density function f and the fluxes are 2D arrays.


!  real(8), allocatable, dimension(:,:) :: f_p    ! Old value of f.
!  real(8), allocatable, dimension(:,:) :: sf     ! Source for phase space density.
!  real(8), allocatable, dimension(:,:) :: flux_r ! flux in r direction:  p*f
!  real(8), allocatable, dimension(:,:) :: flux_p ! flux in p direction:  force*f

! The integrated particle density, current and continuity euqstion are 1D arrays.

  real(8), allocatable, dimension(:) :: rho      ! Particle density in physical space.
  real(8), allocatable, dimension(:) :: rho_p    ! Old value of rho.
  real(8), allocatable, dimension(:) :: avg_rho  ! Average over the cells of the density
  real(8), allocatable, dimension(:) :: curr     ! Particle current in physical space.
  real(8), allocatable, dimension(:) :: curr_p   ! Old value of the current.
  real(8), allocatable, dimension(:) :: cont     ! Continuity equation (should converge to 0).

  real(8), allocatable, dimension(:,:) :: res   ! Residual evaluator.


  end module arrays

