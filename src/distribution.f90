! ===========================================================================
! distribution.f90
! ===========================================================================
!> Initial distribution function in angle-action variables, F0(Q,J,L), for
!! state "aa_quad".
!!
!! The parameter "dftype" chooses which one is used, so the rest of the code
!! never mentions a particular functional form: the quadrature state aa_quad
!! in initial_data.f90 only asks for F0 at each node.
!!
!! All of them are written with the same convention for a width, exp(-x^2/s^2),
!! and are evaluated on Q in [0,2pi) and J >= 0. Every profile in (Q,J) is
!! multiplied by the same factor in angular momentum, C(L) = exp(-(L-l0)^2/sl^2);
!! "king" also depends on L through the energy E(J,L).

module distribution

  use parameters

  implicit none

! --- "bimodal": two groups of actions, three angular harmonics ---
!
!   F0 = [1 + b1 cos Q + b2 cos 2Q] * [exp(-(J-Ja)^2/sa^2) + wb exp(-(J-Jb)^2/sb^2)]
!
! The two groups sit at different J, so they mix at different rates and h_1
! beats at their frequency difference. The angular part is a finite Fourier
! series: it has no harmonic above k=2, so h_3 and h_4 must vanish exactly.

  real(8), parameter :: bim_b1 = 0.8d0    ! amplitude of cos(Q)
  real(8), parameter :: bim_b2 = 0.3d0    ! amplitude of cos(2Q)
  real(8), parameter :: bim_Ja = 0.10d0   ! first group: centre in J
  real(8), parameter :: bim_sa = 0.025d0  !              width in J
  real(8), parameter :: bim_Jb = 0.24d0   ! second group: centre in J
  real(8), parameter :: bim_sb = 0.035d0  !               width in J
  real(8), parameter :: bim_wb = 0.7d0    ! weight of the second group

! --- "spiral": not separable, wound up already at t=0 ---
!
!   F0 = exp(-(J-J0)^2/s0^2) * exp(-sin((Q - beta J)/2)^2/sq^2)
!
! The angular profile is centred on Q = beta*J, so the distribution starts as
! a spiral in the (Q,J) plane. Phase mixing carries the centre to
! Q = beta*J + omega(J) t, and the spread of that across J closes at
! t = -beta/omega'(J0) > 0, where the distribution unwinds and |h_k| peaks
! before decaying for good. With a spread in L the frequency also varies
! across L at fixed J, and that spread is not wound up at t=0, so it mixes
! the spiral away independently of the unwinding. Computed with
! tools/hk_exacto.py for l0=2, t* ~ 710: |h_1(t*)|/|h_1(0)| is 2.2 for L
! fixed, 2.0 with sl=0.02 and 2.2e-3 with sl=0.2 on [1.6,2.4], where the
! edges of the truncated support leave an algebraic tail. Use a narrow L
! range to see the unwinding itself.

  real(8), parameter :: spi_J0   = 0.15d0  ! centre in J
  real(8), parameter :: spi_s0   = 0.04d0  ! width in J
  real(8), parameter :: spi_beta = 50.0d0  ! winding of the initial spiral
  real(8), parameter :: spi_sq   = 0.5d0   ! width in Q

! --- "king": equilibrium in J with a single-harmonic perturbation ---
!
!   F0 = f_eq(J) * [1 + eps cos Q],
!   f_eq(J) = exp(-E(J)/sE2) - exp(-E(Jt)/sE2)   for J < Jt,  0 beyond,
!
! with E(J) = -1/(2(J+c)^2) the energy of the radial orbit. This is the
! lowered isothermal (King) profile written in action space, perturbed by one
! angular harmonic: the textbook setting for phase mixing of a linear
! perturbation. Unlike the other two it is not smooth, it ends in a corner at
! J = Jt, which is what makes its quadrature converge algebraically instead of
! faster than any power.

  real(8), parameter :: kin_Jt  = 0.35d0   ! action where the profile is cut
  real(8), parameter :: kin_sE2 = 0.01d0   ! energy scale (sigma_E^2)
  real(8), parameter :: kin_eps = 0.5d0    ! amplitude of the cos(Q) perturbation

 contains

!> Energy of the orbit of radial action J and angular momentum L in the isochrone.

  function energy_of_J(J,L) result(E)

    implicit none

    real(8), intent(in) :: J,L
    real(8) :: E

    E = -1.0d0/(2.0d0*(J + 0.5d0*(L+sqrt(L**2+4.0d0)))**2)

  end function energy_of_J


!> Lowered isothermal profile in action space, zero beyond the cut.

  function king_feq(J,L) result(feq)

    implicit none

    real(8), intent(in) :: J,L
    real(8) :: feq

    if (J >= kin_Jt) then
       feq = 0.0d0
    else
       feq = exp(-energy_of_J(J,L)/kin_sE2) - exp(-energy_of_J(kin_Jt,L)/kin_sE2)
    end if

  end function king_feq

!> F0 at one point of (Q,J,L).

  function df0(Q,J,L) result(F)

    implicit none

    real(8), intent(in) :: Q,J,L
    real(8) :: F

    select case (dftype)

    case ('gauss')
!      The expression of the original state "aa", factor by factor.
       F = exp(-sin(0.5d0*Q)**2/sp**2)*exp(-J**2/sr**2)*J**2*exp(-(L-l0)**2/sl**2)

    case ('bimodal')
       F = (1.0d0 + bim_b1*cos(Q) + bim_b2*cos(2.0d0*Q)) &
           *(exp(-(J-bim_Ja)**2/bim_sa**2) + bim_wb*exp(-(J-bim_Jb)**2/bim_sb**2)) &
           *exp(-(L-l0)**2/sl**2)

    case ('spiral')
       F = exp(-(J-spi_J0)**2/spi_s0**2) &
           *exp(-sin(0.5d0*(Q-spi_beta*J))**2/spi_sq**2)*exp(-(L-l0)**2/sl**2)

    case ('king')
       F = king_feq(J,L)*(1.0d0 + kin_eps*cos(Q))*exp(-(L-l0)**2/sl**2)

    case default
       F = 0.0d0
       call df0_unknown

    end select

  end function df0


  subroutine df0_unknown

    implicit none

    print *
    print *, 'Unknown dftype: ',trim(dftype)
    print *, 'Aborting ...'
    print *
    stop 1

  end subroutine df0_unknown


!> Window in J that contains the whole distribution, for the quadrature grid.
!! Each profile declares its own, so the quadrature no longer borrows the width
!! of one particular distribution.

  subroutine df0_Jrange(Jlo,Jhi)

    implicit none

    real(8), intent(out) :: Jlo,Jhi

!   Six widths past the last feature leaves exp(-36) of the peak outside, which
!   is below rounding. The lower end is J=0 itself: the quadrature grid uses
!   midpoints, so the first node sits at dJ/2 and the degenerate circular orbit
!   at J=0 is never evaluated.

    select case (dftype)

    case ('gauss')
       Jlo = 0.0d0
       Jhi = 6.0d0*sr

    case ('bimodal')
       Jlo = 0.0d0
       Jhi = bim_Jb + 6.0d0*bim_sb

    case ('spiral')
       Jlo = 0.0d0
       Jhi = spi_J0 + 6.0d0*spi_s0

    case ('king')
!      The profile is exactly zero from the cut on, so the cut is the window.
       Jlo = 0.0d0
       Jhi = kin_Jt

    case default
       Jlo = 0.0d0
       Jhi = 0.0d0
       call df0_unknown

    end select

  end subroutine df0_Jrange


!> Print which distribution is in use, with its constants, so that the run log
!! records the shape the particles were drawn from.

  subroutine df0_report

    implicit none

    print *
    print *, 'Initial distribution: ',trim(dftype)

    select case (dftype)
    case ('gauss')
       print '(a,es10.3,a,es10.3)', '   sigma_Q = ',sp,'   sigma_J = ',sr
    case ('bimodal')
       print '(a,f5.2,a,f5.2)', '   b1 = ',bim_b1,'   b2 = ',bim_b2
       print '(a,f6.3,a,f6.3)', '   Ja = ',bim_Ja,'   sa = ',bim_sa
       print '(a,f6.3,a,f6.3)', '   Jb = ',bim_Jb,'   sb = ',bim_sb
       print '(a,f5.2)',        '   wb = ',bim_wb
    case ('spiral')
       print '(a,f6.3,a,f6.3)', '   J0 = ',spi_J0,'   s0 = ',spi_s0
       print '(a,f7.2,a,f5.2)', '   beta = ',spi_beta,'   sigma_Q = ',spi_sq
    case ('king')
       print '(a,f6.3,a,es10.3)', '   Jt = ',kin_Jt,'   sigma_E^2 = ',kin_sE2
       print '(a,f5.2)',          '   eps = ',kin_eps
    end select

    print *

  end subroutine df0_report


end module distribution
