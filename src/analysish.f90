  subroutine analysish

    use parameters
    use arrays

    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8),dimension(1:Npart) :: energy,s, s1, s2, er1, er2,argaux
    complex(8) :: ii                          !imaginary unit
    complex(8) :: expv                        !exp(-ii*Qr(j)) for the current particle
    complex(8),dimension(0:4) :: hk           !h_k mode
    real(8),dimension(0:4) :: abs_hk       !Magnitude h_k mode
    integer  :: i,j,k                         !Counters
    integer  :: mode = 4                          !Number of modes

    real(8) :: eta
    real(8) :: smallpi

!   Simpson intervals for the angular integral (must be even). With sp=0.1
!   A(Q) is narrow and 20 intervals gave a_k with relative errors of
!   1.2e-2 (k=0) to 3.2e-2 (k=4); 512 give <= 2.2e-16 against the closed
!   form a_k = exp(-x) I_k(x), x = 1/(2 sp**2). It is evaluated once per run.
    integer, parameter :: nquad = 512
    real(8) :: quadQ(0:nquad),quadW(0:nquad)  !Quadrature nodes/weights on [0,pi]
    real(8) :: hstep
    real(8), save :: ak(0:4)                  !Angular coefficients a_k, constant for the whole run
    logical, save :: ak_ready = .false.
    real(8) :: w                              !f * L * B(J) * C(L) for the current particle

    character(20) filestatus


    ! Constants
    smallpi =  acos(-1.0d0)

    energy = -1.0/(1.0D0+dsqrt(1.0D0+r_part**2)) + 0.5d0*l_part**2/(r_part**2) + 0.5D0*p_part**2
    er1 = dsqrt((1.d0+energy*(2.d0+l_part**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+l_part**2)))/(2.d0*energy**2))
    er2 = dsqrt((1.d0+energy*(2.d0+l_part**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+l_part**2)))/(2.d0*energy**2))
    s1 = 1.d0 + sqrt(1.d0+er1**2)
    s2 = 1.d0 + sqrt(1.d0+er2**2)
    s  = 1.d0 + sqrt(1.d0+r_part**2)
    argaux = (s1+s2-2.0*s)/(s2-s1)
    Jr = 1.d0/sqrt(-2.d0*energy)-0.5d0*(l_part+sqrt(l_part**2+4.d0))

    do i=1,Npart

          if (p_part(i)>=0.d0) then
            eta = dacos(sign(min(abs(argaux(i)),1.0),argaux(i)))

          else
            eta = dacos(-sign(min(abs(argaux(i)),1.0),argaux(i)))+smallpi
          end if

          Qr(i) = eta - sqrt((-2.d0*energy(i))**3) &
                  *sqrt(-l_part(i)**2-2.d0*energy(i)-2.d0-0.5D0/energy(i))/(-2.d0*energy(i))*sin(eta)
    end do
    ii = (0.d0,1.d0)

! The test function is separable, Phi(Q,J,L) = A(Q) B(J) C(L), with
!
!   A(Q) = exp(-sin(Q/2)**2/sp**2),  B(J) = J**2 exp(-J**2/sr**2),
!   C(L) = exp(-(L-l0)**2/sl**2),
!
! so its angular Fourier coefficients
!
!   a_k = (1/pi) Int_0^pi A(Q) cos(kQ) dQ      (A is even in Q)
!
! do not depend on the particle, nor on time. They are computed once, with
! Simpson's rule, and each particle only contributes B(J) C(L) exp(-ikQ):
!
!   h_k = 8 pi**2 drc dpc dlc Sum_j f_j L_j B(J_j) C(L_j) a_k exp(-ikQ_j).

    if (.not. ak_ready) then
      hstep = smallpi/dble(nquad)
      do k=0,nquad
        quadQ(k) = dble(k)*hstep
      end do
      quadW(0)     = 1.d0
      quadW(nquad) = 1.d0
      do k=1,nquad-1,2
        quadW(k) = 4.d0
      end do
      do k=2,nquad-2,2
        quadW(k) = 2.d0
      end do
      do i=0,mode
        ak(i) = 0.d0
        do k=0,nquad
          ak(i) = ak(i) + quadW(k)*exp(-sin(0.5d0*quadQ(k))**2/sp**2)*cos(dble(i)*quadQ(k))
        end do
        ak(i) = hstep/(3.d0*smallpi)*ak(i)
      end do
      ak_ready = .true.
    end if

    hk = (0.d0,0.d0)

! Unbound particles (E >= 0) have no action-angle variables: the formulas
! above give NaN, and a single one turned every h_k into NaN for the whole
! run. They are left out, which is the continuous extension of the test
! function: J -> infinity as E -> 0-, and B(J) = J**2 exp(-J**2/sr**2) -> 0.

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,i,w,expv) REDUCTION(+:hk)
    do j=1,Npart

      if (energy(j) >= 0.d0) cycle

      w = f(j)*l_part(j)*Jr(j)**2*exp(-Jr(j)**2/sr**2)*exp(-(l_part(j)-l0)**2/sl**2)

      hk(0) = hk(0) + w*ak(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        hk(i) = hk(i) + w*ak(i)*expv**i
      end do

    end do
    !$OMP END PARALLEL DO

    hk = drc*dpc*dlc*hk
    abs_hk = 8.0*smallpi**2*abs(hk)


! *****************
! *** SAVE DATA ***
! *****************

! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim("hk")//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim("hk")//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(7ES16.8)") t,abs_hk(:)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


  end subroutine analysish
