! Projects the particle distribution onto two test functions Phi_1, Phi_2
! and records the first five angular Fourier modes of each,
!
!   h_k = Int F(r,p,L) Phi(Q,J,L) e^{-ikQ} 8 pi**2 L dr dp dL,
!
! with (Q,J) the radial angle-action pair of the isochrone of unit mass and
! scale for the particle's own L. Each test function is separable,
!
!   Phi_n(Q,J,L) = A_n(Q) B_n(J) C_n(L),
!   A_n(Q) = exp(-sin(Q/2)**2/sq_n**2),
!   B_n(J) = J**2 exp(-(J-j_n)**2/sj_n**2),
!   C_n(L) = exp(-(L-lt_n)**2/slt_n**2),
!
! with its own parameters (j1, sj1, sq1, lt1, slt1 and the same with 2),
! independent of the initial distribution. Unless given, they take the
! values of the distribution: sq=sp, sj=sr, j=0, lt=l0, slt=sl.
!
! Output: hk1.tl / hk2.tl hold |h_k|, and hk1_complex.tl / hk2_complex.tl
! hold Re and Im separately. The phase tells a static residue from a
! rotating mode, and lets discretisation noise cancel in a complex average
! of several runs.

  subroutine analysish

    use parameters
    use arrays
!$  use omp_lib

    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8),dimension(1:Npart) :: energy,s, s1, s2, er1, er2,argaux
    complex(8) :: ii                          !imaginary unit
    complex(8) :: expv                        !exp(-ii*Qr(j)) for the current particle
    complex(8),dimension(0:4) :: hk1,hk2      !h_k modes of each test function
    real(8),dimension(0:4) :: abs_hk1,abs_hk2
    integer  :: i,j,k                         !Counters
    integer, parameter :: mode = 4            !Highest mode

    real(8) :: eta
    real(8) :: smallpi

!   Simpson intervals for the angular integral (must be even). With sq=0.1
!   A(Q) is narrow: 20 intervals give a_k with relative errors of 1.2e-2
!   (k=0) to 3.2e-2 (k=4), 512 give <= 2.2e-16 against the closed form
!   a_k = exp(-x) I_k(x), x = 1/(2 sq**2). It is evaluated once per run.
    integer, parameter :: nquad = 512
    real(8) :: quadQ(0:nquad),quadW(0:nquad)  !Quadrature nodes/weights on [0,pi]
    real(8) :: hstep
    real(8), save :: ak1(0:mode),ak2(0:mode)  !Angular coefficients, constant for the whole run
    logical, save :: ak_ready = .false.
    real(8) :: w1,w2                          !f * L * B_n(J) * C_n(L) for the current particle
    integer :: nth,tid                        !Number of threads, thread index
    complex(8) :: loc1(0:mode),loc2(0:mode)   !Partial sums of one thread
    complex(8), allocatable :: part1(:,:),part2(:,:)

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

! The angular Fourier coefficients of A_n,
!
!   a_k = (1/pi) Int_0^pi A_n(Q) cos(kQ) dQ      (A_n is even in Q),
!
! do not depend on the particle, nor on time. They are computed once, with
! Simpson's rule, and each particle only contributes B_n(J) C_n(L) exp(-ikQ):
!
!   h_k = 8 pi**2 drc dpc dlc Sum_j f_j L_j B_n(J_j) C_n(L_j) a_k exp(-ikQ_j).

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
        ak1(i) = 0.d0
        ak2(i) = 0.d0
        do k=0,nquad
          ak1(i) = ak1(i) + quadW(k)*exp(-sin(0.5d0*quadQ(k))**2/sq1**2)*cos(dble(i)*quadQ(k))
          ak2(i) = ak2(i) + quadW(k)*exp(-sin(0.5d0*quadQ(k))**2/sq2**2)*cos(dble(i)*quadQ(k))
        end do
        ak1(i) = hstep/(3.d0*smallpi)*ak1(i)
        ak2(i) = hstep/(3.d0*smallpi)*ak2(i)
      end do
      ak_ready = .true.
    end if

    hk1 = (0.d0,0.d0)
    hk2 = (0.d0,0.d0)

! Unbound particles (E >= 0) have no action-angle variables: the formulas
! above give NaN, and a single one would turn every h_k into NaN. They are
! left out, which is the continuous extension of the test function:
! J -> infinity as E -> 0-, and B(J) -> 0.

! Each thread sums its own share of particles into a local accumulator,
! and the shares are added afterwards in thread order. An OpenMP reduction
! would combine them in whatever order the threads finish, which varies
! from run to run; this way h_k is reproducible to the last bit for a given
! number of threads.

    nth = 1
!$  nth = omp_get_max_threads()
    allocate(part1(0:mode,0:nth-1),part2(0:mode,0:nth-1))
    part1 = (0.d0,0.d0)
    part2 = (0.d0,0.d0)

    !$OMP PARALLEL PRIVATE(j,i,w1,w2,expv,tid,loc1,loc2)
    loc1 = (0.d0,0.d0)
    loc2 = (0.d0,0.d0)
    !$OMP DO SCHEDULE(STATIC)
    do j=1,Npart

      if (energy(j) >= 0.d0) cycle

      w1 = f(j)*l_part(j)*Jr(j)**2*exp(-(Jr(j)-j1)**2/sj1**2)*exp(-(l_part(j)-lt1)**2/slt1**2)
      w2 = f(j)*l_part(j)*Jr(j)**2*exp(-(Jr(j)-j2)**2/sj2**2)*exp(-(l_part(j)-lt2)**2/slt2**2)

      loc1(0) = loc1(0) + w1*ak1(0)
      loc2(0) = loc2(0) + w2*ak2(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        loc1(i) = loc1(i) + w1*ak1(i)*expv**i
        loc2(i) = loc2(i) + w2*ak2(i)*expv**i
      end do

    end do
    !$OMP END DO
    tid = 0
!$  tid = omp_get_thread_num()
    part1(:,tid) = loc1
    part2(:,tid) = loc2
    !$OMP END PARALLEL

    do tid=0,nth-1
      hk1 = hk1 + part1(:,tid)
      hk2 = hk2 + part2(:,tid)
    end do
    deallocate(part1,part2)

    hk1 = drc*dpc*dlc*hk1
    hk2 = drc*dpc*dlc*hk2
    abs_hk1 = 8.0*smallpi**2*abs(hk1)
    abs_hk2 = 8.0*smallpi**2*abs(hk2)
    hk1 = 8.0*smallpi**2*hk1
    hk2 = 8.0*smallpi**2*hk2


! *****************
! *** SAVE DATA ***
! *****************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/hk1.tl',form='formatted',status=filestatus)
     open(102,file=trim(directory)//'/hk2.tl',form='formatted',status=filestatus)
     open(103,file=trim(directory)//'/hk1_complex.tl',form='formatted',status=filestatus)
     open(104,file=trim(directory)//'/hk2_complex.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/hk1.tl',form='formatted',status=filestatus,position='append')
     open(102,file=trim(directory)//'/hk2.tl',form='formatted',status=filestatus,position='append')
     open(103,file=trim(directory)//'/hk1_complex.tl',form='formatted',status=filestatus,position='append')
     open(104,file=trim(directory)//'/hk2_complex.tl',form='formatted',status=filestatus,position='append')
  end if

! Columns: t, |h_0| ... |h_4|; and t, Re h_0, Im h_0, ..., Re h_4, Im h_4.

  write(101,"(7ES24.16)") t,abs_hk1(:)
  write(102,"(7ES24.16)") t,abs_hk2(:)
  write(103,"(11ES24.16)") t,(real(hk1(k)),aimag(hk1(k)),k=0,mode)
  write(104,"(11ES24.16)") t,(real(hk2(k)),aimag(hk2(k)),k=0,mode)

  close(101)
  close(102)
  close(103)
  close(104)

  end subroutine analysish
