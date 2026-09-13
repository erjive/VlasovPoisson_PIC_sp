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

    integer, parameter :: nquad = 20          !Simpson intervals for the phik integral (must be even)
    real(8) :: quadQ(0:nquad),quadW(0:nquad)  !Quadrature nodes/weights on [0,pi], built once
    real(8) :: hstep,quadnorm
    real(8) :: gval(0:nquad)                  !mode-independent integrand, per quadrature node
    real(8) :: contrib(0:4)                   !phik(.,mode) for mode=0..4, current particle

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

! Quadrature nodes/weights for Simpson's rule on [0,pi] with nquad
! (even) sub-intervals -- the same grid for every particle and every
! mode, so build it once instead of inside phik() on every one of the
! (mode+1)*Npart calls it used to get (this loop's own cost is
! negligible, O(nquad), done once per analysish() call).

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
    quadnorm = hstep/(3.d0*smallpi)

    hk = (0.d0,0.d0)

! phik(J,l,mode) = (1/pi) * Simpson[ g(J,l,Q)*cos(mode*Q) dQ, Q=0..pi ],
! with g(J,l,Q) = exp(-sin(Q/2)^2/sp^2)*exp(-J^2/sr^2)*J^2*exp(-(l-l0)^2/sl^2)
! the part of the original phi(J,Q,l,l0,mode,sp,sr,sl) integrand that
! does NOT depend on mode.  The previous implementation called phik()
! once per (particle, mode) pair -- (mode+1)=5 times per particle --
! and each call recomputed g(J,l,Q) at all 21 quadrature points from
! scratch, i.e. 5x more exp() evaluations than necessary, since g
! only depends on the particle (through J=Jr(j), l=l_part(j)), not on
! which of the 5 modes is being accumulated.  Computing gval(:) once
! per particle and reusing it for all 5 modes removes that
! redundancy.  It also replaces the original complex phi(Q) =
! g(Q)*exp(-i*mode*Q) with g(Q)*cos(mode*Q) directly: phik's own
! result was always real anyway (the previous code built it from
! "real(auxsum)*2", silently discarding the imaginary part of the
! Simpson sum every time), so the sin(mode*Q) part it implicitly threw
! away is simply never computed now.
!
! Each particle's contribution to hk(0:4) is independent and only
! summed, so this parallelizes over particles with a plain reduction
! on the (tiny, 5-element) hk accumulator -- no atomics needed.

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,k,i,gval,contrib,expv) REDUCTION(+:hk)
    do j=1,Npart

      do k=0,nquad
        gval(k) = exp(-sin(0.5d0*quadQ(k))**2/sp**2)*exp(-Jr(j)**2/sr**2)*Jr(j)**2*exp(-(l_part(j)-l0)**2/sl**2)
      end do

      do i=0,mode
        contrib(i) = 0.d0
        do k=0,nquad
          contrib(i) = contrib(i) + quadW(k)*gval(k)*cos(dble(i)*quadQ(k))
        end do
        contrib(i) = quadnorm*contrib(i)
      end do

      hk(0) = hk(0) + f(j)*l_part(j)*contrib(0)

      expv = exp(-ii*Qr(j))
      do i=1,mode
        hk(i) = hk(i) + f(j)*l_part(j)*contrib(i)*expv**i
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
