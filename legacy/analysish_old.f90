  subroutine analysish

    use parameters
    use arrays
    !use utils


    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8) :: energy,s, s1, s2, er1, er2,argaux
    complex(8),dimension (1:Npart)  :: exp_vals,conj_phik
    complex(8) :: ii                          !imaginary unit
    complex(8),dimension(0:4) :: hk1,hk2           !h_k mode
    real(8),dimension(0:4) :: abs_hk1,abs_hk2       !Magnitude h_k mode
    integer  :: i,j                           !Counters
    integer  :: mode = 4                          !Number of modes

    real(8) :: eta
    real(8) :: smallpi

    complex(8) :: phik

    character(20) filestatus


    ! Constants
    smallpi =  acos(-1.0d0)


    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,energy,s, s1, s2, er1, er2,argaux) SHARED(Qr,Jr)
    do i = 1,Npart
      energy = -1.0/(1.0D0+dsqrt(1.0D0+r_part(i)**2)) + 0.5d0*Lfix**2/(r_part(i)**2) + 0.5D0*p_part(i)**2
      er1 = dsqrt((1.d0+energy*(2.d0+Lfix**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
      er2 = dsqrt((1.d0+energy*(2.d0+Lfix**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+Lfix**2)))/(2.d0*energy**2))
      s1 = 1.d0 + dsqrt(1.d0+er1**2)
      s2 = 1.d0 + dsqrt(1.d0+er2**2)
      s  = 1.d0 + dsqrt(1.d0+r_part(i)**2)
      argaux = (s1+s2-2.0*s)/(s2-s1)
      Jr(i) = 1.d0/dsqrt(-2.d0*energy)-0.5d0*(Lfix+dsqrt(Lfix**2+4.d0))

  !  do i=1,Npart

          if (p_part(i)>=0.d0) then
            eta = dacos(sign(min(abs(argaux),1.0D0),argaux))

          else
            eta = dacos(-sign(min(abs(argaux),1.0D0),argaux))+smallpi
          end if

          Qr(i) = eta - dsqrt((-2.d0*energy)**3)*sqrt(-Lfix**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*dsin(eta)
      
    end do
    !$OMP END PARALLEL DO
   
 ii = (0.d0,1.d0)
    !exp_vals = exp(-ii*Qr)
    !!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,conj_phik,exp_vals)

! 
    do i = 0,mode

      do j=1,Npart
        conj_phik(j) = conjg(phik(Jr(j), i, j1,sq1, sj1))
      end do


      if (i==0) then
        hk1(i) = drc*dpc*sum(f*conj_phik)
      else 
        exp_vals = exp(-ii*Qr*i)
        hk1(i) = drc*dpc*sum(f*exp_vals*conj_phik)
      end if

    end do
    !!$OMP END PARALLEL DO
    abs_hk1 = abs(hk1)

    do i = 0,mode

      do j=1,Npart
        conj_phik(j) = conjg(phik(Jr(j), i, j2,sq2, sj2))
      end do


      if (i==0) then
        hk2(i) = drc*dpc*sum(f*conj_phik)
      else 
        exp_vals = exp(-ii*Qr*i)
        hk2(i) = drc*dpc*sum(f*exp_vals*conj_phik)
      end if

    end do
    !!$OMP END PARALLEL DO
    abs_hk2 = abs(hk2)


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
     open(101,file=trim(directory)//'/'//trim("hk1")//'.tl',form='formatted',status=filestatus)
     open(102,file=trim(directory)//'/'//trim("hk2")//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim("hk1")//'.tl',form='formatted',status=filestatus,position='append')
     open(102,file=trim(directory)//'/'//trim("hk2")//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(7ES16.8)") t,abs_hk1(:)
  write(102,"(7ES16.8)") t,abs_hk2(:)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)
  close(102)


  end subroutine analysish

  function phik(J, l, J0,sq, sj)


    implicit none
    
    integer :: n
    integer :: i
    integer :: l
    complex(8) :: phik
    real(8) :: J,J0,sq,sj
    real(8) :: a, b, h,auxsum
    real(8) :: smallpi

    ! Constants

    smallpi = acos(-1.0d0)

    ! Define the integration limits
    a = 0.0
    b = smallpi
    
    ! Number of intervals (must be even)
    n = 512
    
    ! Calculate the step size
    h = dble((b - a) / real(n))
    
    ! Perform the integration
    auxsum = phi(J,a,l,J0,sq,sj) + phi(J,b,l,J0,sq,sj)

    do i = 1, n-1, 2
        auxsum = auxsum + 4.0d0 * phi(J,a + dble(i) * h,l,J0,sq,sj)
    end do
    do i = 2, n-2, 2
        auxsum = auxsum + 2.0d0 * phi(J,a + dble(i) * h,l,J0,sq,sj)
    end do
    
    ! Calculate the result

    phik = dble(h / 3.0d0* (real(auxsum)*2.0d0))

    phik = 0.5d0/smallpi*phik

    contains
    
    ! Define the function to be integrated
    function phi(J, Q, l, J0, sq, sj)
        integer l
        real(8) :: J,Q,J0,sq,sj
        complex(8) :: phi
        complex(8) :: ii

        ii = (0.d0,1.0d0)

        phi = (exp(-sin(0.5d0*Q)**2/sq**2)*exp(-(J-J0)**2/sj**2)*J**2*exp(-ii*l*Q))

       
    end function phi
    

  end function phik
