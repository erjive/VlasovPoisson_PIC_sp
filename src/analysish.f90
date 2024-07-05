  subroutine analysish

    use parameters
    use arrays
    !use utils


    implicit none

    real(8),dimension(1:Npart) :: Qr,Jr      !Action angle arrays
    real(8),dimension(1:Npart) :: energy,s, s1, s2, er1, er2,argaux
    complex(8),dimension (1:Npart)  :: exp_vals,conj_phik
    complex(8) :: ii                          !imaginary unit
    complex(8),dimension(0:4) :: hk           !h_k mode
    real(8),dimension(0:4) :: abs_hk       !Magnitude h_k mode
    integer  :: i,j                           !Counters
    integer  :: mode = 4                          !Number of modes

    real(8) :: raux,paux

    real(8) :: eta
    real(8) :: smallpi

    complex(8) :: phik

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

          Qr(i) = eta - sqrt((-2.d0*energy(i))**3)*sqrt(-l_part(i)**2-2.d0*energy(i)-2.d0-0.5D0/energy(i))/(-2.d0*energy(i))*sin(eta)
    end do
    ii = (0.d0,1.d0)
    exp_vals = exp(-ii*Qr)
    !!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,conj_phik,exp_vals)

    do i = 0,mode

      do j=1,Npart
        conj_phik(j) = conjg(phik(Jr(j),l_part(j),l0,mode, sp, sr, sl))
      end do

      if (i==0) then
        hk(i) = drc*dpc*dlc*sum(f*conj_phik*l_part)
      else 
        hk(i) = drc*dpc*dlc*sum(f*exp_vals**i*conj_phik*l_part)
      end if

    end do
    !!$OMP END PARALLEL DO
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

  function phik(J,l,l0,mode, sp, sr, sl)


    implicit none
    
    integer :: n
    integer :: i
    integer :: mode
    complex(8) :: phik
    real(8) :: J,l,l0,sp,sr,sl
    real(8) :: a, b, h,auxsum
    real(8) :: smallpi

    ! Constants

    smallpi = acos(-1.0d0)

    ! Define the integration limits
    a = 0.0
    b = smallpi
    
    ! Number of intervals (must be even)
    n = 20
    
    ! Calculate the step size
    h = (b - a) / real(n)
    
    ! Perform the integration
    auxsum = phi(J,a,l,l0,mode,sp,sr,sl) + phi(J,b,l,l0,mode,sp,sr,sl)

    do i = 1, n-1, 2
        auxsum = auxsum + 4.0d0 * phi(J,a + real(i) * h,l,l0,mode,sp,sr,sl)
    end do
    do i = 2, n-2, 2
        auxsum = auxsum + 2.0d0 * phi(J,a + real(i) * h,l,l0,mode,sp,sr,sl)
    end do
    
    ! Calculate the result

    phik = h / 3.0d0* (real(auxsum)*2.0d0)

    phik = 0.5d0/smallpi*phik
    contains
    
    ! Define the function to be integrated
    function phi(J, Q, l,l0,mode, sp, sr, sl)
        integer mode
        real(8) :: J,Q,l,l0,sp,sr,sl
        complex(8) :: phi
        complex(8) :: ii

        ii = (0.d0,1.0d0)

        phi = (exp(-sin(0.5d0*Q)**2/sp**2)*exp(-J**2/sr**2)*J**2*exp(-(l-l0)**2/sl**2)*exp(-ii*mode*Q))

       
    end function phi
    

  end function phik
