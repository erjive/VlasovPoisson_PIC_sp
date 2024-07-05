
subroutine density


! ******************************************************
! ***   FIND DENSITY AND CURRENT IN PHYSICAL SPACE   ***
! ******************************************************
!
! This subroutine integrates over momentum space
! to find rho and curr.  These quantities are
! defined as:
!
!                /
! rho  =  1/r**2 | f dp 
!                /
!
!                /
! curr =  1/r**2 | p f dp 
!                /
!

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor,average_rho,mass

  character(20) filename ! Name of outupt file.


  smallpi = acos(-1.0d0)


! Zero Angular Momentum
  !if (Lfix == 0.0d0) then

     !factor = 1.0d0

! Include Angular Momentum
  !else

     factor = 2.0*smallpi*drc*dpc*dlc!0.25D0/smallpi*drc*dpc

  !endif  

  rho = 0.D0
  avg_rho = 0.D0
  curr = 0.D0

  !!$OMP PARALLEL DO SCHEDULE(GUIDED) private (j) collapse(2)
  !$OMP PARALLEL DO SCHEDULE(GUIDED) collapse(2)

  do i=1,Nr
    do j=1,Npart
      if (abs(r(i)-r_part(j))<=(bsplineorder)*drc) then

        rho(i) = rho(i) + f(j)*l_part(j)*Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
        curr(i) = curr(i)+f(j)*l_part(j)*p_part(j)*Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)
      end if

      if (abs(r(i)-r_part(j))<=(bsplineorder)*dr) then

        avg_rho(i) = avg_rho(i) + f(j)/(dr+dr**3/(12.d0*r(i)**2))*l_part(j)*Wn(bsplineorder,(r(i)-r_part(j))/dr)

      end if

    end do
  end do
  !$OMP END PARALLEL DO
  do i=1,ghost
      rho(i-1) = rho(i)
      avg_rho(i-1) = avg_rho(i)
  end do


  rho = factor*m0*rho/r**2
  avg_rho = factor*m0*avg_rho/r**2

  !mass = 0.D0

  !do i=1,Nr-1
  !    mass = mass + 0.5D0*(rho(i)*r(i)**2+rho(i+1)*r(i+1)**2)
  !end do
  !    mass = mass*4.D0*smallpi*dr

  !print *, "Total mass=", mass

  average_rho = 0.D0

! Integrate with the trapezoidal rule. Second order accurate.

!  do i=1,Nr
    do j=1,Npart
!      if (r1<=r(i) .and. r(i)<=r2) then
      if (r_part(j)>=r1 .and. r_part(j)<= r2) then
      
 !       average_rho = average_rho + 1.0D0/(r2**2-r1**2)*0.5D0*(avg_rho(i)*r(i)**2+avg_rho(i+1)*r(i+1)**2)*dr
        average_rho = average_rho + 1.0D0/(r2**2-r1**2)*f(j)*0.25D0/smallpi

      end if
    end do
!  end do

  filename = 'vlasov_rhomix'
  call save0Ddata(directory,filename,t,average_rho)


end subroutine density



subroutine avg_density

! In order to make the Poisson subroutine more efficient, 
! I will separate the subroutine that calculates 
! the average density from the rest of integrals.

! ******************************************************
! ***   FIND DENSITY AND CURRENT IN PHYSICAL SPACE   ***
! ******************************************************
!
! This subroutine integrates over momentum space
! to find only rho which is defined as:
!                /
! rho  =  1/r**2 | f dp 
!                /

! Include modules.

  use parameters
  use arrays
  use functions
  use utils
! Declare variables.

  implicit none

  integer i,j
  real(8) :: smallpi,factor,diff,contribution

  smallpi = acos(-1.0d0)


! Zero Angular Momentum
!  if (Lfix == 0.0d0) then

!     factor = 1.0d0

! Include Angular Momentum
!  else

     factor = 2.0*smallpi*drc*dpc*dlc

!  endif  

  avg_rho = 0.D0

  !!$OMP PARALLEL DO SCHEDULE(GUIDED) private (j) collapse(2)
!!$OMP PARALLEL DO SCHEDULE(GUIDED) collapse(2)

!  do i=1,Nr
!    do j=1,Npart
!      if (abs(r(i)-r_part(j))<=(bsplineorder+1)*drc) then


!        avg_rho(i) = avg_rho(i) + f(j)/drc*Wn(bsplineorder,(r(i)-r_part(j))/drc)


!      end if
!    end do
!  end do
!  !$OMP END PARALLEL DO

!  !$OMP PARALLEL DO 
!$OMP PARALLEL DO SCHEDULE(GUIDED) SHARED(avg_rho, r, r_part, f, bsplineorder, drc) PRIVATE(i, j, diff, contribution)

  do i = 1, Nr
    do j = 1, Npart
        diff = abs(r(i) - r_part(j))
        if (diff <= (bsplineorder + 1) * dr) then
            contribution = f(j) / (r(i)**2*dr+dr**3/12.d0) *l_part(j)* Wn(bsplineorder, diff / dr)
           !$OMP ATOMIC
            avg_rho(i) = avg_rho(i) + contribution
        end if
    end do
  end do
  !$OMP END PARALLEL DO


  do i=1,ghost
      avg_rho(i-1) = avg_rho(i)
  end do

  avg_rho = factor*m0*avg_rho

end subroutine avg_density
