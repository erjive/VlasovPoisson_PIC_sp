! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays
    use utils

    implicit none

    logical :: accepted
    integer :: i,j,k,indx
    real(8) :: smallpi,f_max
    real(8) :: raux,paux,laux
    real(8) :: rand3(3)
    real(8) :: gaussian
    real(8) :: w,x,y,z
    real(8) :: energy 

!   Auxiliary variables for a distribution function 
!   that depends on action-angle varialbes
    real(8) :: J3, Q3, w3, s, s1, s2, er1, er2, eta, argaux

    smallpi = acos(-1.0d0)

! 
! For a fixed value of L, We generate particles for an 
! arbitrary distribution function f(r,pr,L) via an acceptance-rejection method.
! Let fmax the maximum value of f. We generate arbitrary (x,y,z) numbers 
! in the range of (rmin,rmax), (pmin,pmax), (0,fmax) respectively. 
! Then evaluate W=f(x,y,L), if z<=W, accept the point, 
! otherwise, repeat until the condition in fulfilled.


! Initial data for the density function. Notice that
! we add two copies of the function in order to guarantee 
! that the ! boundary condition f(-r,-p) = f(r,p) is satisfied.

! Find the size of the cell

    drc = (rmaxc-rminc)/dble(Nrc)
    dpc = (pmaxc-pminc)/dble(Npc)
    dlc = (lmaxc-lminc)/dble(Nlc)

    print *, "(drc,dpc,dlc)=",drc,dpc,dlc
    if(state.eq."gaussian1") then


!      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      indx = 1 
      do k=1,Nlc
      do i=1,Nrc
        do j=1,Npc
          
            raux = rminc+(dble(i))*drc
            paux = pminc+dble(j)*dpc
            laux = lminc+dble(k)*dlc
            r_part(indx)  = raux
            p_part(indx)  = paux
            l_part(indx)  = laux
            f(indx)       = gaussian(1.0D0,r0,p0,l0,raux,paux,laux,sr,sp,sl)

            indx = indx+1
          end do
        end do
      end do
!      !$OMP END PARALLEL DO

      f_max = maxval(f)

      indx = 1
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc
        
             if (f(indx)<= cutoff*f_max) then
                f(indx)=0.0D0
                r_part(indx) = 1000000.D0
             end if
!          !print *, Qr,Jr
             indx = indx+1
          end do
        end do 
      end do
 
      call reduce_arrays


      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc
    
    else if (state.eq."Plummer") then 

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux,energy)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i))*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux

          energy = -1.0/(1.0D0+sqrt(1.0D0+raux**2)) + 0.5d0*l_part(i)**2/(raux**2 + eps*eps) + 0.5D0*paux**2



          if (energy<0.d0) then
            f((i-1)*Npc+j) = (-energy)**3.5
          else
            f((i-1)*Npc+j) = 0.0D0
          end if
        end do
      end do
      !$OMP END PARALLEL DO
    
    else if(state == "compact") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = a0/(32.D0*smallpi**2*sr*sp)*(1.D0+dcos(smallpi/sr*(raux-r0))) * &
                                                              (1.D0+dcos(smallpi/sp*(paux-p0))) 
        end do
      end do
      !$OMP END PARALLEL DO

      f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)

    else if(state == "compact2") then

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      do i=1,Nrc
        do j=1,Npc
          raux = rminc+(dble(i)-0.5D0)*drc
          paux = pminc+dble(j)*dpc

          r_part((i-1)*Npc+j) = raux
          p_part((i-1)*Npc+j) = paux          
          f((i-1)*Npc+j)      = 2.0D0*a0/(9.D0*smallpi**2*sr*sp)*(dcos(0.5D0*smallpi/sr*(raux-r0)))**4* &
                                                            (dcos(0.5D0*smallpi/sp*(paux-p0)))**4 
        end do
      end do
      !$OMP END PARALLEL DO

      f = f*drc*dpc*8.D0*smallpi**2
      print *, "Initial total mass=",sum(f)

    else if(state .eq."aa") then

      !!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(j,raux,paux)
      indx = 1
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc
          
            raux = rminc+dble(i)*drc
            paux = pminc+dble(j)*dpc
            laux = lminc+dble(k)*dlc

            r_part(indx) = raux
            p_part(indx) = paux     
            l_part(indx) = laux

            energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*laux**2/(raux**2) + 0.5D0*paux**2
            er1 = dsqrt((1.d0+energy*(2.d0+laux**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+laux**2)))/(2.d0*energy**2))
            er2 = dsqrt((1.d0+energy*(2.d0+laux**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+laux**2)))/(2.d0*energy**2))
            s1 = 1.d0 + sqrt(1.d0+er1**2)
            s2 = 1.d0 + sqrt(1.d0+er2**2)
            s  = 1.d0 + sqrt(1.d0+raux**2)
            argaux = (s1+s2-2.0*s)/(s2-s1)

            if (paux>=0.d0) then
              eta = dacos(sign(min(abs(argaux),1.0),argaux))

            else
              eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
            end if

            Q3 = eta - sqrt((-2.d0*energy)**3)*sqrt(-laux**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*sin(eta)

            J3 = 1.d0/sqrt(-2.d0*energy)-0.5d0*(laux+sqrt(laux**2+4.d0))

            !w3 = 8.0D0/(laux+2.0D0*J3+sqrt(laux**2+4.0D0))**3

            Q3 = Q3

            f(indx) = exp(-sin(0.5d0*Q3)**2/sp**2)*exp(-J3**2/sr**2)*J3**2*exp(-(laux-l0)**2/sl**2)

            if (f(indx) /= f(indx) ) then
              f(indx) = 0.D0
              r_part(indx) = 10000.D0
            end if
            indx = indx+1
          end do
        end do
      end do
!      !$OMP END PARALLEL DO
      
      f_max = maxval(f)

      indx = 1
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc
        
             if (f(indx)<= cutoff*f_max) then
                f(indx)=0.0D0
                r_part(indx) = 1000000.D0
             end if
!          !print *, Qr,Jr
             indx = indx+1
          end do
        end do 
      end do
 
      call reduce_arrays


      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc


    else if(state.eq."checkpoint") then !NO IMPLEMENTED

!       open(101,file=CheckPointFile)

!       do i=0,Nr
!          do j=0,Np
!             read(101,*) aux1, aux2, f(i,j)
!          end do
!       end do

!       close(101)

    else if(state.eq."other3") then !NO IMPLEMENTED
       f = 0.0d0       
    endif

  end subroutine initial_data

  function gaussian(a0,r0,p0,l0,x,y,z,sr,sp,sl)

  implicit none

  real(8) a0
  real(8) gaussian
  real(8) r0,p0,l0,x,y,z,sr,sp,sl
  real(8) smallpi

  smallpi = acos(-1.0d0)

  gaussian = dble(a0)*&
             (dexp(-(x-r0)**2/sr**2)*dexp(-(y-p0)**2/sp**2)+ &
              dexp(-(x+r0)**2/sr**2)*dexp(-(y+p0)**2/sp**2))*&
              dexp(-(z-l0)**2/sl**2)

  end function gaussian


