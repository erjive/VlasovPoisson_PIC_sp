! Ramas de initial_data.f90 retiradas del código activo (ver legacy/README.md).
! Son estados del código con L fija: llenan solo Nrc*Npc de las Nrc*Npc*Nlc
! partículas y nunca asignan l_part, así que con distribución en L daban
! densidad, energía y h_k idénticamente cero.

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
