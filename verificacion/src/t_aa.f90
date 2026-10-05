! U6. Mapa angulo-accion del isocrono (utils.f90).
!  (a) ida y vuelta (Q,J,L) -> invert_QJ_to_rp -> init_action_angle -> (Q',J')
!      y energia H(r,p,L) = -1/(2(J+c)^2), para L en {0, 0.01, 0.5, 1, 2, 5},
!      J en [1e-4, 10] y Q en [0, 2pi);
!  (b) kepler_eta, que resuelve la ecuacion de Kepler para invert_QJ_to_rp,
!      converge para toda excentricidad e < 1 y todo Q.
! Detecta los defectos D3 (Newton diverge para e >~ 0.98) y D4 (L=0:
! radicando negativo por redondeo, NaN y r=0) de la auditoria. Hasta
! corregir D3 la parte (b) probaba una copia del Newton sin salvaguarda,
! que era lo que tenia invert_QJ_to_rp: fallaba en 190 de 10000 casos.
program t_aa
  use parameters
  use arrays
  use utils
  implicit none
  integer :: i,iq,ij,il,n,nbad,nfailQ,ntot
  real(8) :: Ls(6) = [0.0d0,0.01d0,0.5d0,1.0d0,2.0d0,5.0d0]
  real(8) :: Q,J,L,rv,pv,E,c,eQ,eJ,eE,pi,ecc,eta,g,Qm,resmax,dQ
  pi = acos(-1.0d0)
  rmin=0; rmax=20; dr=0.1; Nrc=40; Npc=40; Nlc=6; autointeraction=.false.
  call set_grid_size(); call alloc_mem_set0(); call construct_grid()
  n = 0
  do il=1,6
    do ij=1,40
      do iq=1,40
        n = n+1
        L = Ls(il); J = 1.0d-4*10.0d0**(dble(ij-1)*5.0d0/39.0d0); Q = (dble(iq)-0.5d0)*2*pi/40
        call invert_QJ_to_rp(Q,J,L,rv,pv)
        r_part(n)=rv; p_part(n)=pv; l_part(n)=L; f(n)=1
      end do
    end do
  end do
  call init_action_angle()
  n=0; ntot=0
  do il=1,6
    eQ=0; eJ=0; eE=0; nfailQ=0
    do ij=1,40
      do iq=1,40
        n=n+1
        L = Ls(il); J = 1.0d-4*10.0d0**(dble(ij-1)*5.0d0/39.0d0); Q=(dble(iq)-0.5d0)*2*pi/40
        c = 0.5d0*(L+sqrt(L*L+4))
        E = -1.0d0/(1+sqrt(1+r_part(n)**2)) + 0.5d0*L*L/max(r_part(n)**2,1d-300) + 0.5d0*p_part(n)**2
        eE = max(eE, abs(E + 0.5d0/(J+c)**2)/(0.5d0/(J+c)**2))
        eJ = max(eJ, abs(j0_part(n)-J)/J)
        dQ = abs(modulo(q0_part(n)-Q+pi,2*pi)-pi)
        eQ = max(eQ, dQ)
        if (.not. dQ < 1.0d-8) nfailQ = nfailQ+1
      end do
    end do
    print '(a,f5.2,a,es9.2,a,es9.2,a,es9.2,a,i5,a)','  L=',Ls(il),'  max|dQ|=',eQ,'  max|dJ|/J=',eJ, &
          '  max|dE/E|=',eE,'  casos con |dQ|>1e-8:',nfailQ,'/1600'
    ntot = ntot + nfailQ
  end do
  nbad=0; resmax=0
  do i=0,49
    ecc = 1.0d0 - 10.0d0**(-dble(i)/49.0d0*8.0d0)
    do iq=0,199
      Qm = dble(iq)/200.0d0*2*pi
      eta = kepler_eta(Qm,ecc,1.0d-14)
      g = eta - ecc*sin(eta) - Qm
      resmax = max(resmax,abs(g))
      if (.not. abs(g) <= 1.0d-12) nbad = nbad+1
    end do
  end do
  print '(a,i5,a,es9.2)','  Kepler (kepler_eta de utils.f90), e en [0, 1-1e-8]: sin converger',nbad,' de 10000; max|g|=',resmax
  if (ntot == 0 .and. nbad == 0) then
    print '(a)', 'PASA U6 mapa angulo-accion'
  else
    print '(a)', 'FALLA U6 mapa angulo-accion'
  end if
end program
