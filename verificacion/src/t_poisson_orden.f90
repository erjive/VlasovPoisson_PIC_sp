! U4. Solucion analitica de Poisson y orden de convergencia en la malla.
! rho = rho0 (1-r^2)^3 para r<1 (C2 en el borde), M=1:
!   M(r)   = (315/16)(r^3/3 - 3r^5/5 + 3r^7/7 - r^9/9)
!   Phi(r) = -M(r)/r - (315/16)(1-r^2)^4/8   (r<1),   -1/r   (r>=1)
! Muestreo cuasi-continuo (64 cascaras por dr, masa exacta de cada
! subcascara) para que el error de cuadratura de particulas sea despreciable.
! Se espera orden 2 en pot y force de la malla para n = 1,2,3.
program t_poisson_orden
  use parameters
  use arrays
  use utils
  implicit none
  integer :: n,idr,i,nsub,nfail
  real(8) :: drs(3),eF(3),eP(3),a,b,pi,oF,oP
  data drs /0.02d0,0.01d0,0.005d0/
  pi = acos(-1.0d0)
  nfail = 0
  do n=1,3
    do idr=1,3
      dr = drs(idr); rmin=0; rmax=1.5d0; autointeraction=.true.; BGtype='null'
      bsplineorder=n; Npc=1; Nlc=1; drc=1; dpc=1; dlc=1
      nsub = nint(1.0d0/dr)*64; Nrc = nsub
      call set_grid_size(); call alloc_mem_set0(); call construct_grid()
      do i=1,nsub
        a = dble(i-1)/dble(nsub); b = dble(i)/dble(nsub)
        r_part(i) = 0.5d0*(a+b); p_part(i)=0; l_part(i)=1
        f(i) = (Menc(b)-Menc(a))/(8.0d0*pi**2)
      end do
      call poisson_rk()
      eF(idr) = 0; eP(idr) = 0
      do i=1,Nr
        eF(idr) = max(eF(idr),abs(force(i)+Menc(r(i))/r(i)**2))
        eP(idr) = max(eP(idr),abs(pot(i)-Phi(r(i))))
      end do
      call deallocate_mem()
    end do
    oF = log(eF(2)/eF(3))/log(2.0d0); oP = log(eP(2)/eP(3))/log(2.0d0)
    print '(a,i1,a,3es10.2,a,f5.2,a,3es10.2,a,f5.2)', '  n=',n,'  Linf(F):',eF,'  orden',oF,'   Linf(Phi):',eP,'  orden',oP
    if (oF < 1.9d0 .or. oP < 1.9d0) nfail = nfail+1
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U4 Poisson analitico, orden 2'
  else
    print '(a)', 'FALLA U4 Poisson analitico, orden 2'
  end if
contains
  real(8) function Menc(x)
    real(8) :: x,y
    y = min(x,1.0d0)
    Menc = (315.0d0/16.0d0)*(y**3/3-3*y**5/5+3*y**7/7-y**9/9)
  end function
  real(8) function Phi(x)
    real(8) :: x
    if (x>=1.0d0) then
      Phi = -1.0d0/x
    else
      Phi = -Menc(x)/x - (315.0d0/16.0d0)*(1-x**2)**4/8.0d0
    end if
  end function
end program
