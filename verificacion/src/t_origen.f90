! U2. Una particula aislada a +x y a -x (una particula puede estar en r<0
! dentro de un paso, antes de la reflexion de main.f90). Por la simetria
! f(r,p) = f(-r,-p):
!   (a) la masa depositada en la malla es la de la particula;
!   (b) la densidad depositada desde -x es la misma que desde +x;
!   (c) la fuerza de autogravedad en -x es la opuesta a la de +x.
! Se prueba x/dr en 0.3 .. 3.7 para n = 1,2,3. Detecta el defecto D2 de
! AUDITORIA_CIENTIFICA_2026-09-21.md (lista de celdas y ventana de la
! autofuerza calculadas con r y no con |r|).
program t_origen
  use parameters
  use arrays
  use utils
  implicit none
  integer :: n,k,nfail
  real(8) :: xs(9),m,pi,mdep,rhop(1:12),fp,em,er,ef
  data xs /0.3d0,0.5d0,0.8d0,1.0d0,1.2d0,1.6d0,1.9d0,2.4d0,3.7d0/
  pi = acos(-1.0d0)
  nfail = 0
  do n=1,3
    rmin=0; rmax=4.0d0; dr=0.1d0; autointeraction=.true.; BGtype='null'
    bsplineorder=n; Nrc=1; Npc=1; Nlc=1; drc=1; dpc=1; dlc=1
    call set_grid_size(); call alloc_mem_set0(); call construct_grid()
    m = 1.0d-2
    l_part = 1.0d0; f = m/(8.0d0*pi**2); p_part = 0
    do k=1,size(xs)
      r_part = xs(k)*dr
      call poisson_rk()
      rhop = avg_rho(1:12); fp = force_part(1)
      r_part = -xs(k)*dr
      call poisson_rk()
      mdep = sum(avg_rho(1:Nr)*4*pi*(r(1:Nr)**2*dr + dr**3*dble(n+1)/12))
      em = abs(mdep-m)/m; er = maxval(abs(avg_rho(1:12)-rhop))/maxval(abs(rhop)); ef = abs(force_part(1)+fp)
      if (em > 1.0d-12 .or. er > 1.0d-12 .or. ef > 1.0d-12*m/(xs(k)*dr)**2) then
        nfail = nfail + 1
        print '(a,i1,a,f4.1,a,es9.2,a,es9.2,a,es9.2)', '  n=',n,' x/dr=',xs(k), &
          ': |dM|/m=',em,'  |rho(-x)-rho(x)|/rho=',er,'  |F(-x)+F(x)|=',ef
      end if
    end do
    call deallocate_mem()
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U2 simetria en el origen'
  else
    print '(a,i3,a)', 'FALLA U2 simetria en el origen (',nfail,' casos)'
  end if
end program
