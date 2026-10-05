! U8. Una particula aislada de masa m a radio x: la masa que ve el campo,
! M = -F(r_N) r_N**2 en el ultimo nodo (teorema de las capas: fuera de la
! capa el campo es -m/r**2).
!  n = 1: debe ser m para todo x. La recta que integra poisson_rk pesa el
!         nodo k con 4 pi dr (r_k**2 + dr**2/6), que es el volumen del peso
!         lineal con que se deposito.
!  n = 2, 3: se informa, no se exige. El volumen del deposito es
!         4 pi dr (r_k**2 + (n+1) dr**2/12) y el campo ve menos masa cerca
!         del origen. Igualarla reescalando los valores nodales (como hace
!         vlasov-poisson_PIC) baja el orden de Poisson de 2 a 1 (U4) y
!         deja un piso en el colapso frio; ver poisson_rk.f90.
program t_masa
  use parameters
  use arrays
  use utils
  implicit none
  integer :: n,k,nfail
  real(8) :: xs(7),m,pi,mcampo(7)
  data xs /0.0d0,0.3d0,1.0d0,2.3d0,5.3d0,20.3d0,35.7d0/
  pi = acos(-1.0d0)
  nfail = 0
  print '(a,7f10.1)', '  x/dr:        ',xs
  do n=1,3
    rmin=0; rmax=4.0d0; dr=0.1d0; autointeraction=.true.; BGtype='null'
    bsplineorder=n; Nrc=1; Npc=1; Nlc=1; drc=1; dpc=1; dlc=1
    call set_grid_size(); call alloc_mem_set0(); call construct_grid()
    m = 1.0d-2
    l_part = 1.0d0; f = m/(8.0d0*pi**2); p_part = 0
    do k=1,size(xs)
      r_part = xs(k)*dr
      call poisson_rk()
      mcampo(k) = -force(Nr)*r(Nr)**2/m
    end do
    if (n == 1) then
      print '(a,es9.2)', '  n=1  max |M_campo - m|/m =',maxval(abs(mcampo-1.0d0))
      if (maxval(abs(mcampo-1.0d0)) > 1.0d-12) nfail = nfail + 1
    else
      print '(a,i1,a,7f10.6)', '  n=',n,'  M_campo/m:',mcampo
    end if
    call deallocate_mem()
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U8 masa vista por el campo (n = 1)'
  else
    print '(a)', 'FALLA U8 masa vista por el campo (n = 1)'
  end if
end program
