! U1. Funciones de forma W_n (functions.f90): particion de la unidad, area 1,
! segundo momento (n+1)/12 y simetria. El volumen V_i del deposito
! (density.f90) se basa en estas dos ultimas.
program t_forma
  use functions
  implicit none
  integer :: n,i,k,nfail
  real(8) :: x,s,area,m2,errpu,errsym,h
  nfail = 0
  do n=1,3
    errpu = 0; errsym = 0
    do i=0,1000
      x = dble(i)/1000.0d0 - 0.5d0
      s = 0
      do k=-5,5
        s = s + Wn(n,x-dble(k))
      end do
      errpu = max(errpu,abs(s-1.0d0))
      errsym = max(errsym,abs(Wn(n,x)-Wn(n,-x)))
    end do
    h = 1.0d-4; area = 0; m2 = 0
    do i=-40000,40000
      x = dble(i)*h
      area = area + Wn(n,x)*h
      m2 = m2 + x*x*Wn(n,x)*h
    end do
    if (errpu > 1.0d-14 .or. errsym > 0.0d0 .or. abs(area-1) > 1.0d-12 .or. &
        abs(m2-dble(n+1)/12) > 1.0d-8) nfail = nfail + 1
    print '(a,i1,a,es9.2,a,es9.2,a,es9.2,a,es9.2)', '  n=',n,'  |sum W-1|=',errpu, &
      '  |W(x)-W(-x)|=',errsym,'  |area-1|=',abs(area-1),'  |m2-(n+1)/12|=',abs(m2-dble(n+1)/12)
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U1 funciones de forma'
  else
    print '(a)', 'FALLA U1 funciones de forma'
  end if
end program
