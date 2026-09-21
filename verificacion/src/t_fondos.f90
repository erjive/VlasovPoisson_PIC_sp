! U5. Fondos de grav_force.f90 y termino centrifugo: fuerza = -dPhi/dr
! (derivada numerica de 6o orden del propio pot_part) y paridad exacta
! Phi(-r) = Phi(r), F(-r) = -F(r). "sphere" se evalua lejos del quiebre r=1.
program t_fondos
  use parameters
  use arrays
  use utils
  implicit none
  character(12) :: bgs(7) = [character(12) :: 'Isochrone','Central','sphere','iso','isotrun','nfw','burkert']
  integer :: b,k,m,nfail
  real(8) :: rs(8), h, d, emax, epar, eodd, x, pp(-3:3), Fnum, Pref
  data rs /0.05d0,0.3d0,0.7d0,0.95d0,1.05d0,3.0d0,7.0d0,15.0d0/
  nfail = 0
  do b=1,7
    rmin=0; rmax=20; dr=0.1; autointeraction=.false.; BGtype=bgs(b)
    bsplineorder=1; Nrc=1; Npc=1; Nlc=1
    call set_grid_size(); call alloc_mem_set0(); call construct_grid()
    emax=0; epar=0; eodd=0
    do k=1,size(rs)
      x = rs(k); h = 1.0d-3*x
      do m=-3,3
        r_part = x+m*h; l_part = 0.7d0
        call grav_force()
        pp(m) = pot_part(1)
      end do
      Fnum = -( (pp(3)-pp(-3))/60 - 3*(pp(2)-pp(-2))/20 + 3*(pp(1)-pp(-1))/4 )/h
      r_part = x; call grav_force(); d = force_part(1); Pref = pot_part(1)
      emax = max(emax, abs(d-Fnum)/abs(d))
      r_part = -rs(k); call grav_force()
      epar = max(epar, abs(pot_part(1)-Pref))
      eodd = max(eodd, abs(force_part(1)+d))
    end do
    print '(a,a12,a,es9.2,a,es9.2,a,es9.2)', '  ',bgs(b),'  |F+dPhi/dr|/|F|=',emax,'  |Phi(-r)-Phi(r)|=',epar,'  |F(-r)+F(r)|=',eodd
    if (emax > 1.0d-9 .or. epar > 0 .or. eodd > 0) nfail = nfail+1
    call deallocate_mem()
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U5 fondos: fuerza=-dPhi/dr y paridad'
  else
    print '(a)', 'FALLA U5 fondos: fuerza=-dPhi/dr y paridad'
  end if
end program
