! U3. Integrador radial de poisson_rk. Con n=1 y una particula en cada nodo,
! avg_rho(i) es el valor nodal deseado (W(0)=1, W(+-1)=0). El integrador
! supone rho lineal entre nodos (constante en [0,r1]) y dice integrar en
! forma cerrada: M y Phi deben coincidir con la integral exacta de ese
! interpolante hasta el redondeo. Vuelca r, rho, pot, force a nodal_<Nr>.dat;
! la referencia en precision extendida la calcula py/ref_nodal.py.
program t_poisson_nodal
  use parameters
  use arrays
  use utils
  implicit none
  integer :: i
  real(8) :: pi,rho_t,vol,Rsup
  character(32) :: arg
  pi = acos(-1.0d0)
  call get_command_argument(1,arg); read(arg,*) rmax
  call get_command_argument(2,arg); read(arg,*) dr
  Rsup = 0.8d0*rmax
  rmin=0; autointeraction=.true.; BGtype='null'
  bsplineorder=1; Npc=1; Nlc=1; drc=1; dpc=1; dlc=1
  Nrc = int((rmax-rmin)/dr) + 1
  call set_grid_size(); call alloc_mem_set0(); call construct_grid()
  do i=1,Nr
    if (r(i) < Rsup) then
      rho_t = exp(-(r(i)/(0.4d0*Rsup))**2)*(1.0d0+0.5d0*sin(3.0d0*r(i)))
    else
      rho_t = 0.0d0
    end if
    vol = 4.0d0*pi*(r(i)**2*dr + dr**3*2.0d0/12.0d0)
    r_part(i) = r(i); p_part(i) = 0; l_part(i) = 1.0d0
    f(i) = rho_t*vol/(8.0d0*pi**2)
  end do
  call poisson_rk()
  write(arg,'(i0)') Nr
  open(10,file='nodal_'//trim(arg)//'.dat')
  do i=1,Nr
    write(10,'(4es26.17)') r(i),avg_rho(i),pot(i),force(i)
  end do
  close(10)
end program
