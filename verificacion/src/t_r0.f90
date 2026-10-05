! U7. Particula exactamente en r = 0 con L = 0 (orbita radial por el centro).
! El limite del termino centrifugo es 0 y el isocrono y la esfera son
! regulares en el origen: pot = -1/2 y -3/2, fuerza 0. Detecta el defecto N2
! de AUDITORIA_FISICA_2026-09-21.md (0*0/0 = NaN en L**2 r/r**4).
program t_r0
  use parameters
  use arrays
  use utils
  implicit none
  character(12) :: bgs(2) = [character(12) :: 'Isochrone','sphere']
  real(8) :: pex(2) = [-0.5d0, -1.5d0]
  integer :: b,nfail
  nfail = 0
  do b=1,2
    rmin=0; rmax=4; dr=0.1; autointeraction=.false.; BGtype=bgs(b)
    bsplineorder=1; Nrc=1; Npc=1; Nlc=1
    call set_grid_size(); call alloc_mem_set0(); call construct_grid()
    r_part = 0.0d0; l_part = 0.0d0; p_part = 0.3d0; f = 1
    call grav_force()
    print '(a,a10,a,2es12.4,a,2es12.4)', '  ',bgs(b),'  r=0, L=0: pot, F =',pot_part(1),force_part(1),'   esperado',pex(b),0.0d0
    if (.not. (abs(pot_part(1)-pex(b)) < 1.0d-14 .and. abs(force_part(1)) < 1.0d-14)) nfail = nfail+1
    call deallocate_mem()
  end do
  if (nfail == 0) then
    print '(a)', 'PASA U7 particula radial en r=0'
  else
    print '(a)', 'FALLA U7 particula radial en r=0'
  end if
end program
